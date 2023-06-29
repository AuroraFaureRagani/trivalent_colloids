#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 2
#define N 500
#define NBIN 1000
#define MAX_RINGSIZE 7

/* Initialization variables */
const int mc_steps = 300000;
const int output_steps = 100;
const int n_patches = 3;
const double diameter = 1.0;
const double cos_max = 0.99;
const double delta = 0.1;
const double delta_theta = 0.02;
const double u0 = 1.0;

const double density = 0.1;
const double temp = 0.15; // temp 1/beta

/* Simulation variables */
int n_particles = 100;
double e_cut;
double radius;
double particle_volume;
double r[N][NDIM];
double theta[N];
double box[NDIM];
int nhis[NBIN];
double g[NBIN];
int nn[N];
double energy;
double avg_nn = 0.0;
double r_cut  = 2.5;
int bonds[N][N];
int ring_number[4];

/* Functions */

double distance(double vector1[NDIM], double vector2[NDIM], double box_size[NDIM]) {
    // calculate the coordinates differences for each direction
    double dx = fabs(vector2[0] - vector1[0]);
    double dy = fabs(vector2[1] - vector1[1]);
    
    // apply periodic boundary conditions
    if (dx > box_size[0]/2) { dx = box_size[0] - dx; }
    if (dy > box_size[1]/2) { dy = box_size[1] - dy; }
    //calculate the euclidean distance
    double dist = sqrt(dx*dx + dy*dy);
    return dist;
}

void read_data(char* init_filename){

   /* open the file in 'reading-mode' */
   FILE *read_coord;
   read_coord = fopen(init_filename, "r");

   /* check if the file was successfully opened */
   if (read_coord == NULL)
	{
		fprintf(stderr,"Error with file opening: %s!\n",init_filename);
		exit(EXIT_FAILURE);
	}

   /* read the number of particles */
   fscanf(read_coord, "%d\n", &n_particles);
   n_particles = n_particles/4;
   /* read the size of the box*/
    double min, max;
    for(int i = 0; i<NDIM; i++){
     fscanf(read_coord, "%lf %lf\n", &min, &max);
     box[i] = max-min;
    }
    double a, b;
    fscanf(read_coord, "%lf %lf\n", &a, &b);
   
   /* read the coordinates of the spheres */
    for(int i = 0; i < n_particles; i++){
        for(int j = 0; j < NDIM; j++){
            fscanf(read_coord, "%lf ", &(r[i][j]));
        }
        fscanf(read_coord, "%*lf %*lf %*d %lf\n", &(theta[i]));
        fscanf(read_coord, "%*lf %*lf %*lf %*lf %*d\n");
        fscanf(read_coord, "%*lf %*lf %*lf %*lf %*d\n");
        fscanf(read_coord, "%*lf %*lf %*lf %*lf %*d\n");

    }

   /* close the file */
   fclose(read_coord);

}


void write_bonds(int step){
    char buffer[128];
    sprintf(buffer, "bonds/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    char buffer2[128];
    sprintf(buffer2, "angles/coords_step%07d.dat", step);
    FILE* fp2 = fopen(buffer2, "w");

    if (buffer == NULL || buffer2 == NULL)
	{
		fprintf(stderr,"Error with file opening!\n");
		exit(EXIT_FAILURE);
	}
    int d, n;
    fprintf(fp, "%lf\n", fabs(2*energy) + fabs(2*energy)*n_patches);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n", 0.0, box[d]);
    }
    fprintf(fp, "%lf %lf\n", 0.0, 0.0);
    for(int i = 0; i < n_particles; ++i){
        for(int j = i+1; j < n_particles; ++j){
            if (bonds[i][j] == 1){
                /* se c'è legame, scrivi coordinate */
                for(d = 0; d < NDIM; ++d) {
                    fprintf(fp, "%f\t", r[i][d]);
                    fprintf(fp2, "%f\t", r[i][d]);
                }
                fprintf(fp2, "%f\t%d\n", theta[i], i);
                fprintf(fp, "%lf\t%lf\t%d\n", 0.0, diameter, 1);
                for(int p=0; p<n_patches; p++){
                    double patchX = r[i][0] + radius*cos(theta[i]+ p*2*M_PI/n_patches);
                    double patchY = r[i][1] - radius*sin(theta[i] + p*2*M_PI/n_patches);
                    fprintf(fp, "%f\t", patchX);
                    fprintf(fp, "%f\t", patchY);
                    fprintf(fp, "%f\t", 0.0);
                    fprintf(fp, "%lf\t%d\n", diameter/2.5, 2);
                }

                for(d = 0; d < NDIM; ++d) {
                    fprintf(fp, "%f\t", r[j][d]);
                    fprintf(fp2, "%f\t", r[j][d]);
                }
                fprintf(fp2, "%f\t%d\n", theta[j], j);
                fprintf(fp, "%lf\t%lf\t%d\n", 0.0, diameter, 1);
                for(int p=0; p<n_patches; p++){
                    double patchX = r[j][0] + radius*cos(theta[j]+ p*2*M_PI/n_patches);
                    double patchY = r[j][1] - radius*sin(theta[j] + p*2*M_PI/n_patches);
                    fprintf(fp, "%f\t", patchX);
                    fprintf(fp, "%f\t", patchY);
                    fprintf(fp, "%f\t", 0.0);
                    fprintf(fp, "%lf\t%d\n", diameter/2.5, 2);
                }
            }
        
        }
        
    }
    fclose(fp);
}

double potential(){
    double V = 0.0; 
    double dx, dy, patch_i_x, patch_i_y, patch_j_x, patch_j_y, px, py;
    double rx, ry, uix, uiy, ujx, ujy;
    double dist, norm;
    int boundary_cond_x = 0, boundary_cond_y = 0;

    for(int i = 0; i < n_particles; i++){
        for(int j = 0; j < n_particles; j++){
            bonds[i][j] = 0;
        }
    }
    for(int i = 0; i < n_particles; i++){
        for(int j = i+1; j < n_particles; j++){
            dist = distance(r[i], r[j], box);
            if(dist < diameter){ 
                printf("ERROR: check for overlaps did not succeed\n"); 
                }
            if(dist >= diameter && dist <= r_cut){ 
                /* check if 2 patches intersect r_ij */
                dx = r[j][0]-r[i][0];
                if ( (int)(2.0 * dx / box[0]) * box[0] != 0){
                    dx -= (int)(2.0 * dx / box[0]) * box[0];
                    boundary_cond_x = 1;
                }
                
                dy = r[j][1]-r[i][1];
                if ( (int)(2.0 * dy / box[1]) * box[1] != 0){
                    dy -= (int)(2.0 * dy / box[1]) * box[1];
                    boundary_cond_y = 1;
                }

                rx = dx/sqrt(dx*dx + dy*dy);
                ry = dy/sqrt(dx*dx + dy*dy);
                for(int p = 0; p < n_patches; p++){
                    for (int q = 0; q < n_patches; q++){
                        patch_i_x = r[i][0] + radius*cos(theta[i]+ p*2*M_PI/n_patches);
                        patch_i_y = r[i][1] - radius*sin(theta[i] + p*2*M_PI/n_patches);

                        uix = patch_i_x - r[i][0];
                        uiy = patch_i_y - r[i][1];
                        norm = sqrt(uix*uix + uiy*uiy);
                        uix = uix/norm;
                        uiy = uiy/norm;

                        patch_j_x = r[j][0] + radius*cos(theta[j]+ q*2*M_PI/n_patches);
                        patch_j_y = r[j][1] - radius*sin(theta[j] + q*2*M_PI/n_patches);

                        ujx = patch_j_x - r[j][0];
                        ujy = patch_j_y - r[j][1];
                        norm = sqrt(ujx*ujx + ujy*ujy);
                        ujx = ujx/norm;
                        ujy = ujy/norm;
                        
                        if( rx*uix + ry*uiy > cos_max && - rx*ujx - ry*ujy > cos_max  ){
                            V += -u0;
                            bonds[i][j] = 1;
                            bonds[j][i] = 1;
                        }
                    }
                    
                }
                
            }
        }
        
    }
    return V;
}

void find_hexagones(){
    int p_bonds = 0;
    int k = 0;
    int candidate = 0, node_to_find = 0;
    int ring[MAX_RINGSIZE];
    int bonded_part[3];
    int find_ring = 1;

    /* quadrati */
    ring_number[0]=0;
    /* pentagoni */
    ring_number[1]=0;
    /* esagoni */
    ring_number[2]=0;
    /* ettagoni */
    ring_number[3]=0;
    
    for (int p = 0; p < n_particles; p++){


        /* find first three particles for ring */
        p_bonds = 0; k = 0;
        for (int z = 0; z<3; z++){bonded_part[z] = -1;}
        for(int n = 0; n < n_particles; n++){
            if (bonds[n][p] == 1 && p != n) {
                bonded_part[p_bonds] = n;
                p_bonds ++;
            }
        }
        while (p_bonds >= 2 & k < p_bonds){
            /* initialize ring  */
            for (int i = 0; i < MAX_RINGSIZE; i++){
                ring[i] = -1;
            }
            /* use p as reference particle to find ring */
            ring[0] = p;
            ring[1] = bonded_part[k];
            ring[MAX_RINGSIZE-1] = bonded_part[(k+1) % n_patches];
            node_to_find = 2;

            /* create hexagone */
            find_ring = 1;
            while(find_ring && node_to_find < MAX_RINGSIZE){
                candidate = 0;

                /* find next node */
                while(ring[node_to_find] == -1 && candidate < n_particles ){
                    /* check if the candidate is valid */
                    if(bonds[ring[node_to_find-1]][candidate] && candidate != ring[node_to_find-2] && distance(r[candidate], r[ring[(node_to_find+4) % 7]], box) <= 2*r_cut*1.15){
                        ring[node_to_find] = candidate;
                    } else { candidate ++; }
                }
                /* if node wasn't find after checking all the particles*/
                if(ring[node_to_find] == -1){
                    find_ring = 0;
                /* if node was found, check if the ring is closed or find next node */
                } else if (bonds[ring[node_to_find]][ring[MAX_RINGSIZE-1]]) {
                    ring_number[node_to_find-2]++;
                    find_ring = 0;
                } else {
                    node_to_find ++;
                }
            }
            k++;
        }
    }
}

void compute_radial_distr(){
    
    double dist;
    int b;
    for(b = 0; b < n_particles; b ++){
        nn[b] = 0;
    }

    /* fill the histograms bins */
    for (int i = 0; i < n_particles; i++){
        for(int j = i+1; j < n_particles; j++){

            /* calculate the distance for each different pair of particles */
            dist = distance(r[i], r[j], box);
            if(dist <= 1.2*diameter){
                nn[i] ++;
                nn[j] ++;
            }
            
        }
    }

    /* compute average neighbours */
    avg_nn = 0.0;
    for(int n = 0; n < n_particles; n++){
        avg_nn += nn[n];
    }
    avg_nn /= n_particles;
}


int main(int argc, char* argv[]){

    double temperatures[] = {0.04, 0.06, 0.08, 0.10, 0.15, 0.18};
    double densities[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

    char buffer[128];
    //sprintf(buffer, "last_config/%.2lf_%.1lf.dat", temp, density);
    //sprintf(buffer, "%.2lfcoords%.1lf/coords_step0290000.dat", temp, density);

    
    FILE* data_file = fopen("data.csv", "w");
    
    if (data_file == NULL)
	{
		fprintf(stderr,"Error with file opening!\n");
		exit(EXIT_FAILURE);
	}

    dsfmt_seed(time(NULL));
    

    radius = 0.5 * diameter;

    if(NDIM == 3) { particle_volume = M_PI * pow(diameter, 3.0) / 6.0; }
    else if(NDIM == 2) { particle_volume = M_PI * pow(radius, 2.0); }
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }
    r_cut = 1.1*diameter;

    double box_vol = 1.0;
    for (int i = 0; i < NDIM; i++){
        box_vol *= box[i];
    }

    double rmax = box[0]/2;
    double dr = rmax/NBIN;
    
    for(int t = 0; t < 6; t++){
        for(int d = 0; d < 6; d++){
            sprintf(buffer, "last_config/%.2lf_%.1lf.dat", temperatures[t], densities[d]);
            read_data(buffer);
            energy = potential();
            find_hexagones();
            compute_radial_distr();

            double sum = 0.0;
            for (int i = 0; i < 4; i++){ sum += ring_number[i]/(i+4.0);}
            printf("pot = %lf, rings %lf, squares %lf, pentagones %lf, hexagones %lf, heptagones %lf nn %lf \n", energy, sum, ring_number[0]/4.0, ring_number[1]/5.0, ring_number[2]/6.0, ring_number[3]/7.0, avg_nn);
    
    
            fprintf(data_file, "%lf; %lf; %lf; %lf; %lf; %lf; %lf; %lf; %lf \n", temperatures[t], densities[d], energy, sum, ring_number[0]/4.0, ring_number[1]/5.0, ring_number[2]/6.0, ring_number[3]/7.0, avg_nn);
    
        }
    }
    
    
    fclose(data_file);

    return 0;
}