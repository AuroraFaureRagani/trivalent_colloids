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

/* Initialization variables */
const int mc_steps = 500000;
const int output_steps = 100;
const int n_patches = 3;
const double density = 0.75;
const double diameter = 1.0;
const double cos_max = 0.99;
//0.174533;  //<<- da rivedere
const double delta = 0.1;
const double delta_theta = 0.02;
const double u0 = 1.0;
const double temp = 0.2; // temp 1/beta
const char* init_filename = "coords/coords_step0006800.dat";

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
int hexa = 0;

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

void read_data(void){

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

void initialize(){
    // int rows = 32;    // Number of rows in the square grid
    // int cols = 32;    // Number of columns in the square grid
    // for (int i = 0; i < n_particles; i++) {
    //         r[i][0] = (i % cols);
    //         r[i][1] = (i / cols);
    //         theta[i] = 0.0;
    //         //dsfmt_genrand()*2.0/3.0*M_PI;
    // }
    // box[0] = 32.0;
    // box[1] = 32.0;

    // int rows = 23;    // Number of rows in the square grid
    // int cols = 23;    // Number of columns in the square grid
    // double spacing = 1.2; // Spacing between circles

    // for (int i = 0; i < n_particles; i++) {
    //     r[i][0] = 1+ (i % cols) * spacing;
    //     r[i][1] = 1+ (i / cols) * spacing;
    //     theta[i] = dsfmt_genrand()*2.0/3.0*M_PI;
    // }
    // box[0] = rows * spacing;
    // box[1] = cols * spacing;

    int rows = 10;
    int cols = 10;
    double spacing = 2.0;
    double shift = 1.0;

    for (int i = 0; i < n_particles; i++) {
         r[i][0] = 1 + (i % cols) * spacing + ((i / cols) % 2) * shift;
         r[i][1] = 1 + (i / cols) * spacing;
         theta[i] = dsfmt_genrand()*2.0/3.0*M_PI;
    }
    box[0] = rows * spacing;
    box[1] = cols * spacing;
    // r[1][0] = r[0][0]+1.1;
    // theta[0] = 0.0;
    // theta[1] = 0.0;

    for(int i = 0; i < n_particles; i++){
        for(int j = 0; j < n_particles; j++){
            bonds[i][j] = 0;
        }
    }

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
            if(dist < diameter){ printf("ERROR: check for overlaps did not succeed\n"); }
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


int move_particle(void){

    double U_end, dU;
    double old_pos[NDIM];
    double random_displ[NDIM];

    /* select a random particle */
    int random_particle = (int)(dsfmt_genrand() * n_particles); //generates a random integer between 0 and n_particles-1
    
    /* select a random displacement in [-delta, delta] */
    for(int i = 0; i < NDIM; i++) {
        old_pos[i] = r[random_particle][i];

        /* move the particle */
        r[random_particle][i] += dsfmt_genrand()* 2 * delta - delta;;
        //r[random_particle][i] -= (r[random_particle][i] / box[i]) * box[i];
        if (r[random_particle][i] - box[i] > 0) {
            r[random_particle][i] -= box[i];
        }
        if (r[random_particle][i] < 0) {
            r[random_particle][i] += box[i];
        }
    }
    double d;

    /* check for overlaps */
    for (int i = 0; i < n_particles; i++){
        if(i != random_particle){

             /* calculate the distance d - considering boundary conditions*/
             double d = distance(r[random_particle], r[i], box);
             /* if the distance is not acceptable (overlapping occurs) reject the move */
            if (d < diameter) {
                 for(int j = 0; j<NDIM; j++){
                    r[random_particle][j] = old_pos[j];
                 }
                 return 0;
            }
        }
    }

    /* if no overlaps, compute delta energy and accept accourding to the boltzman factor*/
    U_end = potential();
    dU = U_end - energy;

    if(dU < 0.0 || dsfmt_genrand() < exp(-dU/temp)){
        /* accept move */
        energy = U_end;
        return 1;
    }

    /* if condition above not satisfied, reject move */
    for(int j = 0; j < NDIM; j++){
        r[random_particle][j] = old_pos[j];
    }
    return 0;

}


int rotate_particle(void){

    double alpha1, alpha2;
    double old_theta;
    double U_end, dU;
    /* select a random particle */
    int random_particle = (int)(dsfmt_genrand() * n_particles); //generates a random integer between 0 and n_particles-1
    
    /* select a random rotation in [-delta_theta, delta_theta] */
    old_theta = theta[random_particle];
    theta[random_particle] += dsfmt_genrand()* 2 * delta_theta - delta_theta;;
    theta[random_particle] -= (int)(theta[random_particle] / (2.0*M_PI/n_patches)) * (2.0*M_PI/n_patches);

    U_end = potential();
    dU = U_end - energy;

    if(dU < 0.0 || dsfmt_genrand() < exp(-dU/temp)){
        /* accept move */
        energy = U_end;
        return 1;
    }

    /* if condition above not satisfied, reject move */
    theta[random_particle] = old_theta;
    return 0;

}


void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    if (buffer == NULL)
	{
		fprintf(stderr,"Error with file opening!\n");
		exit(EXIT_FAILURE);
	}
    int d, n;
    fprintf(fp, "%d\n", n_particles + n_particles*n_patches);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n", 0.0, box[d]);
    }
    fprintf(fp, "%lf %lf\n", 0.0, 0.0);
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\t%lf\t%d\t%lf\n", 0.0, diameter, 1, theta[n]);
        for(int p=0; p<n_patches; p++){
            double patchX = r[n][0] + radius*cos(theta[n]+ p*2*M_PI/n_patches);
            double patchY = r[n][1] - radius*sin(theta[n] + p*2*M_PI/n_patches);
            fprintf(fp, "%f\t", patchX);
            fprintf(fp, "%f\t", patchY);
            fprintf(fp, "%f\t", 0.0);
            fprintf(fp, "%lf\t%d\n", diameter/2.5, 2);
        }
        
    }
    fclose(fp);
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

void set_density(void){
    /* scale the coordinates of the particles and the size of the box in order to reach the set packing fraction */

    /* calculate the volume of the box */
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    /* calculate new volume and scale factor */
    double target_volume = (n_particles * particle_volume) / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    /* update positions and sizes */
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) { box[d] *= scale_factor; }
}

void compute_radial_distr(double dr, double vol, double rmax){
    
    double dist;
    int b;

    /* initialize the array that represents the histogram */
    for(b = 0; b < NBIN; b ++){
        nhis[b] = 0;
    }
    for(b = 0; b < n_particles; b ++){
        nn[b] = 0;
    }

    /* fill the histograms bins */
    for (int i = 0; i < n_particles; i++){
        for(int j = i+1; j < n_particles; j++){

            /* calculate the distance for each different pair of particles */
            dist = distance(r[i], r[j], box);
            /* particular and very rare case */
            if(dist <= rmax){
                if(dist/dr - floor(dist/dr) == 0.0){ 
                b = dist/dr -1;
                nhis[b]+=2;
                } else {
                    // add 2 particle pairs (i,j) and (j,i) to the right bin
                    b = (int) floor(dist/dr);
                    nhis[b]+=2;
                }
            }
            if(dist <= 1.1*diameter){
                nn[i] ++;
                nn[j] ++;
            }
            
        }
    }

    /* compute the radial distribution g */
    double rr;
    double rho = n_particles / vol;
    double nid;
    
    for (b = 0; b < NBIN; b++){
        /* independent variable: distance from the reference particle */
        rr = (b + 0.5) * dr;
        /* particles at distance between r and r+dr in an ideal gas */
        nid = 4 * M_PI * rho * (pow(rr + dr, 3.0) - pow(rr, 3.0)) / 3.0;
        /* radial distribution function */
        g[b] = nhis[b] / (nid * n_particles);
    }

    /* compute average neighbours */
    avg_nn = 0.0;
    for(int n = 0; n < n_particles; n++){
        avg_nn += nn[n];
    }
    avg_nn /= n_particles;
}

void write_radial_distr(FILE* file){
    // writing the values of g in a line of the file adressed by the pointer file
    for (int b = 0; b < NBIN; b++){
        fprintf(file, "%lf; ", g[b]);
    }
    fprintf(file, "\n");
}

void find_hexagones(){
    int p_bonds = 0;
    int candidate = 0, node_to_find = 0;
    int hexagone[6];
    int bonded_part[3];
    int not_done_with_this_particle = 0;
    int find_hexagone = 1;

    /* initialize hexagone */
    hexa = 0;
    for (int i = 0; i < 6; i++){
        hexagone[i] = -1;
    }
    
    for (int p = 0; p < n_particles; p++){
        /* find first three particles for hexagone */
        p_bonds = 0;
        for(int n = 0; n < n_particles; n++){
            if (bonds[n][p] == 1 && p != n) {
                bonded_part[p_bonds] = n;
                p_bonds ++;
            }
        }
        if (p_bonds >= 2){
            hexagone[0] = p;
            hexagone[1] = bonded_part[0];
            hexagone[5] = bonded_part[1];
            if(p_bonds == 3){ not_done_with_this_particle = 1;}
            node_to_find = 2;

            /* create hexagone */
            find_hexagone = 1;
            while(find_hexagone && node_to_find < 5){
                candidate = 0;

                /* find next node */
                while(hexagone[node_to_find] == -1 && candidate < n_particles ){
                    /* check if the candidate is valide */
                    if(bonds[hexagone[node_to_find-1]][candidate] && candidate != hexagone[node_to_find-2] && distance(r[candidate], r[hexagone[(node_to_find+3) % 6]], box) <= 2*r_cut*1.05){
                        hexagone[node_to_find] = candidate;
                    } else { candidate ++; }
                }
                /* if node wasn't find */
                if(hexagone[node_to_find] == -1){
                    find_hexagone = 0;
                /* if node wasn find, go to the next */
                } else { 
                    node_to_find ++;
                }
            }
            if( find_hexagone && node_to_find == 5 ){ 
               if(bonds[hexagone[5], hexagone[4]]){
                        hexa ++; 
                    }
            }

            if(p_bonds == 3){
                hexagone[5] = bonded_part[2];
                hexagone[2] = hexagone[3] = hexagone[4] = -1;
                node_to_find = 2;

                /* create hexagone */
                find_hexagone = 1;
                while(find_hexagone && node_to_find < 5){
                    candidate = 0;

                    /* find next node */
                    while(hexagone[node_to_find] == -1 && candidate < n_particles ){
                        /* check if the candidate is valide */
                        if(bonds[hexagone[node_to_find-1]][candidate] == 1 && candidate != hexagone[node_to_find-2] && distance(r[candidate], r[hexagone[(node_to_find+3) % 6]], box) <= 2*r_cut*1.05){
                            hexagone[node_to_find] = candidate;
                        } else { candidate ++; }
                    }
                    /* if node wasn't find */
                    if(hexagone[node_to_find] == -1){
                        find_hexagone = 0;
                    /* if node wasn find, go to the next */
                    } else { 
                        node_to_find ++;
                    }
                }
                if( find_hexagone && node_to_find == 5 ){ 
                    if(bonds[hexagone[5], hexagone[4]]){
                        hexa ++; 
                    }
                }
            }
            for (int i = 0; i < 6; i++){ hexagone[i] = -1; }
        }
    }
}

int main(int argc, char* argv[]){

    //assert(density != 0.3);
    assert(diameter > 0.0);
    //e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));
    radius = 0.5 * diameter;

    if(NDIM == 3) { particle_volume = M_PI * pow(diameter, 3.0) / 6.0; }
    else if(NDIM == 2) { particle_volume = M_PI * pow(radius, 2.0); }
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }
    

    dsfmt_seed(time(NULL));
    initialize();
    //read_data();

    double box_vol = 1.0;
    for (int i = 0; i < NDIM; i++){
        box_vol *= box[i];
    }
    r_cut = 1.1*diameter;

    set_density();
    write_data(1);
    energy = potential();

    /* radial distr file */
    char buffer[128];
    sprintf(buffer, "%.2lfhisto%.1lf.csv", temp, density);
    FILE* hist_file = fopen(buffer, "w");
    
    if (hist_file == NULL)
	{
		fprintf(stderr,"Error with file opening!\n");
		exit(EXIT_FAILURE);
	}

    /* data file */
    char buffer2[128];
    sprintf(buffer2, "%.2lfmeasurement%.1lf.csv", temp, density);
    FILE* data_file = fopen(buffer2, "w");
    
    if (data_file == NULL)
	{
		fprintf(stderr,"Error with file opening!\n");
		exit(EXIT_FAILURE);
	}

    double rmax = box[0]/2;
    double dr = rmax/NBIN;
    // prima riga del file: x del grafico (i.e. r)
    double r;
    for (int b=0; b<NBIN; b++){
        r = (b + 0.5) * dr;
        fprintf(hist_file, "%lf; ", r);
    }
    fprintf(hist_file, "\n");

    // monte carlo simulation
    int k, h;
    int accepted_trans = 0;
    int accepted_rot = 0;
    int step, n;
    for(step = 0; step < mc_steps; ++step){

        for(n = 0; n < n_particles; ++n){
            k = move_particle();
            h = rotate_particle();
            accepted_trans += k;
            accepted_rot += h;
            }
        
        if(step % (output_steps*10) == 0){
            
            //write_data(step);
            //write_bonds(step);
            compute_radial_distr(dr, box_vol, rmax);
            //find_hexagones();

            write_radial_distr(hist_file);
            fprintf(data_file, "%d; %lf; %lf; %d\n", step, energy, avg_nn, hexa);

            printf("Step %d - Move acceptance: trans %.2lf \t rot %.2lf \t en: %.0lf \t avg nn: %.2lf \t hexa: %d\n", step, (double)accepted_trans / (n_particles * output_steps), (double)accepted_rot / (n_particles * output_steps), energy, avg_nn, hexa);
            accepted_trans = 0;
            accepted_rot = 0;
            
        }
        if(step % (output_steps*10) == 0){
            write_data(step);
        }
    }
    //fclose(hist_file);
    return 0;
}