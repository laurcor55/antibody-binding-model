#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_SUBSTRATES 2

float sq(float input);

float sq(float input){
  return input*input;
}

void center_distance(float *distances, float *ligand_location, float substrate_locations[3][N_SUBSTRATES]);

void center_distance(float *distances, float *ligand_location, float substrate_locations[3][N_SUBSTRATES]){
  int ii, jj;
  for (ii=0; ii<3; ii++){
    for (jj=0; jj<N_SUBSTRATES; jj++){
      distances[jj] = sq((ligand_location[0] - substrate_locations[0][jj])) + sq((ligand_location[1] - substrate_locations[1][jj])) + sq((ligand_location[2] - substrate_locations[2][jj]));
    }
  } 
}


int main(){
  float ligand_location[3] = {0};
  float substrate_locations[3][N_SUBSTRATES] = {0};
  float distances[N_SUBSTRATES] = {0};
  int ii, jj;
  for (ii=0; ii<3; ii++){
    for (jj=0; jj<N_SUBSTRATES; jj++){
      printf("ligand: %f, substrate: %f\n", ligand_location[ii], substrate_locations[ii][jj]);
    }
  } 
  center_distance(distances, ligand_location, substrate_locations);

  for (ii=0; ii<N_SUBSTRATES; ii++){
    printf("distance %f\n", distances[ii]);
  } 
  return 0;
}