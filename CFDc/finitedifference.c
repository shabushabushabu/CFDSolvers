#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectorutils.h"

float** central_ddx(float** f, float dx, int m, int n){
  float** diff = init_zero_vector(m,n);
  
  for (int i=1; i<m-1; ++i){
    for (int j=1; j<n-1; ++j){
      diff[i][j] = (f[i][j+1] - f[i][j-1])/(2.0*dx);
    } 
  }
  return diff;
}


float** central_ddy(float** f, float dy, int m, int n){
  float** diff = init_zero_vector(m,n);
  
  for (int i=1; i<m-1; ++i){
    for (int j=1; j<n-1; ++j){
      diff[i][j] = (f[i+1][j] - f[i-1][j])/(2.0*dy);
    } 
  }
  return diff;
}


float** laplace(float** f, float dx, float dy, int m, int n){
  float** diff = init_zero_vector(m,n);

  for (int i=1; i<m-1; ++i){
    for (int j=1; j<n-1; ++j){
      diff[i][j] = (f[i][j+1] - 2*f[i][j] + f[i][j-1])/pow(dx, 2) + (f[i+1][j] - 2*f[i][j] + f[i-1][j])/pow(dy, 2);
    } 
  }
  return diff;
}
