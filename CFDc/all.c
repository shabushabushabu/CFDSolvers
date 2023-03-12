#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float** init_zero_vector(int m, int n){
  float* values = calloc(m*n, sizeof(float));
  float** vec = malloc(m*sizeof(float*));
  
  for (int i=0; i<m; ++i){
    vec[i] = values + i*n;
  }
  return vec;
}

void copy_vector(float** newV, float** previousV, int m, int n){
  size_t size = sizeof(float)*n;
  for (int i=0; i<m; ++i){
    for (int j=0; j<n; ++j){
      newV[i][j] = previousV[i][j];
    }
  }
}

void print_vector(float** vec, int m, int n){
  for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            printf("%f,", vec[i][j]);
        printf("\n");
    }
}

int write_result(char* fileName, float** vec, int m, int n){
  FILE *outFile;
  outFile=fopen(fileName,"w+");
  if (outFile==NULL){
        printf("[Error] Failed to open output file %s", fileName);
        return -1;
    }

  for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            fprintf(outFile, "%f,", vec[i][j]);
        fprintf(outFile, "\n");
    }
  fclose(outFile);

  return 0;
}

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


float** compute_b(float dx,float dy,float rho,float dt, float** u, float** v, int m, int n){
  float** b = init_zero_vector(m,n);
  float** dudx = central_ddx(u, dx, m, n);
  float** dvdy = central_ddy(v, dy, m, n);
  
  float** dudy = central_ddy(u, dy, m, n);
  float** dvdx = central_ddx(v, dx, m, n);

  for (int i=1; i<m-1; ++i){
      for (int j=1; j<n-1; ++j){
        b[i][j] = -rho*(pow(dudx[i][j],2) + 2*dudy[i][j]*dvdx[i][j] + pow(dvdy[i][j],2)) + rho/dt*(dudx[i][j] + dvdy[i][j]);
      } 
    }

  free(dudx);
  free(dvdy);
  free(dudy);
  free(dvdx);

  return b;
}


float** pressure_poisson(float** p, float dx, float dy, float rho, float dt, float** u, float** v, int m, int n){
  float** pn = init_zero_vector(m,n);
  int maxIteration = 50;

  float** b = compute_b(dx,dy,rho,dt,u,v,m,n);

  for (int it=0; it<maxIteration; ++it){
    copy_vector(pn, p, m, n);

    for (int i=1; i<m-1; ++i){
      for (int j=1; j<n-1; ++j){
        p[i][j] = (pow(dy,2)*(pn[i][j+1] + pn[i][j-1]) + pow(dx,2)*(pn[i+1][j] + pn[i-1][j]) - b[i][j]*pow(dy,2)*pow(dx,2))/(2*pow(dx,2) + 2*pow(dy,2));

      } 
    }
    // B.C.
    for (int i=0; i<m; ++i){
      p[i][n-1] = p[i][n-2]; // dp/dx = 0 at x = 2 (right)
      p[i][0] = p[i][1]; // dp/dx = 0 at x = 0
    }
  
    for (int j=0; j<n; ++j){
      p[0][j] = p[1][j]; // dp/dy = 0 at y = 0
      p[m-1][j] = 0.0; // p = 0   at y = 2 (top)
    }
  }
  free(pn);

  return p;
}


int cavity_flow(int nt, float** u, float** v, float dt, float dx, float dy, float** p, float rho, float nu, float ut, int m, int n){
  float** un = init_zero_vector(m,n);
  float** vn = init_zero_vector(m,n);

  for (int it=0; it<nt; ++it){
    copy_vector(un, u, m, n);
    copy_vector(vn, v, m, n);

    p = pressure_poisson(p, dx, dy, rho, dt, un, vn, m, n);

    float** dudx = central_ddx(un, dx, m, n);
    float** dudy = central_ddy(un, dy, m, n);
    float** dpdx = central_ddx(p, dx, m, n);
    float** laplaceU = laplace(un, dx, dy, m, n);

    float** dvdx = central_ddx(vn, dx, m, n);
    float** dvdy = central_ddy(vn, dy, m, n);
    float** dpdy = central_ddy(p, dy, m, n);
    float** laplaceV = laplace(vn, dx, dy, m, n);
    for (int i=1; i<m-1; ++i){
        for (int j=1; j<n-1; ++j){
          u[i][j] = un[i][j] + dt*(-un[i][j]*dudx[i][j]-vn[i][j]*dudy[i][j]-1/rho*dpdx[i][j]+nu*laplaceU[i][j]);
          v[i][j] = vn[i][j] + dt*(-un[i][j]*dvdx[i][j]-vn[i][j]*dvdy[i][j]-1/rho*dpdy[i][j]+nu*laplaceV[i][j]);
        }
    }

    free(dudx);
    free(dudy);
    free(dpdx);
    free(laplaceU);
    free(dvdx);
    free(dvdy);
    free(dpdy);
    free(laplaceV);
    
    // B.C.
      for (int j=0; j<n; ++j){
      u[0][j] = 0.0;   // bottom
      v[0][j] = 0.0;
      u[m-1][j] = ut;  // top
      v[m-1][j] = 0.0; 
    }
    for (int i=0; i<m; ++i){
      u[i][0] = 0.0;   // left
      v[i][0] = 0.0;
      u[i][n-1] = 0.0; // right
      v[i][n-1] = 0.0; 
    }
  }

  free(un);
  free(vn);

  return 0;
}


int main(void) {

  // Mesh and parameter set-up
  int N_GRID_X = 41;
  int N_GRID_Y = 41;
  float LENGTH_X = 2.0;
  float LENGTH_Y = 2.0;

  float dx = LENGTH_X / (N_GRID_X - 1);
  float dy = LENGTH_Y / (N_GRID_Y - 1);

  int N_TIMESTEP = 100;
  float dt = 0.001;

  float DENSITY = 1.0;
  float KINEMATIC_VISCOSITY = 0.1;
  float U_TOP = 1.0;

  // I.C.
  float** u = init_zero_vector(N_GRID_X, N_GRID_Y);
  float** v = init_zero_vector(N_GRID_X, N_GRID_Y);
  float** p = init_zero_vector(N_GRID_X, N_GRID_Y);

  // Loop through timesteps
  cavity_flow(
    N_TIMESTEP, u, v, dt, dx, dy, p, DENSITY,
    KINEMATIC_VISCOSITY, U_TOP, N_GRID_X, N_GRID_Y);

  write_result("./results/u.txt", u, N_GRID_X, N_GRID_Y);
  write_result("./results/v.txt", v, N_GRID_X, N_GRID_Y);
  write_result("./results/p.txt", p, N_GRID_X, N_GRID_Y);

  return 0;
}
