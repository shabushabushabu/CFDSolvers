float** compute_b(float dx,float dy,float rho,float dt, float** u, float** v, int m, int n);
float** pressure_poisson(float** p, float dx, float dy, float rho, float dt, float** u, float** v, int m, int n);
int cavity_flow(int nt, float** u, float** v, float dt, float dx, float dy, float** p, float rho, float nu, float ut, int m, int n);