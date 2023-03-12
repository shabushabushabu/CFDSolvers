#include <stdio.h>
#include <stdlib.h>
#include "vectorutils.h"
#include "fluiddynamics.h"

int main(void) {
    printf("Solving Cavity Fluid Flow using Finite-Difference Approach\n");

    // Mesh and parameter set-up
    const int N_GRID_X = 41;
    const int N_GRID_Y = 41;
    const float LENGTH_X = 2.0;
    const float LENGTH_Y = 2.0;

    const float dx = LENGTH_X / (N_GRID_X - 1);
    const float dy = LENGTH_Y / (N_GRID_Y - 1);

    const int N_TIMESTEP = 100;
    const float dt = 0.001;

    const float DENSITY = 1.0;
    const float KINEMATIC_VISCOSITY = 0.1;
    const float U_TOP = 1.0;

    // I.C.
    float** u = init_zero_vector(N_GRID_X, N_GRID_Y);
    float** v = init_zero_vector(N_GRID_X, N_GRID_Y);
    float** p = init_zero_vector(N_GRID_X, N_GRID_Y);

    // Loop through timesteps
    cavity_flow(
        N_TIMESTEP, u, v, dt, dx, dy, p, DENSITY,
        KINEMATIC_VISCOSITY, U_TOP, N_GRID_X, N_GRID_Y);

    int u_txt = write_result("./results/u.txt", u, N_GRID_X, N_GRID_Y);
    int v_txt = write_result("./results/v.txt", v, N_GRID_X, N_GRID_Y);
    int p_txt = write_result("./results/p.txt", p, N_GRID_X, N_GRID_Y);

    if (u_txt!=0 || u_txt!=0 || u_txt!=0) {
        printf("[ERROR] Failed to write output files\n");
    } else {
        printf("[SUCCESS] Wrote output files\n");
        return 0;
    }

}
