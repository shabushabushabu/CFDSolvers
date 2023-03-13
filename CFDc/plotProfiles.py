import numpy as np
import matplotlib.pyplot as plt


def fill_profile(filename, vec):
    with open(filename) as f:
        for i, line in enumerate(f):
            data = line.split(",")[:-1]
            for j, value in enumerate(data):
                vec[i][j] = float(value)


def plot_flow(lx,ly,nx,ny,u,v,p):
    x = np.linspace(0, lx, nx)
    y = np.linspace(0, ly, ny)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(12,8), dpi=100)
    plt.contourf(X, Y, p, alpha=0.3)  
    plt.streamplot(X, Y, u, v)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('2D-Cavity Flow')
    return fig

def main():
    N_GRID_X = 41
    N_GRID_Y = 41
    LENGTH_X = 2.0
    LENGTH_Y = 2.0

    u = np.zeros((N_GRID_X, N_GRID_Y))
    v = np.zeros((N_GRID_X, N_GRID_Y))
    p = np.zeros((N_GRID_X, N_GRID_Y))

    fill_profile("./results/u.txt", u)
    fill_profile("./results/v.txt", v)
    fill_profile("./results/p.txt", p)

    cavity_flow_plot = plot_flow(
        lx=LENGTH_X,
        ly=LENGTH_Y,
        nx=N_GRID_X,
        ny=N_GRID_Y,
        u=u,
        v=v,
        p=p)

    cavity_flow_plot.savefig("./results/cavityFlow.png", dpi=300)

main()