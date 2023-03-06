import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def central_ddx(f, dx):
    """
    Compute df/dx with central difference 
    """
    diff = np.zeros_like(f)
    diff[1:-1,1:-1] = (f[1:-1,2:] - f[1:-1,:-2])/(2*dx)
    return diff

def central_ddy(f, dy):
    """
    Compute df/dy with central difference 
    """
    diff = np.zeros_like(f)
    diff[1:-1,1:-1] = (f[2:,1:-1] - f[:-2,1:-1])/(2*dy)
    return diff

def laplace(f, dx, dy):
    diff = np.zeros_like(f)
    diff[1:-1,1:-1] = (f[1:-1,2:] -2*f[1:-1,1:-1] +f[1:-1,:-2])/dx**2 + \
        (f[2:,1:-1] -2*f[1:-1,1:-1] +f[:-2,1:-1])/dy**2 
    return diff


def pressure_poisson(p, dx, dy, rho, dt, u, v):
    pn = np.empty_like(p)
    max_it = 50

    b = np.empty_like(p)
    b = - rho*(central_ddx(u, dx)**2 + 
               2*central_ddy(u, dy)*central_ddx(v, dx) + 
               central_ddy(v, dy)**2) + \
                rho/dt*(central_ddx(u, dx) + central_ddy(v, dy))

    for _ in range(max_it):
        pn = p.copy()
  
        p[1:-1,1:-1] = (dy**2*(pn[1:-1,2:] + pn[1:-1,:-2]) + \
                        dx**2*(pn[2:,1:-1] + pn[:-2,1:-1]) - \
                            b[1:-1, 1:-1]*(dy**2)*(dx**2)) /(2*dx**2 + 2*dy**2)
        
        p[:, -1] = p[:, -2] # dp/dx = 0 at x = 2 (right)
        p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
        p[-1, :] = 0.0      # p = 0     at y = 2 (top)

    return p


def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, ut):
    un = np.empty_like(u)
    vn = np.empty_like(v)

    for _ in range(nt):
        un = u.copy()
        vn = v.copy()

        p = pressure_poisson(p, dx, dy, rho, dt, un, vn)

        u = un + dt*(-central_ddx(un*un, dx) -central_ddy(vn*un, dy) - 1/rho*central_ddx(p, dx) + \
            nu*laplace(un,dx,dy))
        v = vn + dt*(-central_ddx(un*vn, dx) -central_ddy(vn*vn, dy) - 1/rho*central_ddy(p, dy) + \
            nu*laplace(vn,dx,dy))
        
        # B.C.
        u[0,:] = 0.0  # bottom
        v[0,:] = 0.0     
            
        u[-1,:] = ut # top
        v[-1,:] = 0.0
        
        u[:,0] = 0.0 # left
        v[:,0] = 0.0
        
        u[:,-1] = 0.0 # right
        v[:,-1] = 0.0
        
    return u, v, p


def plot_flow(lx,ly,nx,ny,u,v,p):
    x = np.linspace(0, lx, nx)
    y = np.linspace(0, ly, ny)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(12,8), dpi=100)
    # plotting the pressure field as a contour
    plt.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)  
    plt.colorbar()
    # plotting the pressure field outlines
    plt.contour(X, Y, p, cmap=cm.viridis)  
    # plotting velocity field
    # plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2]) 
    plt.streamplot(X, Y, u, v)
    plt.xlabel('X')
    plt.ylabel('Y')
    return fig



def main():
    # Mesh and parameter set-up
    N_GRID_X = 41
    N_GRID_Y = 41
    LENGTH_X = 2.0
    LENGTH_Y = 2.0

    dx = LENGTH_X / (N_GRID_X - 1)
    dy = LENGTH_Y / (N_GRID_Y - 1)

    N_TIMESTEP = 100
    dt = 0.001

    DENSITY = 1
    KINEMATIC_VISCOSITY = 0.1
    U_TOP = 1.0

    # Initial 
    u = np.zeros((N_GRID_Y, N_GRID_X))
    v = np.zeros((N_GRID_Y, N_GRID_X))
    p = np.zeros((N_GRID_Y, N_GRID_X))

    # Loop
    u, v, p = cavity_flow(
        nt=N_TIMESTEP, 
        u=u, 
        v=v, 
        dt=dt, 
        dx=dx, 
        dy=dy, 
        p=p, 
        rho=DENSITY, 
        nu=KINEMATIC_VISCOSITY, 
        ut=U_TOP)
    
    cavity_flow_plot = plot_flow(
        lx=LENGTH_X,
        ly=LENGTH_Y,
        nx=N_GRID_X,
        ny=N_GRID_Y,
        u=u,
        v=v,
        p=p)
    
    cavity_flow_plot.savefig("cavityFlow.png", dpi=300)


main()

    

