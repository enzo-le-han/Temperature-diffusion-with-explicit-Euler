import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

path='' #### path to save image for then construct video


# Box parameters and time
# x
Lx = 20.0  # Box size x [0, L]
Nx = 60  # Number of points on the spatial grid, EVEN
dx = Lx / Nx  # Spatial step
x = np.linspace(0, Lx, Nx)

# y
Ly = 80.0  # Box size y [0, L]
Ny = 240  # Number of points on the spatial grid, EVEN
dy = Ly / Ny  # Spatial step
y = np.linspace(0, Ly, Ny)

# Time
duration = 8.0  # Duration of simulation
Nt = 200  # Number of time points
dt = duration / Nt  # Temporal step

# Diffusivity coefficients
lambda_ = 400  # Thermal conductivity (copper=400 W/mK)
lambda_source = 45  # Thermal conductivity of normalized source object (skin=0.45W/mK)
rho = 8930  # Density (copper=8930 kg/m^3)
Cp = 382  # Specific heat capacity (copper=382 J/kgK)

a = 1  # lambda / (Cp * rho) # Problem if too high
b = a * lambda_source / lambda_

# Temperature profile
sigmax = 10
sigmay = 10
x0 = int(np.around(Lx / 2))
y0 = int(np.around(Ly / 2))
xm = int(np.around(x0 - sigmax / 2))
xM = int(np.around(x0 + sigmax / 2))
ym = int(np.around(y0 - sigmay / 2))
yM = int(np.around(y0 + sigmay / 2))
T_init_max = 100

X, Y = np.meshgrid(x, y)

def gaussian_2d(X, Y, x0, y0, sigmax, sigmay):
    return T_init_max * np.exp(-((((X - x0) ** 2) / (2 * sigmax ** 2)) + (((Y - y0) ** 2) / (2 * sigmay ** 2))))

# Initialization
T0 = gaussian_2d(X, Y, x0, y0, sigmax, sigmay)
Tp = np.zeros((Ny, Nx), dtype=float)
vmax=0


# Explicit Euler
# stable if dt/dx² <= 1/2
alphax = a * dt / (dx ** 2)  # a = lambda (material)    /    rho.Cp (material)
alphay = a * dt / (dy ** 2)
alphas = b * dt / dx  # b = lambda (source)   /   rho.Cp (material)

print("alphax = ", alphax)
print("alphay = ", alphay)
print("alphas = ", alphas)

# 2D
Trace = np.zeros((Ny, Nx), dtype=float)
file_n = 0

s = 1
p = 0
for n in range(0, Nt):

    if n>Nt/2 : s=0 # Stop source temperature at the middle time


    T = Trace 
    for r in range(1,Nx-1):
        for z in range(1,Ny-1):
            Trace[z, r] = T[z, r] + alphay * (T[z-1, r]-2*T[z, r]+T[z+1, r]) + alphax * (T[z, r-1]-2*T[z, r]+T[z, r+1])

            if int(np.around(xm/dx)) <= r <= int(np.around(xM/dx)) : 
                if int(np.around(ym/dy)) <= z < int(np.around(yM/dy)) : Trace[z, r] += s*alphas*(T0[z, r]-T[z,r]) + p*alphas*(Tp[z, r]-T[z,r])

            if r==1:
                Trace[z, 0] = T[z, 0] + alphay * (T[z-1, 0]-2*T[z, 0]+T[z+1, 0]) + alphax * (-T[z, 0]+T[z, 1])
                if xm==0:
                    if int(np.around(ym/dy)) <= z < int(np.around(yM/dy)):
                        Trace[z, 0] += s*alphas*(T0[z, 0]-T[z,0]) + p*alphas*(Tp[z, 0]-T[z,0])
            if r==Nx-2:
                Trace[z, Nx-1] =T[z, Nx-1] + alphay * (T[z-1, Nx-1]-2*T[z, Nx-1]+T[z+1, Nx-1]) + alphax * (T[z, Nx-2]-T[z, Nx-1]) 
                if int(np.around(xM/dx)) == Nx :
                    if int(np.around(ym/dy)) <= z < int(np.around(yM/dy)):
                        Trace[z, Nx-1] += s*alphas*(T0[z, Nx-1]-T[z,Nx-1]) + p*alphas*(Tp[z, Nx-1]-T[z,Nx-1])

        Trace[0, r] = T[0, r] + alphay * (-T[0, r]+T[1, r]) + alphax * (T[0, r-1]-2*T[0, r]+T[0, r+1])
        Trace[Ny-1, r] =T[Ny-1, r] + alphay * (T[Ny-2, r]-T[Ny-1, r]) + alphax * (T[Ny-1, r-1]-2*T[Ny-1, r]+T[Ny-1, r+1])  
        if int(np.around(xm/dx)) <= r <= int(np.around(xM/dx)): 
            if ym==0 : Trace[0, r] += s*alphas*(T0[0, r]-T[0,r]) + p*alphas*(Tp[0, r]-T[0,r])       
            if yM==Ly : Trace[Ny-1, r] += s*alphas*(T0[Ny-1, r]-T[Ny-1,r]) + p*alphas*(Tp[Ny-1, r]-T[Ny-1,r])
        

    Trace[0, 0] = Trace[1,1]
    Trace[Ny-1, 0] = Trace[Ny-2, 1]
    Trace[0, Nx-1] = Trace[1, Nx-2]
    Trace[Ny-1, Nx-1] = Trace[Ny-2, Nx-2]

    Tmax = np.max(Trace)
    Tmin = np.min(Trace)
    vmin = 0
    if Tmax == 0:
        vmax = 5
    elif Tmax>vmax:
        vmax = Tmax


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(0.85*Ly/4, Lx/2), gridspec_kw={'height_ratios': [1, 0.84]})


    im = ax1.imshow(Trace.T, extent=[0, Ly, 0, Lx], aspect='auto', origin='lower', cmap='hot', vmin=0, vmax=vmax)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("top", aspect=vmax/30, size="7%", pad=0.5)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal', label='Temperature', ticks=[vmin, vmax])

    ax1.text(0.02, 0.95, f'Time = {n * dt:.2f} s', transform=ax1.transAxes, color='white', fontsize=12,verticalalignment='top')
    if s==1 : ax1.text(0.02, 0.90, f'Tsource = {T_init_max:.2f}°C', transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')
    else : ax1.text(0.02, 0.90, "Tsource = nothing", transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')
    ax1.text(0.02, 0.85, f'Tmax = {Tmax:.2f}°C', transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')
    ax1.text(0.02, 0.80, f'Tmin = {Tmin:.2f}°C', transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')
    ax1.text(0.02, 0.75, f'TsideL = {T[0, int(Nx/2)]:.2f}', transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')
    ax1.text(0.02, 0.70, f'TsideR = {T[Ny-1, int(Nx/2)]:.2f}', transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')
    ax1.text(0.02, 0.65, f'Tbottom = {T[int(Ny/2), 0]:.2f}', transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')
    ax1.text(0.02, 0.60, f'Thight = {T[int(Ny/2), Nx-1]:.2f}', transform=ax1.transAxes, color='white', fontsize=12, verticalalignment='top')


    # Gradient
    gradient = np.gradient(-Trace)
    arrow_spacing = 3
    arrow_scale = 2

    X, Y = np.meshgrid(y[::arrow_spacing], x[::arrow_spacing])

    V = gradient[1].T[::arrow_spacing, ::arrow_spacing]
    U = gradient[0].T[::arrow_spacing, ::arrow_spacing]

    ax2.quiver(X, Y, U, V, scale=arrow_scale, color='red', scale_units='xy', angles='xy')
    ax2.set_xlim(0, Ly)
    ax2.set_ylim(0, Lx)



    plt.tight_layout()        
    plt.savefig(path+'/'+("0"*(4-len(str(file_n))))+str(file_n)+".png")
    file_n += 1
    plt.close('all')

    print(str(n+1)+'/'+str(Nt))




# Compilation 
import os
from moviepy.editor import ImageSequenceClip

print(os.getcwd())
image_folder = path + '/'

os.chdir(image_folder)
 
images = [img for img in os.listdir(image_folder)
        if  img.endswith(".jpg") or
            img.endswith(".jpeg") or
            img.endswith("png")]
     
print(images) 
  
clip = ImageSequenceClip(images, fps = 1/dt)
clip.ipython_display(width = 360)




