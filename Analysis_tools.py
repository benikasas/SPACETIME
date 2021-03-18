import sys, os, string, numba
import numpy as np
import pygal
from generate import betas, Nt, Nx, Ny, Nz, action, epsilon, thermal, border, magnitude_1
from plaq_v_cfg import Nstart, Nend
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import gauge_latticeqcd
import matplotlib as cm
import tools_v1
import Relative_tools
from pygal.style import Style

### NEEDS MORE COMMENTS, IT'S VERY MESSY BUT IT WORKS

beta=betas[0]
## Uses the settings to generate the name for the selected data from relevant file.

## Function to graph the plaquette values against configuration number
def grapher_plaq():
    U_infile = './dats/plaquette_from_1_to_440_v_epsilon_0.2_v_cfg_570_11x11x11x11_W_border_10_magnitude_1000000000000000000000000000000000000.dat'
    U_infile_2 = './dats/plaquette_from_1_to_440_v_epsilon_0.2_v_cfg_570_11x11x11x11_W_border_4_magnitude_10000000000000000000000000000000000000.dat'
    data=np.loadtxt(U_infile,dtype=float,delimiter=' ',skiprows=0)
    data_2=np.loadtxt(U_infile_2,dtype=float,delimiter=' ',skiprows=0)
    rendername = './images/plaquette_from_' + str(Nstart) + '_to_' + str(Nend) + '_v_epsilon_' + str(epsilon) + '_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1))
    c_style=Style(
        background='white',
        plot_background='white',
        foreground='black',
        label_font_size=14,
        major_label_font_size=14,
        title_font_size=25,
        legend_font_size=20,
        show_y_guides=True
    )
    lineplot=pygal.XY(stroke=False, show_legend=True, background='white', range=(0.5, 0.72), legend_at_bottom=True, style=c_style)
    lineplot.add('l-QCD', data[:])
    lineplot.add('l-QCD with deformations', data_2[:], stroke=False)
    # lineplot.xrange=(0, 800)
    lineplot.x_title='Configuration number'
    lineplot.y_title='Average plaquette'

    lineplot.render_to_file(rendername + '.svg')




def grapher_ac():
    U_infile_2 = './dats/S_v_cfg_570_11x11x11x11_W_border_4_magnitude_10000000000000000000000000000000000000.dat'
    U_infile = './dats/S_v_cfg_570_11x11x11x11_W_border_10_magnitude_1000000000000000000000000000000000000.dat'
    data=np.loadtxt(U_infile,dtype=float,delimiter=' ',skiprows=0)
    data_2=np.loadtxt(U_infile_2,dtype=float,delimiter=' ',skiprows=0)
    c_style=Style(
        background='white',
        plot_background='white',
        foreground='black',
        label_font_size=14,
        major_label_font_size=14,
        title_font_size=25,
        legend_font_size=20,
        show_y_guides=True
    )

    lineplot=pygal.XY(stroke=False, show_legend=True, background='white', legend_at_bottom=True)
    rendername = './images/S_v_cfg_' + str(Nstart) + '_to_' + str(Nend) + '_v_epsilon_' + str(epsilon) + '_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1))
    lineplot=pygal.XY(stroke=False, show_legend=True, background='white', legend_at_bottom=True, style=c_style)
    # for i in range(len(data)):
    data[:, 1] = data[:, 1] - 224944
    data_2[:, 1]= data_2[:, 1] - 224944
    lineplot.add('l-QCD', data[:])
    lineplot.add('l-QCD with deformations', data_2[:])
    # lineplot.range=(-80000, 10000)
    
    lineplot.x_title='Configuration number'
    lineplot.y_title='Action sum'
    lineplot.render_to_file(rendername + '.svg')

a=grapher_ac()

def grapher_rich():
    end=440
    initial=0
    dir_1 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    bet=betas[0]
    aa=tools_v1.fn_a(bet)
    kappa=Relative_tools.kappa_calc(aa)
    Ricci_sum=0.
    counter=0
    rich_array=np.zeros((140, 2))
    c_style=Style(
        background='white',
        plot_background='white',
        foreground='black',
        label_font_size=14,
        major_label_font_size=14,
        title_font_size=25,
        legend_font_size=20,
        show_y_guides=True
    )
    # print(rich_array)
    for i in range(end - thermal + initial):
        idx=i+thermal-initial
        input_file_1 = dir_1 + U_infile + str(idx)
        Ricci_mat=np.load(input_file_1)
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        Ricci_sum=Ricci_sum+Ricci_mat[t, x, y, z]
                                        counter=counter+1
        mean_Ricci=Ricci_sum/counter
        rich_array[i, 1]=mean_Ricci
        rich_array[i, 0]=i+300
    # print(rich_array)    
    rendername = './images/Mean_Rich_' + str(initial) + '_to_' + str(end) + '_v_epsilon_' + str(epsilon) + '_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1))
    lineplot=pygal.XY(stroke=False, show_legend=False, style=c_style)#, range=(10**(-37), 3*10**(-36)))
    lineplot.add('line', rich_array[:, :])
    # lineplot.title='Mean Ricci scalar curvature against configuration number' 
    lineplot.x_title='Configuration number'   
    lineplot.y_title='Mean Ricci scalar curvature'
    lineplot.render_to_file(rendername + '.svg')
    


## Function to run various analyses on the file
## Very outdated, needs a heavy upgrade
def statistician():
    # U_infile = './dats/plaquette_from_' + str(Nstart) + '_to_' + str(Nend) + '_v_epsilon_' + str(epsilon) + '_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '.dat'
    # U_infile='./dats/S_v_cfg_' + str(int(beta * 100)) + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_' + action + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '.dat'
    U_infile='./dats/plaq-old-edit.dat'
    plaq=np.loadtxt(U_infile,dtype=float,delimiter=' ',skiprows=50,usecols=(1,)) ## Loads plaquette data file
    totsum=np.sum(plaq) ## Finds the sum of the plaquette
    mean_val=totsum/len(plaq) ## Finds average of the plaquette
    varsum=0
    for i in range(len(plaq)):
        varsum=varsum+(plaq[i]-mean_val)**2 ## Finds the sum for the variance
    variance=1/(len(plaq)-1)*varsum ## Finds variance
    sdev=np.sqrt(variance)  ## Standard deviation of the sum
    mean_sdev=sdev/np.sqrt(len(plaq))   ## The standard deviation of the mean
    rel_sdev=sdev/mean_val  ## Relative standard deviation
    rel_error=mean_sdev/mean_val ## Relative error
    
    print('Acceptable value: ', mean_val, ' with standard deviation :', mean_sdev)


### Used to reduce 4D to 3D
### This is done by fixing two coordinates and plotting the Ricci curvature against the other two coordinates

def SP_slicer():
    end=440
    # beta=5.7
    dir_1 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    dir_2 = './Rich_images/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    bet=betas[0]
    aa=tools_v1.fn_a(bet)
    kappa=Relative_tools.kappa_calc(aa)
    for i in range(end - thermal):

        idx=i+thermal
        input_file = dir_1 + U_infile + str(idx)
        SP=np.load(input_file)
        x=range(len(SP))
        y=range(len(SP))
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        X, Y = np.meshgrid(x, y) 
        ax.plot_wireframe(X, Y, SP[border+1][:][:][border+1])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        # norm=plt.Normalize(SP.min(), grid.max())
        colors=plt.viridis()
        surf=ax.plot_surface(X, Y, SP[border+2][:][:][border+2]*kappa, rstride=2, cstride=2, cmap=cm.colors.ListedColormap((0,0,0,0)), antialiased=False, )
        m = plt.cm.ScalarMappable(surf.norm)
        surf.set_edgecolors(m.to_rgba(surf.get_array()))
        ax.set_zlim(0, 9)
        surf.set_facecolor([0]*4)
        ax.set_xlabel('X direction')
        ax.set_ylabel('Y direction')
        ax.set_zlabel('Ricci scalar curvature * ' +  r'$\kappa$')
        ax.view_init(elev=0, azim=150)
        # plt.show()
        if not os.path.exists(dir_2):
            os.mkdir(dir_2)
        out_file= dir_2 + U_infile + '_' + str(idx)
        fig.savefig(out_file)

### Same as above but for LQCD part
def U_slicer():
    end=440
    dir_1 = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    dir_2 = './Plaq_images/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    for i in range(end - thermal+50):

        idx=i+thermal-50
        input_file = dir_1 + U_infile + str(idx)
        U=np.load(input_file)
        x=range(len(U))
        y=range(len(U))
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        X, Y = np.meshgrid(x, y)
        grid=np.zeros((len(U),len(U))) 
        for alpha in range(len(U)):
            for beta in range(len(U)):
                grid[alpha][beta]=gauge_latticeqcd.fn_eval_point_S(U, border+2, alpha, beta, border+2, betas[0], u0=1.)

        # norm = plt.Normalize(grid.min(), grid.max())
        # print(norm)
        # colors = plt.viridis()
        # rcount, ccount, _ = colors.shape
         
        norm=plt.Normalize(grid.min(), grid.max())
        colors=plt.viridis()
        surf=ax.plot_surface(X, Y, grid, rstride=2, cstride=2, cmap=cm.colors.ListedColormap((0,0,0,0)), antialiased=False, )
        m = plt.cm.ScalarMappable(surf.norm)
        surf.set_edgecolors(m.to_rgba(surf.get_array()))
        surf.set_facecolor([0]*4)
        ax.set_zlim(10,35)
        ax.set_xlabel('X direction')
        ax.set_ylabel('Y direction')
        ax.set_zlabel('Action')
        ax.view_init(elev=0, azim=150)
        # surf.set_edgecolors(surf.to_rgba(surf._A))
        # surf.set_facecolors("white")
        # plt.show()
        if not os.path.exists(dir_2):
            os.mkdir(dir_2)
        out_file= dir_2 + U_infile + '_' + str(idx)
        fig.savefig(out_file)

### Combines the above two functions into a single one
def U_SP_fusion():
    end=500     ### Don't really need this, as without it, after running out of configurations, returns an error
    dir_1 = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    dir_2 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    dir_3 = './Fusion/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    initial=100   ### Number of configurations to plot before plotting the spacetime part
    bet=betas[0]
    aa=tools_v1.fn_a(bet)
    kappa=Relative_tools.kappa_calc(aa)
    for i in range(end - thermal + initial):
        idx=i+thermal-initial
        input_file_1 = dir_1 + U_infile + str(idx)
        U=np.load(input_file_1)
        if i >= initial:
            input_file_2 = dir_2 + U_infile + str(idx)
            SP=np.load(input_file_2)

        x=range(len(U))
        y=range(len(U))
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        X, Y = np.meshgrid(x, y)
        grid=np.zeros((len(U),len(U))) 
        for alpha in range(len(U)):
            for beta in range(len(U)):
                grid[alpha][beta]=gauge_latticeqcd.fn_eval_point_S(U, border+1, alpha, beta, border+1, betas[0], u0=1.)     ### Plots the action. Can be made to plot plaquette value at a point
        ax.plot_wireframe(X, Y, grid)
        
        if i>= initial:
            ax.plot_wireframe(X, Y, SP[border+1][:][:][border+1]*kappa, color='red')    ### Unsure about the kappa part
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        # plt.show()
        if not os.path.exists(dir_3):
            os.mkdir(dir_3)
        out_file= dir_3 + U_infile + '_' + str(idx)
        fig.savefig(out_file)
        plt.clf()


def avg_rich():
    end=440
    initial=-50
    dir_1 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    bet=betas[0]
    aa=tools_v1.fn_a(bet)
    kappa=Relative_tools.kappa_calc(aa)
    Ricci_sum=0.
    counter=0
    for i in range(end - thermal + initial):
        idx=i+thermal-initial
        input_file_1 = dir_1 + U_infile + str(idx)
        Ricci_mat=np.load(input_file_1)
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        Ricci_sum=Ricci_sum+Ricci_mat[t, x, y, z]
                                        counter=counter+1

    mean_Ricci=Ricci_sum/counter
    print('Ricci sum: ', Ricci_sum)
    print('Average Ricci scalar curvature: ', mean_Ricci)

    counter=0
    Ricc_var=0.
    for i in range(end - thermal + initial):
        idx=i+thermal-initial
        input_file_1 = dir_1 + U_infile + str(idx)
        Ricci_mat=np.load(input_file_1)
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        Ricc_var=Ricc_var+(Ricci_mat[t, x, y, z]-mean_Ricci)**2 
                                        counter=counter+1
    Ricc_var=Ricc_var/(counter-1)
    print('Ricci variance: ', Ricc_var)
    Ricc_sdev=np.sqrt(Ricc_var)
    print('Ricci standard deviation: ', Ricc_sdev)



def rich_plot():
    end=411
    initial=0
    dir_1 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
    bet=betas[0]
    aa=tools_v1.fn_a(bet)
    kappa=Relative_tools.kappa_calc(aa)
    Ricci_sum=0.
    for i in range(end - thermal + initial):
        idx=i+thermal-initial
        input_file_1 = dir_1 + U_infile + str(idx)
        Ricci_mat=np.load(input_file_1)
        for t in range(Nt):
            for x in range(Nx):
                for y in range(Ny):
                    for z in range(Nz):
                        if t >= border and t < Nt-border:
                            if x >= border and x < Nx-border:
                                if y >= border and y < Ny-border:
                                    if z >= border and z < Nz-border:
                                        Ricci_sum=Ricci_sum+Ricci_mat[t, x, y, z]

    print(Ricci_sum)




def lagrangian_test():
    dir_1 = './logs/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    dir_2 = './Rich/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    dir_3 = './Fusion/' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '/'
    U_infile = 'link_' + action + '_' + str(Nt) + 'x' + str(Nx) + 'x' + str(Ny) + 'x' + str(Nz) + '_b' + str(int(betas[0] * 100)) + '_border_' + str(border) + '_magnitude_' + str(int(magnitude_1)) + '_'
