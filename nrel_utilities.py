import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math
from numpy import loadtxt
from numpy import savetxt

# Turn on interactive mode, if not already on
plt.ion()


def nrel_utilities_helloWorld():
    print('Hello World!')
    return

# Hi I'm a comment!

def paper_figures(path):
    substance = 'CO$_2$ '
    R_unit    = 'kg / m$^3$' 
    mu_unit   = '$\mu$Pa s'  
    lam_unit  = 'W / (m K)'

    refEoS    = 'SW'
    refMuCor  = 'LM'
    refLamCor = 'Huber'

    DV = [''] * 28
    DV[0]  = "DV_compressibility.txt"
    DV[1]  = "DV_density_SRK.txt"
    DV[2]  = "DV_density_SW.txt"
    DV[3]  = "DV_lambda_SRK_Chung.txt"
    DV[4]  = "DV_lambda_SRK_Huber.txt"
    DV[5]  = "DV_lambda_SW_Chung.txt"
    DV[6]  = "DV_lambda_SW_Huber.txt"
    DV[7]  = "DV_mu_SRK_Chung.txt"
    DV[8]  = "DV_mu_SRK_LM.txt"
    DV[9]  = "DV_mu_SW_Chung.txt"
    DV[10] = "DV_mu_SW_LM.txt"
    DV[11] = "DV_dmudrho_SRK.txt"
    DV[12] = "DV_dmudrho_SW.txt"
    DV[13] = "DV_dlamdrho_SRK.txt"
    DV[14] = "DV_dlamdrho_SW.txt"
    DV[15] = "DV_R_Delta.txt"
    DV[16] = "DV_mu_Delta.txt"
    DV[17] = "DV_mu_Delta_EoS.txt"
    DV[18] = "DV_mu_Delta_cor.txt"
    DV[19] = "DV_lam_Delta.txt"
    DV[20] = "DV_lam_Delta_EoS.txt"
    DV[21] = "DV_lam_Delta_cor.txt"
    DV[22] = "DV_mu_Delta_approx1.txt"
    DV[23] = "DV_mu_product.txt"
    DV[24] = "DV_mu_Delta_approx2.txt"
    DV[25] = "DV_lam_Delta_approx1.txt"
    DV[26] = "DV_lam_product.txt"
    DV[27] = "DV_lam_Delta_approx2.txt"

    category = [''] * 28
    category[0]  = 'data'
    category[1]  = 'data'
    category[2]  = 'data'
    category[3]  = 'data'
    category[4]  = 'data'
    category[5]  = 'data'
    category[6]  = 'data'
    category[7]  = 'data'
    category[8]  = 'data'
    category[9]  = 'data'
    category[10] = 'data'
    category[11] = 'derivative'
    category[12] = 'derivative'
    category[13] = 'derivative'
    category[14] = 'derivative'
    category[15] = 'error'
    category[16] = 'error'
    category[17] = 'error'
    category[18] = 'error'
    category[19] = 'error'
    category[20] = 'error'
    category[21] = 'error'
    category[22] = 'error'
    category[23] = 'error'
    category[24] = 'error'
    category[25] = 'error'
    category[26] = 'error'
    category[27] = 'error'

    dataType = [''] * 28
    dataType[0]  = 'compressibility'
    dataType[1]  = 'density'
    dataType[2]  = 'density'
    dataType[3]  = 'conductivity'
    dataType[4]  = 'conductivity'
    dataType[5]  = 'conductivity'
    dataType[6]  = 'conductivity'
    dataType[7]  = 'viscosity'
    dataType[8]  = 'viscosity'
    dataType[9]  = 'viscosity'
    dataType[10] = 'viscosity'
    dataType[11] = 'viscosity'
    dataType[12] = 'viscosity'
    dataType[13] = 'conductivity'
    dataType[14] = 'conductivity'
    dataType[15] = 'density'
    dataType[16] = 'viscosity'
    dataType[17] = 'viscosity'
    dataType[18] = 'viscosity'
    dataType[19] = 'conductivity'
    dataType[20] = 'conductivity'
    dataType[21] = 'conductivity'
    dataType[22] = 'viscosity'
    dataType[23] = 'viscosity'
    dataType[24] = 'viscosity'
    dataType[25] = 'conductivity'
    dataType[26] = 'conductivity'
    dataType[27] = 'conductivity'
        

    title = [''] * 28
    title[0]  = substance + 'Compressibility, $Z$ [dimensionless]'
    title[1]  = substance + '$\\rho_\\text{SRK}$ [' + R_unit + ']'
    title[2]  = substance + '$\\rho_\\text{' + refEoS + '}$ [' + R_unit + ']'
    title[3]  = substance + '$\lambda_\\text{SRK--Chung}$ [' + lam_unit + ']'
    title[4]  = substance + '$\lambda_\\text{SRK--' + refLamCor + '}$ [' + lam_unit + ']'
    title[5]  = substance + '$\lambda_\\text{' + refEoS + '--Chung}$ [' + lam_unit + ']'
    title[6]  = substance + '$\lambda_\\text{' + refEoS + '--' + refLamCor + '}$ [' + lam_unit + ']'
    title[7]  = substance + '$\mu_\\text{SRK--Chung}$ [' + mu_unit + ']'
    title[8]  = substance + '$\mu_\\text{SRK--' + refMuCor + '}$ [' + mu_unit + ']'
    title[9]  = substance + '$\mu_\\text{' + refEoS + '--Chung}$ [' + mu_unit + ']'
    title[10] = substance + '$\mu_\\text{' + refEoS + '--' + refMuCor + '}$ [' + mu_unit + ']'
    title[11] = substance + ('$\partial \mu_\\text{Chung} / \partial \\rho_\\text{SRK}$ [' 
            + mu_unit + ' / (' + R_unit + ')]')
    title[12] = (substance + '$\partial \mu_\\text{Chung} / \partial \\rho_\\text{' 
            + refEoS + '}$ [' + mu_unit + ' / (' + R_unit + ')]')
    title[13] = (substance + '$\partial \lambda_\\text{Chung} / \partial \\rho_\\text{SRK}$ [' 
            + lam_unit + ' / (' + R_unit + ')]')
    title[14] = (substance + '$\partial \lambda_\\text{Chung} / \partial \\rho_\\text{' 
            + refEoS + '}$ [' + lam_unit + ' / (' + R_unit + ')]')
    title[15] = (substance + '$\Delta \\rho = \\rho_\\text{SRK} - \\rho_\\text{' 
            + refEoS + '}$ [' + R_unit + ']')
    title[16] = (substance + '$\Delta \mu = \mu_\\text{Chung}(T,\\rho_\\text{SRK})' 
            ' - \mu_\\text{' + refMuCor + '}(T,\\rho_\\text{' + refEoS + '})$ [' + mu_unit + ']')
    title[17] = (substance + '$\\Delta \mu_\\text{EoS} = \\mu_\\text{' + refMuCor + '}(T,\\rho_\\text{SRK})' 
            ' - \mu_\\text{' + refMuCor + '}(T,\\rho_\\text{' + refEoS + '})$ [' + mu_unit + ']')
    title[18] = (substance + '$\Delta \mu_\\text{corr} = \mu_\\text{Chung}(T,\\rho_\\text{' + refEoS + '})' 
            ' - \mu_\\text{' + refMuCor + '}(T,\\rho_\\text{' + refEoS + '})$ [' + mu_unit + ']')
    title[19] = (substance + '$\Delta \lambda = \lambda_\\text{Chung}(T,\\rho_\\text{SRK})' 
            ' - \lambda_\\text{' + refLamCor + '}(T,\\rho_\\text{' + refEoS + '})$ [' + lam_unit + ']')
    title[20] = (substance + '$\Delta \lambda_\\text{EoS} = \lambda_\\text{' 
            + refLamCor + '}(T,\\rho_\\text{SRK})' ' - \lambda_\\text{' 
            + refLamCor + '}(T,\\rho_\\text{' + refEoS + '})$ [' + lam_unit + ']')
    title[21] = (substance + '$\Delta \lambda_\\text{corr} = \lambda_\\text{Chung}(T,\\rho_\\text{' 
            + refEoS + '})' ' - \lambda_\\text{' + refLamCor + '}(T,\\rho_\\text{' 
            + refEoS + '})$ [' + lam_unit + ']')
    title[22] = (substance + '$\Delta \mu_\\text{EoS} + \Delta \mu_\\text{corr}$ [' + mu_unit + ']')
    title[23] = (substance + '$\Delta \\rho \cdot (\partial \mu_\\text{Chung}' 
            ' / \partial \\rho_\\text{SRK})$ [' + mu_unit + ']')
    title[24] = (substance + '$\Delta \\rho \cdot (\partial \mu_\\text{Chung}' 
            ' / \partial \\rho_\\text{SRK}) + \Delta \mu_\\text{corr}$ [' + mu_unit + ']')
    title[25] = (substance + '$\Delta \lambda_\\text{EoS} + \Delta \lambda_\\text{corr}$ [' + lam_unit + ']')
    title[26] = (substance + '$\Delta \\rho \cdot (\partial \lambda_\\text{Chung}' 
            ' / \partial \\rho_\\text{SRK})$ [' + lam_unit + ']')
    title[27] = (substance + '$\Delta \\rho \cdot (\partial \lambda_\\text{Chung}' 
            ' / \partial \\rho_\\text{SRK}) + \Delta \lambda_\\text{corr}$ [' + lam_unit + ']')

    for x in title: 
        print(x)
        
    for i in range(len(DV)):

        if dataType[i] == 'compressibility':
            bounds = [0, 1.4]

        if dataType[i] == 'viscosity':
            if category[i] == 'data':
                bounds = [0, 140]
            elif category[i] == 'error':
                bounds = [-20, 20]
            elif category[i] == 'derivative':
                bounds = [-0.5, 0.5]
                
        elif dataType[i] == 'conductivity':
            if category[i] == 'data':
                bounds = [0, 2]
            elif category[i] == 'error':
                bounds = [-2, 2]
            elif category[i] == 'derivative':
                bounds = [-1e-4, 1e-4]

        elif dataType[i] == 'density':
            if category[i] == 'data':
                bounds = [0, 1000]
            elif category[i] == 'error':
                bounds = [-120, 120]

        plot_scalar2(path,DV[i],title[i],'CO2',category[i],bounds)

    return


def plot_scalar2(path,
	DV_filename,
    title,
    substance,
    category,
    bounds,
	fig_width=0,
	fig_height=0,
	cbar_shrink=0,
	cbar_aspect=0):

	# Parse data identifiers
    path = check_path(path)
	
	# Read temperature, pressure, and dependent variable
    x_axis = loadtxt(path + "IV_temperature.txt",delimiter=',')
    y_axis = loadtxt(path + "IV_pressure.txt",delimiter=',')
    DV     = loadtxt(path + DV_filename, comments='#',delimiter=",",unpack=False, skiprows=1)

	# Set variables based on state of zoom variable
    if fig_width   == 0: fig_width = 14
    if fig_height  == 0: fig_height = 14/4.5
    if cbar_shrink == 0: cbar_shrink = 0.72
    if cbar_aspect == 0: cbar_aspect = 7		

	# Plot pcolor data
    X,Y = np.meshgrid(x_axis,y_axis)
    fig,ax = plt.subplots(1,1)
    
    if category == 'error':
        '''cp = ax.pcolor(X,Y,np.fix(DV),norm=colors.Normalize(vmin=bounds[0],vmax=bounds[1]),cmap='seismic')'''
        cp = ax.pcolor(X,Y,DV,
                norm=colors.SymLogNorm(linthresh=1,linscale=1,vmin=bounds[0],vmax=bounds[1],base=2),
                cmap='seismic')
    if category == 'data' or category == 'derivative':
        cp = ax.pcolor(X,Y,DV,norm=colors.Normalize(vmin=bounds[0],vmax=bounds[1]))

    plt.title(title)
    plt.xlabel("Temperature [K]")
    plt.ylabel("Pressure [atm]")
    ax.set_aspect('equal', adjustable='box')
    fig.set_size_inches(fig_width,fig_height)    
    fig.colorbar(cp,shrink=cbar_shrink,aspect=cbar_aspect)

	# Save result
    output_filename = substance + '_' + DV_filename.split(".")[0] 
    plt.savefig(output_filename,dpi=600)
	
    return


def max_dmudrho_array():
	path = ['','']
	path[0] = "/Users/tfara/Results/Low_Fidelity_PT_Results-CO2/Data/"
	path[1] = "/Users/tfara/Results/Low_Fidelity_rhoT_Results-CO2/Data/"

	IV = ['','','']
	IV[0] = "CoolProp_IV_temperature.txt"
	IV[1] = "CoolProp_IV_pressure.txt"
	IV[2] = "CoolProp_IV_density.txt"	
	
	DV = ['','','']
	DV[0] = "CoolProp_DV_density.txt"
	DV[1] = "PelePhys_DV_density.txt"
	DV[2] = "PelePhys_DV_dmudrho(analytic).txt"

	# Load temperatures, pressures, and densities
	temperatures_IV = loadtxt(path[0] + IV[0], delimiter=',', skiprows=0)
	pressures_IV = loadtxt(path[0] + IV[1], delimiter=',', skiprows=0)
	density_IV = loadtxt(path[1] + IV[2], delimiter=',', skiprows=0)

	density_SW_array  = loadtxt(path[0] + DV[0], delimiter=',', skiprows=1)
	density_SRK_array = loadtxt(path[0] + DV[1], delimiter=',', skiprows=1)
	dmudrho_array = loadtxt(path[1] + DV[2], delimiter=',', skiprows=1)

	path =  "/Users/tfara/Results/Low_Fidelity_PT_Results-CO2/Data/"
	IV = ['','']
	IV[0] = "CoolProp_IV_temperature.txt"
	IV[1] = "CoolProp_IV_pressure.txt"

	result = np.zeros((np.size(pressures_IV),np.size(temperatures_IV)))
	j = 0
	for P in pressures_IV:
		i = 0
		for T in temperatures_IV:
			result[j][i] = max_dmudrho(T, P, temperatures_IV, pressures_IV, density_IV,density_SW_array,density_SRK_array,dmudrho_array)
			i = i + 1
		j = j + 1
		print(100 * j/np.size(pressures_IV))

	savetxt(path + "PelePhys_DV_dmudrho(max_along_C).txt",result,delimiter=',',header='CO2 dmudrho max along C [(uPa s) / (kg/m^3) ] PelePhys')

	return



def max_dmudrho(T,P,temperatures_IV,pressures_IV,density_IV,density_SW_array,density_SRK_array,dmudrho_array):
	# Locate indices for temperature and pressure
	temp_index = np.where(temperatures_IV == T)[0][0]
	pres_index = np.where(pressures_IV == P)[0][0]

	# Use temp/pres indices to lookup density computed by SW or SRK
	density_SW  = int(np.around(density_SW_array[pres_index][temp_index]))
	density_SRK = int(np.around(density_SRK_array[pres_index][temp_index]))

	# Gather indices of density corresponding to line segment between rho_SW and rho_SRK
	density_index = [0,0]
	density_index[0] = np.where(density_IV == density_SW)[0][0]
	density_index[1] = np.where(density_IV == density_SRK)[0][0]

	if (density_index[1] < density_index[0]):
		temp = density_index[1]
		density_index[1] = density_index[0]
		density_index[0] = temp

	dmudrho = np.zeros(density_index[1] - density_index[0] + 1)
	j = 0
	for i in range(density_index[0], density_index[1]+1):
		dmudrho[j] = dmudrho_array[i][temp_index]
		j = j + 1

	dmudrho_max = max(np.absolute(dmudrho))

	return dmudrho_max
		


########################################################################################
# CHANGING (T,P) <--> (T,RHO)	 
######################################################################################## 

def PT_to_rhoT(path,DV_filename):
	'''PT_to_rhoT(path,DV_filename) plots data gathered as DV(T,P) for some dependent
	variable DV stored in DV_filename and transforms the plot to the (T,rho) domain.
	'''
	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,DV_filename)

	# Read temperature, pressure, and dependent variable
	T  = loadtxt(path + source + "_IV_temperature.txt", comments='#',delimiter=",",unpack=False)
	P  = loadtxt(path + source + "_IV_pressure.txt", comments='#',delimiter=",",unpack=False)
	R  = loadtxt(path + source + "_DV_density.txt", comments='#',delimiter=',',unpack=False, skiprows=1)
	DV = loadtxt(path + DV_filename, comments='#',delimiter=',',unpack=False,skiprows=1)

	# Compute variables
	xlines = 5
	ylines = 10
	xskip  = math.floor(len(T)/ylines)
	yskip  = math.floor(len(P)/xlines)
	xmax   = xskip*ylines
	ymax   = yskip*xlines
	X,Y = np.meshgrid(T,P)
	fig, (ax1,ax2) = plt.subplots(2,1)

	# Plot input grid
	for i in range(xlines+1):
		ax1.plot(T[0:xmax],Y[i*yskip,0:xmax],color='white',linewidth=0.5,alpha=0.5)
	for i in range(ylines+1):
		ax1.plot(X[0:ymax,i*xskip],P[0:ymax],color='white',linewidth=0.5,alpha=0.5)
	cp1 = ax1.pcolor(X[0:ymax,0:xmax],Y[0:ymax,0:xmax],DV[0:ymax,0:xmax])

	# Plot outline of input grid
	ax1.plot(T[0:xmax],Y[0,0:xmax],color='black',linewidth=2,alpha=1)
	ax1.plot(T[0:xmax],Y[ymax,0:xmax],color='black',linewidth=2,alpha=1)
	ax1.plot(X[0:ymax,0],P[0:ymax],color='black',linewidth=2,alpha=1)
	ax1.plot(X[0:ymax,xmax],P[0:ymax],color='black',linewidth=2,alpha=1)
	
	# Plot output grid
	for i in range(xlines+1):
		ax2.plot(T[0:xmax],R[i*yskip,0:xmax],color='white',linewidth=0.5,alpha=0.5)
	for i in range(ylines+1):
		ax2.plot(X[0:ymax,i*xskip],R[0:ymax,i*xskip],color='white',linewidth=0.5,alpha=0.5)
	cp2 = ax2.pcolor(X[0:ymax,0:xmax],R[0:ymax,0:xmax],DV[0:ymax,0:xmax])

	# Plot outline of output grid
	ax2.plot(T[0:xmax],R[0,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(T[0:xmax],R[ymax,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,0],R[0:ymax,0],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,xmax],R[0:ymax,xmax],color='black',linewidth=2,alpha=1)

	# Size and label Plots
	title = species + ' ' + prop + ' ' + unit + ' (' + source + ')'
	fig.suptitle(title)
	fig.set_size_inches(12,7.5)
	ax1.set(xlabel="Temperature [K]",ylabel="Pressure [atm]")
	ax2.set(xlabel="Temperature [K]",ylabel="Density [kg m^-3]")

	margin = 0.02
	T_max = T[len(T[0:xmax])-1]; T_min = T[0]; T_delta = T_max-T_min
	P_max = P[len(P[0:ymax])-1]; P_min = P[0]; P_delta = P_max-P_min
	R_max = np.max(R[0:ymax,0:xmax]); R_min = np.min(R[0:ymax,0:xmax]); R_delta = R_max-R_min
	ax1.set_xlim((T_min - T_delta * margin,   T_max + T_delta * margin))
	ax1.set_ylim((P_min - P_delta * 2*margin, P_max + P_delta * 2*margin))
	ax2.set_xlim((T_min - T_delta * margin,   T_max + T_delta * margin))
	ax2.set_ylim((R_min - R_delta * 2*margin, R_max + R_delta * 2*margin))

	print(R_min); print(R_max);
	fig.colorbar(cp1,aspect=7)
	fig.colorbar(cp2,aspect=7)

	return fig


def PT_to_rhoT_overlay(path):
	'''PT_to_rhoT_overlay(path) transforms a rectangular region in the (T,P) domain
	and transforms it to the (T,rho) domain using the CoolProp and PelePhys transformations,
	and overlays the transformed domains.
	'''
	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,"CoolProp_DV_density.txt")

	# Read temperature, pressure, and dependent variable
	T    = loadtxt(path + "CoolProp_IV_temperature.txt", comments='#',delimiter=",",unpack=False)
	P    = loadtxt(path + "CoolProp_IV_pressure.txt", comments='#',delimiter=",",unpack=False)
	R0 = loadtxt(path + "CoolProp_DV_density.txt", comments='#',delimiter=',',unpack=False, skiprows=1)
	R1 = loadtxt(path + "PelePhys_DV_density.txt", comments='#',delimiter=',',unpack=False, skiprows=1)

	# Compute variables
	xlines = 10
	ylines = 3
	xskip  = math.floor(len(T)/ylines)
	yskip  = math.floor(len(P)/xlines)
	xmax   = xskip*ylines
	ymax   = yskip*xlines
	X,Y = np.meshgrid(T,P)
	fig, (ax1,ax2) = plt.subplots(2,1)

	# Plot input grid
	for i in range(xlines+1):
		ax1.plot(T[0:xmax],Y[i*yskip,0:xmax],color='black',linewidth=1,alpha=1)
	for i in range(ylines+1):
		ax1.plot(X[0:ymax,i*xskip],P[0:ymax],color='black',linewidth=1,alpha=1)

	# Plot outline of input grid
	ax1.plot(T[0:xmax],Y[0,0:xmax],color='black',linewidth=2,alpha=1)
	ax1.plot(T[0:xmax],Y[ymax,0:xmax],color='black',linewidth=2,alpha=1)
	ax1.plot(X[0:ymax,0],P[0:ymax],color='black',linewidth=2,alpha=1)
	ax1.plot(X[0:ymax,xmax],P[0:ymax],color='black',linewidth=2,alpha=1)
	
	# Plot CoolProp output grid
	for i in range(xlines+1):
		ax2.plot(T[0:xmax],R0[i*yskip,0:xmax],color='black',linewidth=1,alpha=1)
	for i in range(ylines+1):
		ax2.plot(X[0:ymax,i*xskip],R0[0:ymax,i*xskip],color='black',linewidth=1,alpha=1)

	# Plot outline of CoolProp output grid
	ax2.plot(T[0:xmax],R0[0,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(T[0:xmax],R0[ymax,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,0],R0[0:ymax,0],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,xmax],R0[0:ymax,xmax],color='black',linewidth=2,alpha=1)

	# Plot PelePhys output grid
	for i in range(xlines+1):
		ax2.plot(T[0:xmax],R1[i*yskip,0:xmax],color='red',linewidth=1,alpha=0.5)
	for i in range(ylines+1):
			ax2.plot(X[0:ymax,i*xskip],R1[0:ymax,i*xskip],color='red',linewidth=1,alpha=0.5)

	# Plot outline of PelePhys output grid
	ax2.plot(T[0:xmax],R1[0,0:xmax],color='red',linewidth=2,alpha=0.5)
	ax2.plot(T[0:xmax],R1[ymax,0:xmax],color='red',linewidth=2,alpha=0.5)
	ax2.plot(X[0:ymax,0],R1[0:ymax,0],color='red',linewidth=2,alpha=0.5)
	ax2.plot(X[0:ymax,xmax],R1[0:ymax,xmax],color='red',linewidth=2,alpha=0.5)

	# Size and label Plots
	title = species + ' ' + prop + ' ' + unit + ', CoolProp (black) vs. PelePhys (red)'
	fig.suptitle(title)
	fig.set_size_inches(12,7.5)
	ax1.set(xlabel="Temperature [K]",ylabel="Pressure [atm]")
	ax2.set(xlabel="Temperature [K]",ylabel="Density [kg m^-3]")
	margin = 0.02
	T_max = T[len(T[0:xmax])-1]; T_min = T[0]; T_delta = T_max-T_min
	P_max = P[len(P[0:ymax])-1]; P_min = P[0]; P_delta = P_max-P_min
	R_max = np.max(R0[0:ymax,0:xmax]); R_min = np.min(R0[0:ymax,0:xmax]); R_delta = R_max-R_min
	ax1.set_xlim((T_min - T_delta * margin,   T_max + T_delta * margin))
	ax1.set_ylim((P_min - P_delta * 2*margin, P_max + P_delta * 2*margin))
	ax2.set_xlim((T_min - T_delta * margin,   T_max + T_delta * margin))
	ax2.set_ylim((R_min - R_delta * 2*margin, R_max + R_delta * 2*margin))

	return fig


def PT_to_rhoT_overlay_only(path):
	'''PT_to_rhoT_overlay_only(path) transforms a rectangular region in the (T,P) domain
	and transforms it to the (T,rho) domain using the CoolProp and PelePhys transformations,
	and overlays the transformed domains.

	PT_to_rhoT_overlay_only plots *only* the transformed region; c.f. PT_to_rhoT_overlay,
	which plots the region in the (T,P) domain and the transformed region in the (T,rho) domain.
				   '''
	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,"CoolProp_DV_density.txt")

	# Read temperature, pressure, and dependent variable
	T  = loadtxt(path + "CoolProp_IV_temperature.txt", comments='#',delimiter=",",unpack=False)
	P  = loadtxt(path + "CoolProp_IV_pressure.txt", comments='#',delimiter=",",unpack=False)
	R0 = loadtxt(path + "CoolProp_DV_density.txt", comments='#',delimiter=',',unpack=False, skiprows=1)
	R1 = loadtxt(path + "PelePhys_DV_density.txt", comments='#',delimiter=',',unpack=False, skiprows=1)

	# Compute variables
	xlines = 12
	ylines = 3
	xskip  = math.floor(len(T)/ylines)
	yskip  = math.floor(len(P)/xlines)
	xmax   = xskip*ylines
	ymax   = yskip*xlines
	X,Y = np.meshgrid(T,P)
	fig, ax2 = plt.subplots(1,1)

	# Plot CoolProp output grid
	for i in range(xlines+1):
		ax2.plot(T[0:xmax],R0[i*yskip,0:xmax],color='black',linewidth=1,alpha=1)
	for i in range(ylines+1):
		ax2.plot(X[0:ymax,i*xskip],R0[0:ymax,i*xskip],color='black',linewidth=1,alpha=1)

	# Plot outline of CoolProp output grid
	ax2.plot(T[0:xmax],R0[0,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(T[0:xmax],R0[ymax,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,0],R0[0:ymax,0],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,xmax],R0[0:ymax,xmax],color='black',linewidth=2,alpha=1)

	# Plot PelePhys output grid
	for i in range(xlines+1):
		ax2.plot(T[0:xmax],R1[i*yskip,0:xmax],color='red',linewidth=1,alpha=0.5)
	for i in range(ylines+1):
		ax2.plot(X[0:ymax,i*xskip],R1[0:ymax,i*xskip],color='red',linewidth=1,alpha=0.5)

	# Plot outline of PelePhys output grid
	ax2.plot(T[0:xmax],R1[0,0:xmax],color='red',linewidth=2,alpha=0.5)
	ax2.plot(T[0:xmax],R1[ymax,0:xmax],color='red',linewidth=2,alpha=0.5)
	ax2.plot(X[0:ymax,0],R1[0:ymax,0],color='red',linewidth=2,alpha=0.5)
	ax2.plot(X[0:ymax,xmax],R1[0:ymax,xmax],color='red',linewidth=2,alpha=0.5)

	# Size and label Plots
	title = species + ' ' + prop + ' ' + unit + ', CoolProp (black) vs. PelePhys (red)'
	fig.suptitle(title)
	fig.set_size_inches(12,7.5)
	ax2.set(xlabel="Temperature [K]",ylabel="Density [kg m^-3]")

	margin = 0.02
	T_max = T[len(T[0:xmax])-1]; T_min = T[0]; T_delta = T_max-T_min
	P_max = P[len(P[0:ymax])-1]; P_min = P[0]; P_delta = P_max-P_min
	R_max = np.max(R0[0:ymax,0:xmax]); R_min = np.min(R0[0:ymax,0:xmax]); R_delta = R_max-R_min
	ax2.set_xlim((T_min - T_delta * margin,   T_max + T_delta * margin))
	ax2.set_ylim((R_min - R_delta * 2*margin, R_max + R_delta * 2*margin))
	ax2.set_aspect('equal', adjustable='box')

	return fig, ax2


def rhoT_to_PT(path,DV_filename):
	'''rhoT_to_PT(path,DV_filename) plots data gathered as DV(T,rho) for some dependent
	variable DV stored in DV_filename and transforms the plot to the (T,P) domain.

	NOTE: This function was copied from PT_to_rhoT and the meanings of all the variables
	were updated, but not the names. Specifically, the variable "P" holds density data and 
	the variable "R" holds pressure data.
	'''

	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,DV_filename)
	if source == "Error": source = "CoolProp"

	# Read temperature, pressure, and dependent variable
	T  = loadtxt(path + source + "_IV_temperature.txt", comments='#',delimiter=",",unpack=False)
	P  = loadtxt(path + source + "_IV_density.txt", comments='#',delimiter=",",unpack=False)
	R  = loadtxt(path + source + "_DV_pressure.txt", comments='#',delimiter=',',unpack=False, skiprows=1)
	DV = loadtxt(path + DV_filename, comments='#',delimiter=',',unpack=False,skiprows=1)

	# Compute variables
	xlines = 5
	ylines = 10
	xskip  = math.floor(len(T)/ylines)
	yskip  = math.floor(len(P)/xlines)
	xmax   = xskip*ylines
	ymax   = yskip*xlines
	X,Y = np.meshgrid(T,P)
	fig, (ax1,ax2) = plt.subplots(2,1)

	# Plot input grid
	for i in range(xlines+1):
		ax1.plot(T[0:xmax],Y[i*yskip,0:xmax],color='white',linewidth=0.5,alpha=0.5)
	for i in range(ylines+1):
		ax1.plot(X[0:ymax,i*xskip],P[0:ymax],color='white',linewidth=0.5,alpha=0.5)
		cp1 = ax1.pcolor(X[0:ymax,0:xmax],Y[0:ymax,0:xmax],DV[0:ymax,0:xmax])

	# Plot outline of input grid
	ax1.plot(T[0:xmax],Y[0,0:xmax],color='black',linewidth=2,alpha=1)
	ax1.plot(T[0:xmax],Y[ymax,0:xmax],color='black',linewidth=2,alpha=1)
	ax1.plot(X[0:ymax,0],P[0:ymax],color='black',linewidth=2,alpha=1)
	ax1.plot(X[0:ymax,xmax],P[0:ymax],color='black',linewidth=2,alpha=1)

	# Plot output grid
	for i in range(xlines+1):
		ax2.plot(T[0:xmax],R[i*yskip,0:xmax],color='white',linewidth=0.5,alpha=0.5)
	for i in range(ylines+1):
		ax2.plot(X[0:ymax,i*xskip],R[0:ymax,i*xskip],color='white',linewidth=0.5,alpha=0.5)
		cp2 = ax2.pcolor(X[0:ymax,0:xmax],R[0:ymax,0:xmax],DV[0:ymax,0:xmax])

	# Plot outline of output grid
	ax2.plot(T[0:xmax],R[0,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(T[0:xmax],R[ymax,0:xmax],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,0],R[0:ymax,0],color='black',linewidth=2,alpha=1)
	ax2.plot(X[0:ymax,xmax],R[0:ymax,xmax],color='black',linewidth=2,alpha=1)

	# Size and label Plots
	title = species + ' ' + prop + ' ' + unit + ' (' + source + ')'
	fig.suptitle(title)
	fig.set_size_inches(12,7.5)
	ax1.set(xlabel="Temperature [K]",ylabel="Density [kg m^-3]")
	ax2.set(xlabel="Temperature [K]",ylabel="Pressure [atm]")

	margin = 0.02
	T_max = T[len(T[0:xmax])-1]; T_min = T[0]; T_delta = T_max-T_min
	P_max = P[len(P[0:ymax])-1]; P_min = P[0]; P_delta = P_max-P_min
	R_max = np.max(R[0:ymax,0:xmax]); R_min = np.min(R[0:ymax,0:xmax]); R_delta = R_max-R_min
	ax1.set_xlim((T_min - T_delta * margin,   T_max + T_delta * margin))
	ax1.set_ylim((P_min - P_delta * 2*margin, P_max + P_delta * 2*margin))
	ax2.set_xlim((T_min - T_delta * margin,   T_max + T_delta * margin))
	ax2.set_ylim((R_min - R_delta * 2*margin, R_max + R_delta * 2*margin))

	print(R_min); print(R_max);
	fig.colorbar(cp1,aspect=7)
	fig.colorbar(cp2,aspect=7)

	return fig

########################################################################################




########################################################################################
# PLOTTING FUNCTIONS	 
######################################################################################## 

def plot_scalar(path,
	DV_filename,
	zoom=0,
	zoom_bound=200,
	fig_width=0,
	fig_height=0,
	cbar_shrink=0,
	cbar_aspect=0,
	title=''):
	'''plot_scalar(path,DV_filename,zoom,zoom_bound,fig_width,fight_height_cbar_shrink_cbar_aspect)
	creates a pcolor plot of a scalar stored as a .txt file at path + DV_filename. plot_scalar
	will automatically determine whether a (T,P) or (T,rho) plot should be created.

	The arguments zoom, zoom_bound, fig_width, fig_height, cbar_shrink, and cbar_aspect are
	all optional. 

	zoom = 0 by default, in which case the full range of data will be plotted. If zoom != 0, 
	then a square zoomed in region near the lower-left of the data (i.e. near the critical point 
	for the typical use case) will be plotted. 

	zoom_bound represents the size of the zoomed in square in units. It's value is set to 200
	by default. E.g. if the step size is 0.25, then a square of 800 x 800 data points will be
	plotted.

	fig_width and fig_height are given in inches. Default values of 14 x 14/4.5 for zoom == 0
	and 8 x 9 for zoom !=0 will be set depending on the state of the zoom variable. These
	defaults were chosen by experimentation for the best looking plot.

	cbar_shrink and cbar_aspect adjust the colorbar. Default values will be chosen depending on
	the state of the zoom variable.
	'''

	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,DV_filename)
	print(prop)
	
	# If the scalar represents an error, then the function's behavior changes
	if source == "Error":  err_flag = 1
	else: err_flag = 0
	if err_flag == 1: source = "CoolProp"

	# Determine Coordinates
	if   os.path.isfile(path + "IV_pressure.txt"): state = "P"
	elif os.path.isfile(path + "IV_density.txt"):  state = "D"
	else: print( "ERROR: Coordinate system (T,P) or (T,rho) cannot be determined.")

	# Read temperature, pressure, and dependent variable
	x_axis = loadtxt(path + "IV_temperature.txt",delimiter=',')
	if   state == "P": y_axis = loadtxt(path + "IV_pressure.txt",delimiter=',')
	elif state == "D": y_axis = loadtxt(path + "IV_density.txt",delimiter=',')
	DV = loadtxt(path + DV_filename, comments='#',delimiter=",",unpack=False, skiprows=1)
	'''DV1 = loadtxt(path + "Error(percent)_DV_muCP-muCPrho.txt",delimiter=',',skiprows=1)
	DV2 = loadtxt(path + "Error(percent)_DV_muCPrho-muPP.txt",delimiter=',',skiprows=1)
	DV = DV1 + DV2'''

	# Set variables based on state of zoom variable
	if zoom == 0 and state == "P":
		if fig_width   == 0: fig_width = 14
		if fig_height  == 0: fig_height = 14/4.5
		if cbar_shrink == 0: cbar_shrink = 0.72
		if cbar_aspect == 0: cbar_aspect = 7		
	elif zoom == 0 and state == "D":
		if fig_width   == 0: fig_width = 14
		if fig_height  == 0: fig_height = 7
		if cbar_shrink == 0: cbar_shrink = 0.93
		if cbar_aspect == 0: cbar_aspect = 10			

	elif zoom != 0:
		zoom_x_index = np.where(x_axis == x_axis[0]+zoom_bound)[0][0]
		zoom_y_index = np.where(y_axis == y_axis[0]+zoom_bound)[0][0]
		x_axis = x_axis[0:zoom_x_index]
		y_axis = y_axis[0:zoom_y_index]
		DV = DV[0:zoom_y_index,0:zoom_x_index]
	if fig_width   == 0: fig_width   = 8
	if fig_height  == 0: fig_height  = 7
	if cbar_shrink == 0: cbar_shrink = 0.92
	if cbar_aspect == 0: cbar_aspect = 15		

	# Plot pcolor data
	X,Y = np.meshgrid(x_axis,y_axis)
	fig,ax = plt.subplots(1,1)
	'''cp = ax.pcolor(X,Y,DV,norm=colors.Normalize(vmin=-20,vmax=20),cmap='tab20b')'''
	cp = ax.pcolor(X,Y,np.fix(DV / 1) * 1,norm=colors.Normalize(vmin=-25,vmax=25),cmap='seismic')
	'''cp = ax.pcolor(X,Y,np.fix(DV / 0.01) * 0.01,norm=colors.Normalize(vmin=-20,vmax=20),cmap='seismic')'''
	'''cp = ax.pcolor(X,Y,np.fix(DV / 0.1) * 0.1,norm=colors.Normalize(vmin=-.05,vmax=.05),cmap='seismic')'''
	'''cp = ax.pcolor(X,Y,DV,norm=colors.Normalize(vmin=-25,vmax=25))'''
	'''cp = ax.pcolor(X,Y,DV,norm=colors.SymLogNorm(linthresh=1,linscale=1,vmin=-25,vmax=25,base=10))'''
	'''cp = ax.pcolor(X,Y,DV,norm=colors.Normalize(vmin=-3,vmax=3),cmap='seismic')'''
	'''cp = ax.pcolor(X,Y,DV)'''
	if title == '':
		if err_flag == 1: title = species + " " + prop + " "
		else:             title = species + " " + prop + " " + unit + " (" + source + ")"
	plt.title(title)
	plt.xlabel("Temperature [K]")
	if state == "P":
		plt.ylabel("Pressure [atm]")
	elif state == "D":
		plt.ylabel("Density [kg m^-3]")
	ax.set_aspect('equal', adjustable='box')
	fig.set_size_inches(fig_width,fig_height)    
	fig.colorbar(cp,shrink=cbar_shrink,aspect=cbar_aspect)

	'''
	# Plot maxima for specific heat
	if prop == "Specific Heat":
		x_max, y_max = capture_row_max(DV,x_axis,y_axis)
		x_max, y_max = capture_row_max(DV[:,:50],x_axis,y_axis)
		plt.plot(x_max,y_max,"white",linewidth=1)
		plt.show()
	'''

	"""
	# Plot supercritical boundary for (T,rho) plot
	if state == "D":
		if   zoom == 0: plot_rhoT_isobar(73,ax,source)
		elif zoom != 0: plot_rhoT_isobar(73,ax,source,zoom_bound*4-1)

	plot_rhoT_isobar(100,ax,source)	
	plot_rhoT_isobar(150,ax,source)	
	plot_rhoT_isobar(200,ax,source)	
	plot_rhoT_isobar(250,ax,source)	
	plot_rhoT_isobar(300,ax,source)	
	plot_rhoT_isobar(350,ax,source)	
	plot_rhoT_isobar(400,ax,source)	
	"""

	# Save result
	if err_flag == 1: output_filename = "Error_" + '_'.join([species] + prop.split(" "))
	else:             output_filename = '_'.join([source,species] + prop.split(" "))
	plt.savefig(output_filename,dpi=600)
	
	return


def plot_all_coolprop_results(path,zoom=0):
	'''plot_all_coolprop_results(path) plots and saves all results for all DV's measured by 
	CoolProp in a given directory.
	''' 
	path = check_path(path)

	# Check for (T,P) or (T,rho) plot
	if   os.path.isfile(path + "CoolProp_IV_pressure.txt"): state = "P"
	elif os.path.isfile(path + "CoolProp_IV_density.txt"):  state = "R"
	else: print( "ERROR: (T,P) vs. (T,rho) plot cannot be determined.")

	# Compile DV filenames
	DV = ['','','','','']
	DV[0] = "compressibility"
	DV[1] = "conductivity"
	DV[2] = "specificHeat"
	DV[3] = "viscosity"
	if   state == "P": DV[4] = "density"
	elif state == "R": DV[4] = "pressure"

	for i in range(len(DV)):
		plot_scalar(path,"CoolProp_DV_" + DV[i] + ".txt",zoom)

	return


def plot_all_pelephys_results(path,zoom=0):
	'''plot_all_coolprop_results(path) plots and saves all results for all DV's measured by 
	CoolProp in a given directory.
	''' 
	path = check_path(path)

	# Check for (T,P) or (T,rho) plot
	if   os.path.isfile(path + "PelePhys_IV_pressure.txt"): state = "P"
	elif os.path.isfile(path + "PelePhys_IV_density.txt"):  state = "R"
	else: print( "ERROR: (T,P) vs. (T,rho) plot cannot be determined.")

	# Compile DV filenames
	DV = ['','','','','']
	DV[0] = "compressibility"
	DV[1] = "conductivity"
	DV[2] = "specificHeat"
	DV[3] = "viscosity"
	if   state == "P": DV[4] = "density"
	elif state == "R": DV[4] = "pressure"

	for i in range(len(DV)):
		plot_scalar(path,"PelePhys_DV_" + DV[i] + ".txt",zoom)

	return


def plot_all_error_results(path,zoom=0,err_type='percent'):
	'''plot_all_error_results(path,err_type) plots and saves all error results saved
	in the directory given by PATH.

	The ERR_TYPE variable allows to specify the types of error that should be 
	plotted. Allowed values are 'percent', 'absolute', 'absolutePercent', or 'all'.
	ERR_TYPE defaults to 'percent'.
	''' 
	path = check_path(path)

	# Check for (T,P) or (T,rho) plot
	if   os.path.isfile(path + "PelePhys_IV_pressure.txt"): state = "P"
	elif os.path.isfile(path + "PelePhys_IV_density.txt"):  state = "R"
	else: print( "ERROR: (T,P) vs. (T,rho) plot cannot be determined.")

	# Compile DV filenames
	DV = ['','','','','']
	DV[0] = "compressibility"
	DV[1] = "conductivity"
	DV[2] = "specificHeat"
	DV[3] = "viscosity"
	if   state == "P": DV[4] = "density"
	elif state == "R": DV[4] = "pressure"

	if err_type == 'percent' or err_type == 'all':
		for i in range(len(DV)): plot_scalar(path,"Error(percent)_DV_" + DV[i] + ".txt",zoom)
	if err_type == 'absolute' or err_type == 'all':
		for i in range(len(DV)): plot_scalar(path,"Error(absolute)_DV_" + DV[i] + ".txt",zoom)
	if err_type == 'percentPercent' or err_type == 'all':
		for i in range(len(DV)): plot_scalar(path,"Error(absolutePercent)_DV_" + DV[i] + ".txt",zoom)

	return



def plot_all_results(path,zoom=0,source='all',err_type='percent'):
	'''plot_all_results(path,zoom,source) plots all results stored in a given directory.
	
	The ZOOM and SOURCE arguments are optional. If zoom == 0, then results will be plotted
	across the entire stored domain. If zoom != 0, then a 200 unit x 200 unit square that
	is zoomed in to the lower left of the data (i.e. zoomed in around the critical point
	for standard usage) will be plotted.

	The SOURCE argument allows the user to specify a specific data source. Allowed values
	are 'all', 'coolprop', 'pelephys', 'both', or 'error'.

	The ERR_TYPE argument allows to specify the type of error reuslts that should be
	plotted. Allowed values are 'percent', 'absolute', 'absolutePercent', or 'all'. By 
	default, err_type == 'percent'.
	'''
	path = check_path(path)

	if source == 'all' or source == 'coolprop' or source == 'both':
		plot_all_coolprop_results(path,zoom)
	if source == 'all' or source == 'pelephys' or source == 'both':
		plot_all_pelephys_results(path,zoom)
	if source == 'all' or source == 'error':
		plot_all_error_results(path,zoom)

	return



def plot_rhoT_isobar(P,ax,source,Tmax_index=6782,rho_max=10000):
	path = "/Users/tfara/Results/High_Fidelity_PT_Results-CO2/Data/"
	T = loadtxt(path + source + "_IV_temperature.txt",delimiter=',')
	R = loadtxt(path + source + "_DV_density.txt",delimiter=',',skiprows=1)
	P_index = (P - 73)*4

	R[R > rho_max] = np.nan

	#ax.plot(T[0:Tmax_index],R[P_index,0:Tmax_index],linewidth=1,color='white')
	ax.plot(T[0:Tmax_index],R[P_index,0:Tmax_index],linewidth=1,color='white')

	return

  
########################################################################################





########################################################################################
# ERROR FUNCTIONS	 
######################################################################################## 

def compute_percent_error(path,DV_type):
	"""compute_percent_error(path,DV) computes percent error between the values of DV
	reported by CoolProp and PelePhys. The result is saved in the directory indicated
	by path. Error is reported as 

			100 * (PelePhys - CoolProp) / CoolProp 

	so that positive errors indicate the PelePhysics value is greater than the corresponding 
	CoolProp value, and vice versa for negative errors.
	"""	
	# Parse data identifiers
	path = check_path(path)
	species, prop, unit, source = read_title(path,"CoolProp_DV_" + DV_type + ".txt")	

	# Read dependent variable
	DV = ['','']
	DV[0] = loadtxt(path + "CoolProp_DV_" + DV_type + ".txt",delimiter=',',skiprows=1)
	DV[1] = loadtxt(path + "PelePhys_DV_" + DV_type + ".txt",delimiter=',',skiprows=1)

	# Compute error
	error = 100 * (DV[1] - DV[0]) / DV[0]

	output_filename = "Error(percent)_DV_" + DV_type + ".txt"
	header = species + " " + DV_type.capitalize() + " " + "Percent Error [dimensionless] Error"
	np.savetxt(output_filename, error, delimiter=',', header = header, comments = '')

	return


def test_function():
	print("Hello World")
	return

def compute_percent_error_general(path,input1,input2):
	"""compute_percent_error2(path,input1,input2) computes percent error between any two
	inputs. Error is reported as

		100 * (input2 - input1) / input1

	so that error is signed.
	"""
	# Parse data identifiers
	path = check_path(path)
	species, prop, unit, source = read_title(path,input1)
	
	# Read dependent variable
	DV = ['','']
	DV[0] = loadtxt(path + input1,delimiter=',',skiprows=1)
	DV[1] = loadtxt(path + input2,delimiter=',',skiprows=1)

	# Compute and save error
	error = 100 * (DV[1] - DV[0]) / DV[0]
	output_filename = "Error(percent)_DV.txt"
	header = species + " " + "Viscosity" + " " + "Percent Error [dimensionless] Error"
	np.savetxt(output_filename, error, delimiter=',', header = header, comments = '')

	return


def compute_absolute_percent_error(path,DV_type):
	"""compute_percent_error(path,DV) computes absolute percent error between the values of DV
	reported by CoolProp and PelePhys. The result is saved in the directory indicated
	by path. Error is reported as 

			|PelePhys - CoolProp|/|CoolProp| 

	so that positive errors indicate the PelePhysics value is greater than the corresponding 
	CoolProp value, and vice versa for negative errors.
	"""	
	# Parse data identifiers
	path = check_path(path)
	species, prop, unit, source = read_title(path,"CoolProp_DV_" + DV_type + ".txt")	

	# Read dependent variable
	DV = ['','']
	DV[0] = loadtxt(path + "CoolProp_DV_" + DV_type + ".txt",delimiter=',',skiprows=1)
	DV[1] = loadtxt(path + "PelePhys_DV_" + DV_type + ".txt",delimiter=',',skiprows=1)

	# Compute error
	error = (DV[1] - DV[0])/DV[0]
	error = np.abs(error)

	output_filename = "Error(absolutePercent)_DV_" + DV_type + ".txt"
	header = species + " " + DV_type.capitalize() + " " + "Absolute Error [dimensionless] Error"
	np.savetxt(output_filename, error, delimiter=',', header = header, comments = '')

	return


def compute_absolute_error(path,DV_type):
	"""compute_percent_error(path,DV) computes absolute error between the values of DV
	reported by CoolProp and PelePhys. The result is saved in the directory indicated
	by path. Error is reported as |PelePhys - CoolProp| 
	"""	
	# Parse data identifiers
	path = check_path(path)
	species, prop, unit, source = read_title(path,"CoolProp_DV_" + DV_type + ".txt")	

	# Read dependent variable
	DV = ['','']
	DV[0] = loadtxt(path + "CoolProp_DV_" + DV_type + ".txt",delimiter=',',skiprows=1)
	DV[1] = loadtxt(path + "PelePhys_DV_" + DV_type + ".txt",delimiter=',',skiprows=1)

	# Compute error
	error = np.abs(DV[1] - DV[0])

	output_filename = "Error(absolute)_DV_" + DV_type + ".txt"
	header = species + " " + DV_type.capitalize() + " " + "Absolute Percent Error [dimensionless] Error"
	np.savetxt(output_filename, error, delimiter=',', header = header, comments = '')

	return


def compute_all_error(path,err_type='percent'):
	"""compute_all_error(path,err_type) computes and saves error for all dependent variables
	stored in the directory indicated by PATH.

	ERR_TYPE is an optional argument that can be set to 'all', 'percent', 'absolute', or
	'absolutePercent'. By default, it is set to 'percent'.
	"""
	path = check_path(path) 

	# Check for (T,P) or (T,rho) plot
	if   os.path.isfile(path + "CoolProp_IV_pressure.txt"): state = "P"
	elif os.path.isfile(path + "CoolProp_IV_density.txt"):  state = "R"
	else: print( "ERROR: (T,P) vs. (T,rho) plot cannot be determined.")
	
	# Compile DV filenames
	DV = ['','','','','']
	DV[0] = "compressibility"
	DV[1] = "conductivity"
	DV[2] = "specificHeat"
	DV[3] = "viscosity"
	if   state == "P": DV[4] = "density"
	elif state == "R": DV[4] = "pressure"

	if err_type == 'percent' or err_type == 'all':
		for i in range(len(DV)): compute_percent_error(path,DV[i])
	if err_type == 'absolute' or err_type == 'all':
		for i in range(len(DV)): compute_absolute_error(path,DV[i])
	if err_type == 'absolutePercent' or err_type == 'all':
		for i in range(len(DV)): compute_absolute_percent_error(path,DV[i])

	return


def grid_convergence(path):
	path = check_path(path)

	
	return
######################################################################################## 





######################################################################################## 
# UTILITY FUNCTIONS
######################################################################################## 

def read_title(path,filename):
    """read_title(path,filename) parses the first line of text from a data file 
    and captures the chemical species, physical property, units, and the source
    of the data (i.e. CoolProp or PelePhysics).

    The first line of text must be formatted as 'species property [units] source'
    """
    with open(path + filename) as file:
        temp = file.readline().strip('\n')
        species = temp.split(" [")[0].split(" ",1)[0]
        prop    = temp.split(" [")[0].split(" ",1)[1]   
        unit    = "[" + temp.split("[")[1].split("]")[0] + "]"
        source  = temp.split("] ")[1]
    return species, prop, unit, source
        

def read_titles(path1,DV_filename1,path2,DV_filename2):
    '''read_titles(path1,DV_filename1,path2,DV_filename2) is like an overloaded version
    of read_title(); it reads the first line of text from two data files and compiles the
    paths, filenames, chemical species, physical property, units, and sources of the data
    as arrays.
    '''
    # Parse data identifiers
    DV_filename = ['']*2; path = ['']*2
    species = ['']*2; prop = ['']*2; unit = ['']*2; source = ['']*2
    path[0] = path1; path[1] = path2;
    DV_filename[0] = DV_filename1; DV_filename[1] = DV_filename2;

    for i in range(2):
        path[i] = check_path(path[i])
        (species[i],prop[i],unit[i],source[i]) = read_title(path[i],DV_filename[i])

    return path, DV_filename, species, prop, unit, source


def check_path(path):
    '''check_path(path) checks that the path begins and ends with the backslash character
    '''
    if path[-1] != "/":
        path = path + "/"
    if path[0] != "/":
        path = "/" + path
    return path

 
def read_T_P_DV(path,DV_filename,source):
    '''read_T_P_DV(path,DV_filename,source) stores temperature T, pressure P, and 
    some dependent variable DV as numpy arrays.
    '''
    if source == "CoolProp":
        T  = loadtxt(path + "CoolProp_temperatures.txt",
                        comments='#',delimiter=",",unpack=False)
        P  = loadtxt(path + "CoolProp_pressures.txt",
                        comments='#',delimiter=",",unpack=False)
        DV = loadtxt(path + DV_filename,
                        comments='#',delimiter=",",unpack=False,skiprows=1)
    elif source == "PelePhysics":
        T  = loadtxt(path + "PelePhys_temperatures.txt",
                        comments='#',delimiter=",",unpack=False)
        P  = loadtxt(path + "PelePhys_pressures.txt",
                        comments='#',delimiter=",",unpack=False)
        DV = loadtxt(path + DV_filename,
                        comments='#',delimiter=",",unpack=False,skiprows=1)
    else:
        print(" ERROR: Source of data cannot be identified! Check input data."); return
    return T,P,DV


def capture_row_max(array,x_labels,y_labels):
	num_rows = len(array)
	num_cols = len(array[0])

	# instantiate storage
	ind   = [0]*num_rows
	x     = [0]*num_rows
	y     = [0]*num_rows

	# capture indices as (x,y) for max along each row
	for i in range(num_rows):
		ind[i] = array[i].argmax()
		x[i]   = x_labels[ind[i]]
		y[i]   = y_labels[i]
	return x,y
		

######################################################################################## 




######################################################################################## 
# MISCELLANEOUS FUNCTIONS
######################################################################################## 

def approx_viscosity_error(path1,path2):
    '''approx_viscosity_error() computes |dmu/drho * density error|
    '''
	# Parse data identifiers
    path, DV_filename, species, prop, unit, source = read_titles(
                         path1, "CoolProp_density.txt", path2, "PelePhys_density.txt")

    # Compute density error
    density_error = compute_error(path[0],DV_filename[0],path[1],DV_filename[1])
    density_error = np.absolute(density_error)

    # Store T, P, dmudrho
    T, P, dmudrho = read_T_P_DV(path[1],"PelePhys_dmudrho.txt",source[1])
    dmudrho = np.absolute(dmudrho)

    # Compute approximate viscosity error
    result = np.multiply(density_error,dmudrho) 

	# Configure plot state
    plot_state = 1
    zoom_in_bound = 150*8;

    if plot_state == 0:
        result = result[0:zoom_in_bound,0:zoom_in_bound]
        fig_width   = 14
        fig_height  = 14/4.5	
        cbar_shrink = 0.72
        cbar_aspect = 7

    elif plot_state == 1:
        T = T[0:zoom_in_bound]
        P = P[0:zoom_in_bound]
        result = result[0:zoom_in_bound,0:zoom_in_bound]
        fig_width   = 8
        fig_height  = 7	
        cbar_shrink = 0.92
        cbar_aspect = 15

    # Plot pcolor data	
    x_axis = T
    y_axis = P
    X,Y = np.meshgrid(x_axis,y_axis)
    fig,ax = plt.subplots(1,1)
    cp = ax.pcolor(X,Y,result)
    title = species[0] + " (dmu/drho) * err(density)"
    plt.title(title)
    plt.xlabel("Temperature [K]")
    plt.ylabel("Pressure [atm]")
    ax.set_aspect('equal', adjustable='box')
    fig.set_size_inches(fig_width,fig_height)    
    fig.colorbar(cp,shrink=cbar_shrink,aspect=cbar_aspect)

	# Save result
    output_filename = "Error_Approximate_Viscosity"
    plt.savefig(output_filename,dpi=600)

    return T, P, result

   
def central_difference(path,DV_filename):
	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,DV_filename)

	# Determine Coordinates
	if   os.path.isfile(path + source + "_IV_pressure.txt"): state = "P"
	elif os.path.isfile(path + source + "_IV_density.txt"):  state = "D"
	else: print( "ERROR: Coordinate system (T,P) or (T,rho) cannot be determined.")

	# Load temperature, pressure/density, and DV
	x = loadtxt(path + source + "_IV_temperature.txt",delimiter=',')
	if   state == "P": y = loadtxt(path + source + "_IV_pressure.txt",delimiter=',')
	elif state == "D": y = loadtxt(path + source + "_IV_density.txt",delimiter=',')
	DV = loadtxt(path + DV_filename, delimiter=",",skiprows=1)

	print(x)
	
	# Compute partial derivatives (uses central differencing)
	delta_x = x[1] - x[0]
	delta_y = y[1] - y[0]
	grad = np.gradient(DV,delta_x,delta_y,edge_order=2)
	temp = grad[0] 		# reverse array orders because numpy takes d/dy first
	grad[0] = grad[1]
	grad[1] = temp 

	# Save derivatives as .txt file
	output_filename = ['','']
	short_filename  = ['','']
	short_filename[0] = source + '_DV_d' + prop + 'dT(numerical).txt'
	if   state == "P": short_filename[1] = source + '_DV_d' + prop + 'dP(numerical).txt'	
	elif state == "D": short_filename[1] = source + '_DV_d' + prop + 'drho(numerical).txt'		
	output_filename[0] = path + short_filename[0]
	output_filename[1] = path + short_filename[1]

	header = ['','']
	header[0] = species + ' d' + prop + 'dT ' + unit + ' ' + source
	if   state == "P": header[1] = species + ' d' + prop + 'dP '  + unit + ' ' + source
	elif state == "D": header[1] = species + ' d' + prop + 'drho ' + unit + ' ' + source

	np.savetxt(output_filename[0], grad[0], delimiter=',', header = header[0], comments = '')
	np.savetxt(output_filename[1], grad[1], delimiter=',', header = header[1], comments = '')

	# Plot derivatives
	if state == "P":
		PT_plot(path, short_filename[0])
		PT_plot(path, short_filename[1])
	elif state == "rho":
		rhoT_plot(path, short_filename[0])
		rhoT_plot(path, short_filename[1])
	
	plot_scalar(path,short_filename[0])
	plot_scalar(path,short_filename[1])
	
	return grad[0]

########################################################################################






########################################################################################
# deprecated functions 
######################################################################################## 

def PT_plot(path,DV_filename):
	"""PT_plot(path,DV_filename) plots a physical property stored in DV_filename
	as a function of temperature (x-axis) and pressure (y-axis).

	The first line of DV_filename must be formated as 'species property [units] source'
	"""
	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,DV_filename)

	# Read temperature, pressure, and dependent variable
	if os.path.isfile(path + "CoolProp_temperatures.txt"):
		T = loadtxt(path + "CoolProp_temperatures.txt", comments='#',delimiter=",",unpack=False)
		P = loadtxt(path + "CoolProp_pressures.txt",    comments='#',delimiter=",",unpack=False)
	elif os.path.isfile(path + "PelePhys_temperatures.txt"):
		T = loadtxt(path + "PelePhys_temperatures.txt", comments='#',delimiter=",",unpack=False)
		P = loadtxt(path + "PelePhys_pressures.txt",    comments='#',delimiter=",",unpack=False)
	DV = loadtxt(path + DV_filename, comments='#',delimiter=",",unpack=False, skiprows=1)

	plot_state = 0
	zoom_in_bound = 150*8;

	# Plot full range of T and P
	if plot_state == 0:
		# Plot pcolor data
		x_axis = T
		y_axis = P
		X,Y = np.meshgrid(x_axis,y_axis)
		fig,ax = plt.subplots(1,1)
		cp = ax.pcolor(X,Y,DV)
		title = species + " " + prop + " (" + source + ")"
		plt.title(title)
		plt.xlabel("Temperature [K]")
		plt.ylabel("Pressure [atm]")
		ax.set_aspect('equal', adjustable='box')
		fig.set_size_inches(14,14/4.5)    
		fig.colorbar(cp,shrink=0.72,aspect=7)

		# Plot maxima for specific heat
		if prop == "Specific Heat":
			x_max, y_max = capture_row_max(DV,T,P)
			plt.plot(x_max,y_max,"black",linewidth=1)
			plt.show()

	# Plot zoomed in near critical point
	elif plot_state == 1:
		# Plot pcolor data
		x_axis = T[0:zoom_in_bound]
		y_axis = P[0:zoom_in_bound]
		X,Y = np.meshgrid(x_axis,y_axis)
		fig,ax = plt.subplots(1,1)
		cp = ax.pcolor(X,Y,DV[0:zoom_in_bound,0:zoom_in_bound])
		title = species + " " + prop + " (" + source + ")"
		plt.title(title)
		plt.xlabel("Temperature [K]")
		plt.ylabel("Pressure [atm]")
		ax.set_aspect('equal', adjustable='box')
		fig.set_size_inches(8,7)   
		fig.colorbar(cp,shrink=0.92,aspect=15)

		# Plot maxima for specific heat
		if prop == "Specific Heat":
			x_max, y_max = capture_row_max(DV[0:zoom_in_bound,0:zoom_in_bound],
				T[0:zoom_in_bound],P[0:zoom_in_bound])
			plt.plot(x_max,y_max,"white",linewidth=1)
			plt.show()


	# Save result
	output_filename = '_'.join([source,species] + prop.split(" "))
	plt.savefig(output_filename,dpi=600)

	return


def rhoT_plot(path,DV_filename):
	"""PT_plot(path,DV_filename) plots a physical property stored in DV_filename
	as a function of temperature (x-axis) and pressure (y-axis).

	The first line of DV_filename must be formated as 'species property [units] source'
	"""
	# Parse data identifiers
	path = check_path(path)
	(species, prop, unit, source) = read_title(path,DV_filename)

	# Read temperature, pressure, and dependent variable
	if os.path.isfile(path + "CoolProp_temperatures.txt"):
		T   = loadtxt(path + "CoolProp_temperatures.txt", comments='#',delimiter=",",unpack=False)
		rho = loadtxt(path + "CoolProp_densities.txt",    comments='#',delimiter=",",unpack=False)
	elif os.path.isfile(path + "PelePhys_temperatures.txt"):
		T   = loadtxt(path + "PelePhys_temperatures.txt", comments='#',delimiter=",",unpack=False)
		rho = loadtxt(path + "PelePhys_densities.txt",    comments='#',delimiter=",",unpack=False)
	DV = loadtxt(path + DV_filename, comments='#',delimiter=",",unpack=False, skiprows=1)

	plot_state = 0
	zoom_in_bound = 150*8;

	# Plot full range of rho and T
	if plot_state == 0:
		# Plot pcolor data
		x_axis = T
		y_axis = rho
		X,Y = np.meshgrid(x_axis,y_axis)
		fig,ax = plt.subplots(1,1)
		cp = ax.pcolor(X,Y,DV)
		title = species + " " + prop + " (" + source + ")"
		plt.title(title)
		plt.xlabel("Temperature [K]")
		plt.ylabel("Density [kg m^-3]")
		ax.set_aspect('equal', adjustable='box')
		fig.set_size_inches(14,14/4.5)    
		fig.colorbar(cp,shrink=0.72,aspect=7)

		# Plot maxima for specific heat
		if prop == "Specific Heat":
			x_max, y_max = capture_row_max(DV,T,rho)
			plt.plot(x_max,y_max,"black",linewidth=1)
			plt.show()

	# Plot zoomed in near critical point
	elif plot_state == 1:
		# Plot pcolor data
		x_axis = T[0:zoom_in_bound]
		y_axis = P[0:zoom_in_bound]
		X,Y = np.meshgrid(x_axis,y_axis)
		fig,ax = plt.subplots(1,1)
		cp = ax.pcolor(X,Y,DV[0:zoom_in_bound,0:zoom_in_bound])
		title = species + " " + prop + " (" + source + ")"
		plt.title(title)
		plt.xlabel("Temperature [K]")
		plt.ylabel("Density [kg m^-3]")
		ax.set_aspect('equal', adjustable='box')
		fig.set_size_inches(8,7)   
		fig.colorbar(cp,shrink=0.92,aspect=15)

		# Plot maxima for specific heat
		if prop == "Specific Heat":
			x_max, y_max = capture_row_max(DV[0:zoom_in_bound,0:zoom_in_bound],
				T[0:zoom_in_bound],rho[0:zoom_in_bound])
			plt.plot(x_max,y_max,"white",linewidth=1)
			plt.show()


	# Save result
	output_filename = '_'.join([source,species] + prop.split(" "))
	plt.savefig(output_filename,dpi=600)

	return




def PT_plot_error(path1,DV_filename1,path2,DV_filename2):
	"""PT_plot_error(path1,DV_filename1,path2,DV_filename2) plots the percent error
	between PelePhysics and CoolProp in some physical property as a function of 
	temperature and pressure. 
	
	Error is reported as (PelePhys - CoolProp)/CoolProp so that positive errors
	indicate the PelePhysics value is greater than the corresponding CoolProp value,
	and vice versa for negative errors.
	"""
	# Parse data identifiers
	DV_filename = ['']*2; path = ['']*2
	species = ['']*2; prop = ['']*2; unit = ['']*2; source = ['']*2
	path[0] = path1; path[1] = path2;
	DV_filename[0] = DV_filename1; DV_filename[1] = DV_filename2;

	for i in range(2):
		path[i] = check_path(path[i])
		(species[i],prop[i],unit[i],source[i]) = read_title(path[i],DV_filename[i])	

	# Read temperature, pressure, and dependent variables
	T = ['','']; P = ['','']; DV = ['','']
	if source[0] == "CoolProp":
		T[0]  = loadtxt(path[0] + "CoolProp_IV_temperature.txt", 
						comments='#',delimiter=",",unpack=False)	
		P[0]  = loadtxt(path[0] + "CoolProp_IV_pressure.txt",    
						comments='#',delimiter=",",unpack=False)
		DV[0] = loadtxt(path[0] + DV_filename[0], 
						comments='#',delimiter=",",unpack=False,skiprows=1)
		P[1]  = loadtxt(path[1] + "PelePhys_IV_pressure.txt",    
						comments='#',delimiter=",",unpack=False)
		T[1]  = loadtxt(path[1] + "PelePhys_IV_temperature.txt", 
						comments='#',delimiter=",",unpack=False)
		DV[1] = loadtxt(path[1] + DV_filename[1], 
						comments='#',delimiter=",",unpack=False,skiprows=1)
	elif source[1] == "CoolProp":
		T[1]  = loadtxt(path[1] + "CoolProp_IV_temperature.txt", 
						comments='#',delimiter=",",unpack=False)	
		P[1]  = loadtxt(path[1] + "CoolProp_IV_pressure.txt",    
						comments='#',delimiter=",",unpack=False)
		DV[1] = loadtxt(path[1] + DV_filename[0], 
						comments='#',delimiter=",",unpack=False,skiprows=1)
		P[0]  = loadtxt(path[0] + "PelePhys_IV_pressure.txt",    
						comments='#',delimiter=",",unpack=False)
		T[0]  = loadtxt(path[0] + "PelePhys_IV_temperature.txt", 
						comments='#',delimiter=",",unpack=False)
		DV[0] = loadtxt(path[0] + DV_filename[1], 
						comments='#',delimiter=",",unpack=False,skiprows=1)
	else: 
		print(" ERROR: CoolProp data cannot be identified! Check input data."); return
		
	# Check for data mismatches
	if species[0] != species[1]:
		print(" ERROR: Chemical species mismatch! Check input data."); return
	if prop[0] != prop[1]:
		print(" ERROR: Physical property mismatch! Check input data."); return
	if len(T[0]) != len(T[1]):
		print(" ERROR: Number of temperature readings mismatched! Check input data."); return
	if len(P[0]) != len(P[1]):
		print(" ERROR: Number or pressure readings mismatched! Check input data."); return
	if max(abs(x) for x in T[0]-T[1]) > 0:
		print(" ERROR: Inconsistent temperature values! Check input data."); return
	if max(abs(x) for x in P[0]-P[1]) > 0:
		print(" ERROR: Inconsistent pressure values! Check input data."); return

    # Store error
	if source[0] == "CoolProp":
		error = (DV[1] - DV[0])/DV[0]
	elif source[1] == "CoolProp":
		error = (DV[0] - DV[1])/DV[1]
	else:
		print(" ERROR: CoolProp data cannot be identified! Check input data."); return

		#	error = np.absolute(error)


	plot_state = 0
	zoom_in_bound = 150*8;

	# Plot full range of T and P
	if plot_state == 0:
		x_axis = T[0]
		y_axis = P[0]
		X,Y = np.meshgrid(x_axis,y_axis)

		# Plot pcolor data
		fig,ax = plt.subplots(1,1)
		cp = ax.pcolor(X,Y,error)
		title = species[0] + " " + prop[0] + " Percent Error"
		plt.title(title)
		plt.xlabel("Temperature [K]")
		plt.ylabel("Pressure [atm]")
		ax.set_aspect('equal', adjustable='box')
		fig.set_size_inches(14,14/4.5)    
		fig.colorbar(cp,shrink=0.72,aspect=7)

		# Plot maxima for specific heat
		if prop[0] == "Specific Heat":
			x_max = ['']*2; y_max = ['']*2;
			for i in range(2):
				x_max[i], y_max[i] = capture_row_max(DV[i],T[i],P[i])
				plt.plot(x_max[i],y_max[i],"black",linewidth=1)
				plt.show()

	elif plot_state == 1:
		x_axis = T[0][0:zoom_in_bound]
		y_axis = P[0][0:zoom_in_bound]
		X,Y = np.meshgrid(x_axis,y_axis)

		# Plot pcolor data
		fig,ax = plt.subplots(1,1)
		cp = ax.pcolor(X,Y,error[0:zoom_in_bound,0:zoom_in_bound])
		title = species[0] + " " + prop[0] + " %Error"
		plt.title(title)
		plt.xlabel("Temperature [K]")
		plt.ylabel("Pressure [atm]")
		ax.set_aspect('equal', adjustable='box')
		fig.set_size_inches(8,7)    
		fig.colorbar(cp,shrink=0.92,aspect=15)

		# Plot maxima for specific heat
		if prop[0] == "Specific Heat":
			x_max = ['']*2; y_max = ['']*2;
			for i in range(2):
				x_max[i], y_max[i] = capture_row_max(DV[i],T[i],P[i])
				plt.plot(x_max[i],y_max[i],"black",linewidth=1)
				plt.show()



	# Save result
	output_filename = '_'.join(["Error",species[0]] + prop[0].split(" "))
	plt.savefig(output_filename,dpi=600)
	
	return


