#! /usr/bin/env python

## This module calculates statistical properties of a ramses output
from functools import partial
import multiprocessing as mp
from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis
from module_visualization import cal_rot_mat
from module_core import normalisation
import argparse
from pymses import RamsesOutput
from pymses.utils import constants as C
from pymses.utils.regions import Sphere, Cylinder
from pymses.filters import CellsToPoints, RegionFilter, PointFunctionFilter, PointRandomDecimatedFilter
import numpy as np
from pymses.analysis.point_sampling import PointSamplingProcessor
#import multiprocessing as mp
import os
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
import matplotlib.pyplot as plt
from module_visualization import save_fig
## module usage:
## boxlen, len_pc, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
#############################################################################################
def var_op(cells, var='density'):
    if var == 'density': return cells["rho"]
    elif var == 'temperature': return cells["P"]/cells["rho"]
    elif var == 'ekin': return cells["rho"]*cells["vel"]**2
#############################################################################################
def weight_op(cells, weight='volume'):
    dx = cells.get_sizes()
    vol = dx**3
    if weight == "mass": return cells["rho"] * vol
    elif weight == "volume": return vol
    else: return np.ones_like(vol)
#############################################################################################
def do_pdf(cells, bins, var='density', weight='volume'):
    data = var_op(cells, var=var)
    weights = weight_op(cells, weight=weight)
    pdf, range = np.histogram(data,bins=bins,weights=weights)
    return pdf    
#############################################################################################
def pdf_var(ro, lmax=20, scale='log', var='density', weight='volume', vmin=None, vmax=None, nbins=40, ncpu=1, filt_frac=0.01): # cpu-wise pdf for large output
    # var: 'density', 'temperature', 'ekin'
    # weight: 'volume', 'mass'
    if var=='density': svar=["rho"] 
    elif var=='temperature': svar=["rho","P"]
    elif var=='ekin': svar=["rho","vel"]
    amr = ro.amr_source(svar)
    cell_source = CellsToPoints(amr, smallest_cell_level=lmax)
    #filt_func = lambda dset: (var_op(dset)>vmin) * (var_op(dset)<vmax)
    #cell_source = PointFunctionFilter(filt_func, cell_source)
    if vmin is None or vmax is None: #undefined histogram range
        sub_cells = PointRandomDecimatedFilter(filt_frac, cell_source) 
        data = var_op(sub_cells.flatten(),var=var)
        if vmin is None: vmin = data.min() / 1.01
        if vmax is None: vmax = data.max() * 1.01
    if scale=='log': bin_edges = np.logspace(np.log10(vmin),np.log10(vmax),nbins+1)
    else: range = np.linspace(vmin,vmax,nbins+1)
    if ncpu<=1: 
        pool = None
        func = lambda x: do_pdf(x, bin_edges, var=var, weight=weight)
        get_f = lambda x: x
    else: 
        pool = mp.Pool(ncpu)
        func = lambda x: pool.apply_async(partial(do_pdf, bins=bin_edges, var=var, weight=weight),x)
        get_f = lambda x: x.get()
    pdf = np.zeros(nbins)
    sub_pdf = map(func, cell_source.iter_dsets())
    for sub in sub_pdf: pdf += get_f(sub)
    return pdf, bin_edges
#############################################################################################
def plot_hists(readdir, outputs, writedir=None, scl=20, scale='log',var='density',weight='volume', nbins=40, overwrite=False, ncpu=1,filt_frac=0.01,cumulated=False):
    if writedir is None: writedir=readdir
    figname = writedir+'/pdf_'+weight+'_'+var
    if os.path.exists(figname+'.pdf') and not overwrite: return
    ind = readdir[::-1].find('/')
    if ind>0: titlepath = readdir[len(readdir)-ind:]
    else: titlepath = readdir

    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
    len_AU = len_pc * 2.e5
    if var=='density': unit_bin = 1. #dens_gcc
    elif var=='temperature': unit_bin = temperature_K
    elif var=='ekin': unit_bin = mass_kg * vel_ms**2 / (len_cm/1.e2)**3 * 0.5
    if weight=='volume': unit_pdf = (len_cm/1.e2)**3
    elif weight=='mass': unit_pdf = mass_kg
    else: unit_pdf = 1.
    f_pdf = plt.figure()
    if cumulated: f_cdf = plt.figure()
    for ioutput in outputs:
        ro = RamsesOutput(readdir, ioutput, order='<')
        time = ro.info["time"]*time_Myr
    #ind = readdir[::-1].find('/')
    #if ind>0: titlepath = readdir[len(readdir)-ind:]
    #else: titlepath = readdir
    #titre = titlepath+' '+ str(time)[0:5] +'(Myr)' 
        pdf, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var=var, weight=weight, nbins=nbins, ncpu=ncpu, filt_frac=filt_frac) 
        bin_centers = np.sqrt(bin_edges[:-1]*bin_edges[1:])
        plt.figure(f_pdf.number)
        plt.loglog(bin_centers*unit_bin,pdf*unit_pdf,label='%i year'%(time*1.e6))
        if cumulated: 
            pdf = np.cumsum(pdf)
            plt.figure(f_cdf.number)
            plt.semilogx(bin_centers*unit_bin,pdf*unit_pdf,label='%i year'%(time*1.e6))
    #plt.ylim([pdf.min()*unit_pdf,pdf.max()*unit_pdf])
    #plt.ylim([1.e30, 3.e30])
    #plt.yticks(np.arange(1.e30, 4.e30, 1.e30))
    plt.figure(f_pdf.number)
    plt.legend(loc='best')
    plt.title(titlepath+' '+weight+'-'+var)
    if not os.path.exists(figname+'.pdf') or overwrite: save_fig(figname)
    if cumulated: 
        figname = figname+'_cum'
        plt.figure(f_cdf.number)
        plt.legend(loc='best')
        plt.title(titlepath+' '+weight+'-'+var)
        if not os.path.exists(figname+'.pdf') or overwrite: save_fig(figname)
    plt.show() 
#############################################################################################

#############################################################################################
## generat sampling point sets on cylinder surface
def cyl_samples(r_over_h,n_h,random=True):
    if random: #random sampling
        angles = 2.*np.pi*np.random.rand(int(n_h*n_h*2*np.pi*r_over_h))
        side_z = np.random.rand(int(n_h*n_h*2*np.pi*r_over_h))-0.5
        face_xy = np.random.rand(int(n_h*r_over_h*2)**2,2)*2.-1.
    else: #regular sampling
        angles = np.linspace(0,2.*np.pi,int(n_h*n_h*2*np.pi*r_over_h),endpoint=False)
        side_z = np.linspace(-0.5,0.5,int(n_h))
        angles, side_z = np.meshgrid(angles,side_z)
        angles = angles.flatten(); side_z = side_z.flatten()
        face_xy = np.linspace(-1,1,int(n_h*r_over_h*2))
        fx,fy = np.meshgrid(face_xy,face_xy)
        face_xy = np.vstack((fx.flatten(),fy.flatten())).T
    face_r = np.linalg.norm(face_xy,axis=1)
    ind = np.where(face_r<=1)[0]; nface = np.sum(ind)
    face_xy = face_xy[ind,:]; face_r = face_r[ind]
    side_xy = np.vstack((np.cos(angles), np.sin(angles))).T
    return side_xy, side_z, face_xy, face_r

#############################################################################################
## calculate flux in shells of cylinder
def cylinder_flux(ro,c_center,c_axis,c_radius,c_height,n_shell=3,n_h=20,random_sample=True,scl=10):
    amr = ro.amr_source(["rho", "vel"])
    amr2 = ro.amr_source(["rho", "vel"])
    psamp = PointSamplingProcessor(amr)
    side_xy, side_z, face_xy, face_r = cyl_samples(c_radius/c_height,n_h,random=random_sample)
    print 'sample size', side_xy.shape[0], face_xy.shape[0]
    r_shell = np.linspace(0,c_radius,n_shell+1)
    flux_r = np.zeros(n_shell)
    flux_z_in = np.zeros(n_shell)
    flux_z_out = np.zeros(n_shell)     

    c_axis = c_axis/np.linalg.norm(c_axis)
    azz, axx, ayy, mat_Rot = cal_rot_mat(c_axis)

    print 'c_axis', c_axis
    print 'mat_Rot', mat_Rot
    #return

    vec_rad = np.hstack((side_xy, np.zeros((side_xy.shape[0],1))))
    vec_rad = np.sum(vec_rad[:,np.newaxis,:]*mat_Rot[np.newaxis,:,:],axis=2) # radial vector for flux calculation
    vec_face = np.hstack((face_xy, np.zeros((face_xy.shape[0],1))))
    vec_face = np.sum(vec_face[:,np.newaxis,:]*mat_Rot[np.newaxis,:,:],axis=2)

    vec_z = side_z[:,np.newaxis] * c_axis * c_height   
    face_surf =  np.pi * (r_shell[1:]**2-r_shell[:-1]**2)
    face_r = face_r * c_radius

    ## vertical flux at faces
    scl_f = min(int(np.log(c_height/n_h)/np.log(0.5)),scl)
    print 'max read level face =', scl_f
    for top_bot in [-1, 1]:
        face = c_center + vec_face * c_radius + top_bot * 0.5 * c_axis * c_height
        face_points = psamp.process(face,add_level=False,max_search_level=scl_f)
        rho_vinf = face_points["rho"]*top_bot*np.sum(face_points["vel"]*c_axis[np.newaxis,:],axis=1)
        print rho_vinf
        out = (rho_vinf>0)
        flux_out, bins = np.histogram(face_r, bins = r_shell, weights=rho_vinf*out)
        flux_in, bins = np.histogram(face_r, bins = r_shell, weights=rho_vinf*(~out))
        flux_points, bins = np.histogram(face_r, bins = r_shell)
        flux_z_out += flux_out / flux_points
        flux_z_in += flux_in / flux_points  

    ## mass inside cylindrical shells
    scl_f = max(int(np.log(c_radius/n_shell)/np.log(0.5)),scl)+2
    print 'max read level mass =', scl_f
    cyl = Cylinder(c_center, c_axis, c_radius, c_height)
    amr2 = RegionFilter(cyl,amr2)         
    cell_source = CellsToPoints(amr2, smallest_cell_level=scl_f)
    cells = cell_source.flatten()
    dm = cells["rho"]*cells.get_sizes()**3
    print 'total mass', np.sum(dm)
    dx2 = np.sum((cells.points-c_center)**2,axis=1)
    dz2 = np.sum((cells.points-c_center)*c_axis[np.newaxis,:],axis=1)**2
    r_shell_mass = r_shell.copy(); r_shell_mass[-1] = r_shell_mass[-1]+0.5**(scl_f)# *2-r_shell_mass[-2]
    M_shell, bins = np.histogram(np.sqrt(dx2-dz2), bins = r_shell_mass, weights=dm)
    print 'total mass shell', np.sum(M_shell)

    rad_surf = 2. * np.pi * r_shell[1:] * c_height
    ## radial flux at shells
    for ir in range(n_shell):
        r = r_shell[ir+1]
        scl_f = min(int(np.log(c_height*r/c_radius/n_h)/np.log(0.5)),scl)
        print 'max read level shell=', scl_f
        side = c_center + vec_rad*r + vec_z
        side_points = psamp.process(side,add_level=False,max_search_level=scl_f)
        rho_vinf = np.mean(side_points["rho"]*np.sum(side_points["vel"]*vec_rad,axis=1))
        flux_r[ir] = rho_vinf  

    return flux_r, flux_z_out, flux_z_in, r_shell, M_shell, face_surf, rad_surf  

#############################################################################################
def plot_flux(readdir,ioutput,writedir=None,fig_ps=False,c_radius=100,c_hieghts=[50,30,20,1],scl=18,ns=24,nh=12,overwrite=False):
    if writedir is None: writedir=readdir
    figname = writedir+'/mass_flux_'+str(ioutput).zfill(5)      
    if os.path.exists(figname+'.pdf') and not overwrite: return

    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
    len_AU = len_pc * 2.e5
    unit_accretion = mass_sol/time_Myr
    unit_flux = dens_gcc * vel_ms * 1.e3 # kg/m^2/s

    ro = RamsesOutput(readdir, ioutput, order='<')
    time = ro.info["time"]*time_Myr
    ind = readdir[::-1].find('/')
    if ind>0: titlepath = readdir[len(readdir)-ind:]
    else: titlepath = readdir
    titre = titlepath+' '+ str(time)[0:5] +'(Myr)'
    c_radius = c_radius/len_AU
    c_heights = [h/len_AU for h in c_hieghts]
    f, axes = plt.subplots(len(c_heights),2, sharex=True, sharey=False, figsize=(6,6))
    f.suptitle(titre)
    n_shell = min(ns, int(np.floor(c_radius * 2.**scl/10))) #at least 10 cells per shell, at most ns shells
    c_center = find_baricenter(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=200./len_AU,scl=lmax)
    c_axis, c_vel = find_axis(ro, c_center=c_center,c_radius=100./len_AU,scl=lmax)
    for i_height in range(len(c_heights)):
        c_height = c_heights[i_height]
        n_h =  min(nh,int(np.ceil(c_height * 2**scl))*2) #between 2 and nh points
        flux_r, flux_z_out, flux_z_in, r_shell, M_shell, face_surf, rad_surf = cylinder_flux(ro,c_center=c_center,c_axis=c_axis,c_radius=c_radius,c_height=c_height,n_shell=n_shell,n_h=n_h,random_sample=False,scl=scl)
        ## plot flux
        ax = axes[i_height,0]
        print 'lmax', lmax, 'lread', scl
        ax.plot(r_shell[1:]*len_pc*2.e5, flux_r*unit_flux,'b',lw=2,label='radial flux')
        ax.plot(r_shell[1:]*len_pc*2.e5, flux_z_out*unit_flux,'r--',label='vertical outflow/inflow')
        ax.plot(r_shell[1:]*len_pc*2.e5, flux_z_in*unit_flux,'r--')
        ax.plot(r_shell[1:]*len_pc*2.e5, (flux_z_out+flux_z_in)*unit_flux,'r',lw=2,label='vertical flux')
        ylim = 1.e-7
        ax.set_ylim([-ylim,ylim])
        ax.set_xlim([0,c_radius*len_pc*2.e5])
        ax.set_yscale('symlog',linthreshy=ylim/1.e3)
        ax.grid(True)

        ## plot mass accretion rate
        ax = axes[i_height,1]
        flux_z_out = np.cumsum(flux_z_out*face_surf)
        flux_z_in = np.cumsum(flux_z_in*face_surf)
        ax.plot(r_shell[1:]*len_pc*2.e5, flux_r*rad_surf*unit_accretion,'b',lw=2)#,label='radial accretion rate')
        ax.plot(r_shell[1:]*len_pc*2.e5, flux_z_out*unit_accretion,'r--')#,label='vertical outflow/inflow')
        ax.plot(r_shell[1:]*len_pc*2.e5, flux_z_in*unit_accretion,'r--')
        ax.plot(r_shell[1:]*len_pc*2.e5, (flux_z_out+flux_z_in)*unit_accretion,'r',lw=2)#,label='vertical accretion rate')
        ax.plot(r_shell[1:]*len_pc*2.e5, (flux_z_out+flux_z_in+flux_r*rad_surf)*unit_accretion,'gray',lw=2,label='total accretion rate')
        ax.plot(r_shell[1:]*len_pc*2.e5, np.cumsum(M_shell)*mass_sol,'y',lw=2,label='mass (M$_\odot$)')
        ylim = 1.e2
        ax.set_ylim([-ylim,ylim])
        ax.text(5,ylim/10., '2H = '+str(int(c_height*len_pc*2.e5))+'AU')
        ax.set_xlim([0,c_radius*len_pc*2.e5])
        ax.set_yscale('symlog',linthreshy=ylim/1.e3)
        ax.grid(True)
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
            #plt.ylabel('mass accretion rate (M$_\odot$/Myr)')

    f.text(0.01, 0.5, 'mass flux (kg m$^{-2}$ s$^{-1}$)', va='center', rotation='vertical')
    f.text(0.96, 0.5, 'mass accretion rate (M$_\odot$ Myr$^{-1}$)', va='center', rotation=270)
    for i in range(2): axes[-1,i].set_xlabel('radius (AU)')
    axes[0,0].text(5,ylim/10.,'axis '+str([round(c,2) for c in c_axis]))
    axes[0,0].text(5,ylim/30.,'center '+str([round(c*len_pc,4) for c in c_center]))
    axes[0,0].legend(loc='best',frameon=False,labelspacing=0.2)
    axes[0,1].legend(loc='best',frameon=False,labelspacing=0.2)
# Fine-tune figure; make subplots close to each other and hide x ticks for all but bottom plot.
    f.subplots_adjust(0.15,0.09,0.85,0.93,0,0)
    for ax in axes[:-1,:].flatten():
        for tic in ax.xaxis.get_major_ticks(): tic.tick1On = tic.tick2On = False
    #plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    save_fig(figname,ps=fig_ps,tight=False)
        #plt.show()


#############################################################################################


#############################################################################################


#############################################################################################

