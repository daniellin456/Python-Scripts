#! /usr/bin/env python
## This module contains basic analyses of core disc


##### pymses dependences ##################################################################
from pymses import RamsesOutput
from pymses.utils import constants as C
from pymses.utils.regions import Sphere, Cylinder
from pymses.filters import CellsToPoints, RegionFilter
from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
##### python dependences ##################################################################
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle as pck
##### module dependences ##################################################################
from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis, read_sink_cvs
from module_visualization import cal_rot_mat, save_fig
#from module_core import normalisation
## module usage:
## boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
#############################################################################################
#### Save pickle file
def save_pickle(filename,vars):
    f = open(filename,'w')
    for var in vars: pck.dump(var,f)
    f.close()
#############################################################################################
#### Read pickle file
def read_pickle(filename,nvar):
    f = open(filename,'rb')
    vars = []
    for i in range(nvar): vars.append(pck.load(f))
    f.close()
    return vars
#############################################################################################
## generate sampling point sets on cylinder surface
## r_over_h: radius/height; n_h: number of points along cynlider height
## side_xy & size_z on th lateral surface, face_xy & face_r on the spherical surface, 
## radius normalized to 1, z between [-0.5,0.5]
def cyl_samples(r_over_h,n_h,random=False):
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
def cylinder_flux(ro,c_center,c_axis,c_vel,c_radius,c_height,n_shell=3,n_h=20,random_sample=False,scl=20,add_j=False):
    amr = ro.amr_source(["rho", "vel"])
    amr2 = ro.amr_source(["rho"])
    psamp = PointSamplingProcessor(amr)
    side_xy, side_z, face_xy, face_r = cyl_samples(c_radius/c_height,n_h,random=random_sample)
    print 'sample size', side_xy.shape[0], face_xy.shape[0]
    r_shell = np.linspace(0,c_radius,n_shell+1)
    flux_r = np.zeros(n_shell)
    flux_z_in = np.zeros(n_shell)
    flux_z_out = np.zeros(n_shell) 
    if add_j:
        jflux_r = np.zeros(n_shell)
        jflux_z_in = np.zeros(n_shell)
        jflux_z_out = np.zeros(n_shell)    

    c_axis = c_axis/np.linalg.norm(c_axis)
    azz, axx, ayy, mat_Rot = cal_rot_mat(c_axis)

    print 'c_axis', c_axis
    print 'mat_Rot', mat_Rot

    vec_rad = np.hstack((side_xy, np.zeros((side_xy.shape[0],1))))
    vec_rad = np.sum(vec_rad[:,np.newaxis,:]*mat_Rot[np.newaxis,:,:],axis=2) # rotated radial vector for flux calculation
    vec_face = np.hstack((face_xy, np.zeros((face_xy.shape[0],1))))
    vec_face = np.sum(vec_face[:,np.newaxis,:]*mat_Rot[np.newaxis,:,:],axis=2) # rotated face vector for flux calculation

    vec_z = side_z[:,np.newaxis] * c_axis * c_height   
    face_surf =  np.pi * np.diff(r_shell**2)
    face_r = face_r * c_radius

    ## vertical flux at faces
    scl_f = min(int(np.log(c_height/n_h)/np.log(0.5)),scl)
    print 'max read level face =', scl_f
    for top_bot in [-1, 1]:
        face = c_center + vec_face * c_radius + top_bot * 0.5 * c_axis * c_height
        face_points = psamp.process(face,add_level=False,max_search_level=scl_f)
        rho_vinf = face_points["rho"]*top_bot*np.sum((face_points["vel"]-c_vel[np.newaxis,:])*c_axis[np.newaxis,:],axis=1)
        out = (rho_vinf>0)
        flux_out, bins = np.histogram(face_r, bins = r_shell, weights=rho_vinf*out)
        flux_in, bins = np.histogram(face_r, bins = r_shell, weights=rho_vinf*(~out))
        flux_points, bins = np.histogram(face_r, bins = r_shell)
        flux_z_out += flux_out / flux_points
        flux_z_in += flux_in / flux_points  
        if add_j:
            rvphi = np.linalg.norm(np.cross( vec_face, (face_points["vel"]-c_vel[np.newaxis,:])), axis=1) * c_radius
            j_vinf = rho_vinf * rvphi
            jflux_out, bins = np.histogram(face_r, bins = r_shell, weights=j_vinf*out)
            jflux_in, bins = np.histogram(face_r, bins = r_shell, weights=j_vinf*(~out))
            jflux_z_out += jflux_out / flux_points
            jflux_z_in += jflux_in / flux_points

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
        rho_vinf = side_points["rho"]*np.sum((side_points["vel"]-c_vel[np.newaxis,:])*vec_rad,axis=1)
        flux_r[ir] = np.mean(rho_vinf)  
        if add_j:
            rvphi = np.linalg.norm(np.cross( vec_rad, (side_points["vel"]-c_vel[np.newaxis,:]) ), axis=1) * r
            j_vinf = rho_vinf * rvphi
            jflux_r[ir] = np.mean(j_vinf)

    if add_j: return flux_r, flux_z_out, flux_z_in, jflux_r, jflux_z_out, jflux_z_in, r_shell, M_shell, face_surf, rad_surf
    else: return flux_r, flux_z_out, flux_z_in, r_shell, M_shell, face_surf, rad_surf  

#############################################################################################
def plot_flux(readdir,ioutput,writedir=None,fig_ps=False,c_radius=80,c_hieghts=[50,30,20,1],scl=18,ns=24,nh=12,overwrite=False):
    if writedir is None: writedir=readdir
    figname = writedir+'/mass_flux_'+str(ioutput).zfill(5)      
    if os.path.exists(figname+'.pdf') and not overwrite: return None, None, None

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
    sinks = read_sink_cvs(ioutput,readdir,deli=',') 
    if sinks.shape[0]>=1: 
        Msink = np.sum(sinks[:,1])
        c_center = sinks[np.argmax(sinks[:,1]),3:6]/boxlen
    else:
        Msink = 0. 
        c_center = find_baricenter(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=200./len_AU,scl=lmax)
    print 'center = ',  c_center
    c_axis, c_vel = find_axis(ro, c_center=c_center,c_radius=20./len_AU,scl=lmax)
    for i_height in range(len(c_heights)):
        c_height = c_heights[i_height]
        n_h =  min(nh,int(np.ceil(c_height * 2**scl))*2) #between 2 and nh points
        flux_r, flux_z_out, flux_z_in, r_shell, M_shell, face_surf, rad_surf = cylinder_flux(ro,c_center=c_center,c_axis=c_axis,c_vel=c_vel,c_radius=c_radius,c_height=c_height,n_shell=n_shell,n_h=n_h,random_sample=False,scl=scl)
        ## plot flux
        ax = axes[i_height,0]
        print 'lmax', lmax, 'lread', scl
        ax.plot(r_shell[1:]*len_AU, flux_r*unit_flux,'b',lw=2,label='radial flux')
        ax.plot(r_shell[1:]*len_AU, flux_z_out*unit_flux,'r--',label='vertical outflow/inflow')
        ax.plot(r_shell[1:]*len_AU, flux_z_in*unit_flux,'r--')
        ax.plot(r_shell[1:]*len_AU, (flux_z_out+flux_z_in)*unit_flux,'r',lw=2,label='vertical flux')
        ylim = 1.e-8
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
        ax.plot(r_shell[1:]*len_pc*2.e5, np.cumsum(M_shell)*mass_sol+Msink,'y',lw=2,label='mass (M$_\odot$)')
        ylim = 1.e1
        ax.set_ylim([-ylim,ylim])
        ax.text(5,ylim/10., '2H = '+str(int(c_height*len_AU))+'AU')
        ax.set_xlim([0,c_radius*len_AU])
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
    return c_center, c_axis, c_vel

#############################################################################################


#############################################################################################


#############################################################################################

