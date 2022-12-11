#! /usr/bin/env python

## The module disc contains analyses of disc properties

##### pymses dependences ##################################################################
import pymses
from pymses.filters import PointFunctionFilter
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
from scipy.interpolate import interp1d, interp2d
from scipy.misc import derivative
from scipy.optimize import curve_fit
import scipy as sp
import pdb
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm, SymLogNorm,LogNorm
from operator import truediv
import glob
#import multiprocessing as mp
##### module dependences ##################################################################
from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis, normalisation, read_sink_cvs, norm_units
from module_visualization import cal_rot_mat, save_fig
from module_analysis import cylinder_flux, save_pickle, read_pickle
#from module_core import normalisation
## module usage:
## boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 

#############################################################################################
#### save file rho_func_ioutput.pkl containing CDFs of density 
def calc_rho_relations(path,num,path_out=None,overwrite=True,order='<',scl=20,region=None):
    if path_out==None: path_out=path
    filename = path_out+'/rho_func_'+str(num).zfill(5)+'.pkl'
    if os.path.exists(filename) and not overwrite: return
    ro = pymses.RamsesOutput(path,num,order=order)
    amr = ro.amr_source(["rho"])
    if region is not None: amr = RegionFilter(region,amr)
    cell_source = CellsToPoints(amr, smallest_cell_level=scl)
    cells = cell_source.flatten()
    dm = cells["rho"]*cells.get_sizes()**3
    sort_d = np.argsort(cells["rho"])[-1::-1]
    M_of_rho = np.cumsum( dm[sort_d] )
    save_pickle(filename,[cells["rho"][sort_d],M_of_rho])#,J_of_rho,r[sort_r],M_of_r])
#############################################################################################
#### read file rho_func_ioutput.pkl and produce interpoaltion functions of density
def rho_functions(filename):
    [rho, M_of_rho] = read_pickle(filename,2)
    func_M_rho = interp1d(rho,M_of_rho)
    func_dm_drho = lambda x: derivative(func_M_rho,x,dx=x*1e-1)
    func_d2m_drho2 = lambda x: derivative(func_M_rho,x,dx=x*1e-1,n=2)
    return func_M_rho, func_dm_drho, func_d2m_drho2, [rho[0],rho[-1]]
#############################################################################################
### calculate for several time steps and save in disc_basic_parameters.npz the parameters that define the disc
def disc_basic_parameters(path,path_out=None,overwrite=True,order='<',scl=20,center_def='None',rotate=None,c_axis=None,c_vel=None,center_img=[0.5,0.5,0.5],soutputs=False):
    if path_out == None: path_out = path
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    outputs = search_ro(path)
    nout = len(outputs)
    fileparas = path_out+'/disc_basic_parameters.npz'
    if os.path.exists(fileparas):
        discload = np.load(fileparas)
        ioutput_list = discload['ioutput_list'].tolist()
        time_list = discload['time_list'].tolist()
        Msink_list = discload['Msink_list'].tolist()
        Mdisc_list = discload['Mdisc_list'].tolist()
        thres_rho_list = discload['thres_rho_list'].tolist()
        center_sink_list = discload['center_sink_list'].tolist()
        center_disc_list = discload['center_disc_list'].tolist()
        cvel_sink_list = discload['cvel_sink_list'].tolist()
        cvel_disc_list = discload['cvel_disc_list'].tolist()
        ax_csink_list = discload['ax_csink_list'].tolist()
        ax_ctot_list = discload['ax_ctot_list'].tolist()
        r_csink_list = discload['r_csink_list'].tolist()
        r_ctot_list = discload['r_ctot_list'].tolist()
        print ioutput_list
    else:
        ioutput_list = []
        time_list = []
        Msink_list = []
        Mdisc_list = []
        thres_rho_list = []
        center_sink_list = []
        center_disc_list = []
        cvel_sink_list = []
        cvel_disc_list = []
        ax_csink_list = []
        ax_ctot_list = []
        r_csink_list = []
        r_ctot_list = []
    nrho = 100
    print nout, outputs
#    for iout in range(nout)[2::10]:
#        ioutput = outputs[iout]
    nsam = nout/6
    if not soutputs: soutputs = outputs[nsam::nsam]
    for ioutput in soutputs:
        if ioutput in ioutput_list: continue
        filename = path_out+'/rho_func_'+str(ioutput).zfill(5)+'.pkl'
        #if not os.path.exists(filename):  continue
        if not os.path.exists(filename): calc_rho_relations(path,ioutput,path_out=path_out,overwrite=False,order='<',scl=scl) 
        ioutput_list.append(ioutput)
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        time_list.append(ro.info['time'] * time_Myr)        
        sinks = read_sink_cvs(ioutput,path,deli=',')
        if sinks.shape[0]>=1:
            ind_star = np.argmax(sinks[:,1])
            Msink = sinks[ind_star,1]
            center_sink = sinks[ind_star,3:6]/boxlen
            cvel_sink = sinks[ind_star,6:9]
        else: Msink = 0.; center_sink = np.zeros(3); cvel_sink = np.zeros(3) 
        Msink_list.append(Msink)
        center_sink_list.append(center_sink)
        cvel_sink_list.append(cvel_sink)
        func_M_rho,func_dm_drho, func_d2m_drho2,rho_lims = rho_functions(filename)
        rho_list = np.logspace(np.log10(rho_lims[0])-1.5, np.log10(rho_lims[1])+1, nrho)
        m_list = func_M_rho(rho_list)
        dm_list = func_dm_drho(rho_list)
        d2m_list = func_d2m_drho2(rho_list)
        dlogm_list = dm_list*rho_list/m_list
        d2logm_list = d2m_list*rho_list**2/m_list + dm_list*rho_list/m_list - (dm_list*rho_list/m_list)**2
        ind_rho = np.argmax(dlogm_list)
        thres_rho = rho_list[ind_rho]
        thres_rho_list.append(thres_rho)
        filt_dens = lambda dset: (dset["rho"] >= thres_rho)
        amr = ro.amr_source(["rho","vel"])
        cell_source = CellsToPoints(amr, smallest_cell_level=scl)
        cell_source = PointFunctionFilter(filt_dens,cell_source)
        cells = cell_source.flatten()
        dm = cells["rho"]*cells.get_sizes()**3
        Mdisc = np.sum(dm)
        Mdisc_list.append(Mdisc*mass_sol)
        center_disc = np.sum(cells.points * dm[:,np.newaxis], axis=0) / Mdisc
        center_disc_list.append(center_disc)
        center_tot = (center_disc * Mdisc + center_sink * Msink/mass_sol) / (Mdisc+Msink/mass_sol)
        print 'center offset', np.linalg.norm(center_sink-center_disc)*len_AU, np.linalg.norm(center_sink-center_tot)*len_AU
        print 'Msink', Msink, 'Mdisc', Mdisc*mass_sol
        cvel_disc = np.sum(dm[:,np.newaxis]*cells["vel"],axis=0)/Mdisc
        cvel_tot = (cvel_disc * Mdisc + cvel_sink * Msink/mass_sol) / (Mdisc+Msink/mass_sol)
        cvel_disc_list.append(cvel_disc)
        if Msink!= 0.:
            positions = cells.points-center_sink[np.newaxis,:]
            velocities = cells["vel"]-cvel_sink[np.newaxis,:]            
            angmom = np.sum(dm[:,np.newaxis]*np.cross(positions,velocities),axis=0)
            axis = angmom / np.linalg.norm(angmom)
            radii = np.sqrt(np.sum(positions**2,axis=1) - np.sum(positions*axis[np.newaxis,:],axis=1)**2)
            rdisc = radii.max()
        else: 
            axis = np.zeros(3); rdisc = 0.        
        ax_csink_list.append(axis)
        r_csink_list.append(rdisc)

        positions = cells.points-center_tot[np.newaxis,:]
        velocities = cells["vel"]-cvel_tot[np.newaxis,:] 
        angmom = np.sum(dm[:,np.newaxis]*np.cross(positions,velocities),axis=0)+Msink/mass_sol*np.cross((center_sink-center_tot),(cvel_sink-cvel_tot))
        axis = angmom / np.linalg.norm(angmom)
        radii = np.sqrt(np.sum(positions**2,axis=1) - np.sum(positions*axis[np.newaxis,:],axis=1)**2)
        rdisc = radii.max()
        ax_ctot_list.append(axis)
        r_ctot_list.append(rdisc)
    isort = np.argsort(ioutput_list)
    reorder = lambda x: [x[i] for i in isort]
    time_list = reorder(time_list)
    Msink_list = reorder(Msink_list)
    Mdisc_list = reorder(Mdisc_list)
    thres_rho_list = reorder(thres_rho_list)
    center_sink_list = reorder(center_sink_list)
    center_disc_list = reorder(center_disc_list)
    cvel_sink_list = reorder(cvel_sink_list)
    cvel_disc_list = reorder(cvel_disc_list)
    ax_csink_list = reorder(ax_csink_list)
    ax_ctot_list = reorder(ax_ctot_list)
    r_csink_list = reorder(r_csink_list)
    r_ctot_list = reorder(r_ctot_list)
    ioutput_list = reorder(ioutput_list)
 
    np.savez(path_out+'/disc_basic_parameters', ioutput_list=ioutput_list, time_list=time_list, Msink_list=Msink_list, Mdisc_list=Mdisc_list, thres_rho_list=thres_rho_list, center_sink_list=center_sink_list, center_disc_list=center_disc_list, cvel_sink_list=cvel_sink_list, cvel_disc_list=cvel_disc_list, ax_csink_list=ax_csink_list, ax_ctot_list=ax_ctot_list, r_csink_list=r_csink_list, r_ctot_list=r_ctot_list)

#############################################################################################
### identify the cylinder that contains the disc, store in disc_prop.pkl
def disc_cylinder(path,path_out=None,overwrite=True,order='<',scl=20,center_def='None',rotate=None,c_axis=None,c_vel=None,center_img=[0.5,0.5,0.5]):
    if path_out == None: path_out = path
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    unit_flux = dens_gcc * vel_ms * 1.e3
    unit_jflux = unit_flux * vel_ms * len_cm * 1.e-2
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    outputs = search_ro(path)
    nout = len(outputs)
    f, ax = plt.subplots(6,1)
    f.suptitle(titlepath)
    nr = 100
    nr2 = 10
    time_list = []
    rr_list = []
    zz_list = []
    Msink_list = []
    Mdisc_list = []
    Vfluxout_list = []
    Vfluxin_list = []
    Rflux_list = []
    jVfluxout_list = []
    jVfluxin_list = []
    jRflux_list = []
    Sigma_list = []
    center_list = []
    cvel_list = []
    axis_list = []
    ioutput_list = []
    thres_rho_list = []
    center_sink_list = []
    center_disc_list = []
    cvel_sink_list = []
    rad_surf_list = []
    for iout in range(nout)[0::100]:
        ioutput = outputs[iout]
        filename = path_out+'/rho_func_'+str(ioutput).zfill(5)+'.pkl'
        if not os.path.exists(filename):  continue
        ioutput_list.append(ioutput)
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        time=ro.info['time'] * time_Myr
        sinks = read_sink_cvs(ioutput,path,deli=',')
        if sinks.shape[0]>=1:
            Msink = np.sum(sinks[:,1])
            c_center = sinks[np.argmax(sinks[:,1]),3:6]/boxlen
            vsink = sinks[np.argmax(sinks[:,1]),6:9]
        else:
            Msink = 0.; vsink = np.zeros(3)
            c_center = find_baricenter(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=200./len_AU,scl=lmax)
        print 'center = ',  c_center
        sink_center = c_center
        center_sink_list.append(sink_center)
        cvel_sink_list.append(vsink)
        #func_M_rho, func_J_rho, func_j_rho, func_dm_drho, func_d2m_drho2, func_M_r, func_rho_r, rho_lims, r_lims = rho_functions(filename)
        func_M_rho,func_dm_drho, func_d2m_drho2,rho_lims = rho_functions(filename)
        rho_list = np.logspace(np.log10(rho_lims[0])-1.5, np.log10(rho_lims[1])+1, nr)
        m_list = func_M_rho(rho_list)
        dm_list = func_dm_drho(rho_list)  
        d2m_list = func_d2m_drho2(rho_list)
        dlogm_list = dm_list*rho_list/m_list
        d2logm_list = d2m_list*rho_list**2/m_list + dm_list*rho_list/m_list - (dm_list*rho_list/m_list)**2
        ind_rho = np.argmax(dlogm_list)
        thres_rho = rho_list[ind_rho]
        thres_rho_list.append(thres_rho)
        ax[1].loglog(rho_list,m_list,label=time)
        ax[1].scatter(rho_list[ind_rho],m_list[ind_rho],s=20)
        #plt.semilogx(rho_list,dm_list*rho_list/m_list,label=time) 
        #plt.semilogx(rho_list,d2m_list*rho_list**2/m_list + dm_list*rho_list/m_list - (dm_list*rho_list/m_list)**2,label='2')   
        filt_dens = lambda dset: (dset["rho"] >= thres_rho*5)
        amr = ro.amr_source(["rho","vel"])
        cell_source = CellsToPoints(amr, smallest_cell_level=scl)
        cell_source = PointFunctionFilter(filt_dens,cell_source)
        cells = cell_source.flatten()
        dm = cells["rho"]*cells.get_sizes()**3
        mdisc = np.sum(dm)
        disc_center = np.sum(cells.points * dm[:,np.newaxis], axis=0) / mdisc
        center_disc_list.append(disc_center)
        tot_center = (disc_center * mdisc + c_center * Msink/mass_sol) / (mdisc+Msink/mass_sol)
        print 'center offset', np.linalg.norm(c_center-disc_center)*len_AU, np.linalg.norm(c_center-tot_center)*len_AU
        print 'Msink', Msink, 'Mdisc', mdisc*mass_sol
        c_center = sink_center
        c_axis, c_vel = find_axis(ro, c_center=c_center,c_radius=20./len_AU,scl=lmax)
        zz = np.sum((cells.points-c_center[np.newaxis,:]) * c_axis[np.newaxis,:], axis=1)/np.linalg.norm(c_axis)
        rr = np.sqrt( np.sum((cells.points-c_center[np.newaxis,:])**2,axis=1) - zz**2)
        ncell = len(rr)
        r_cyl = np.sort(rr)[-ncell/100]*len_AU
        z_cyl = np.sort(abs(zz))[-ncell/100]*len_AU
        print np.max(rr)*len_AU, r_cyl, np.max(abs(zz))*len_AU, z_cyl
        dm = cells["rho"]*cells.get_sizes()**3
        mdisc = np.sum(dm)
        c_vel = np.sum(dm[:,np.newaxis]*cells["vel"],axis=0)/mdisc
        cvel_list.append(c_vel)
        c_vel = (c_vel * mdisc + vsink * Msink) / (mdisc+Msink)
        c_vel = vsink
        c_axis = np.cross((cells.points-c_center),(cells["vel"]-c_vel))
        c_axis = np.sum(dm[:,np.newaxis]*c_axis,axis=0)
        c_axis = c_axis / np.linalg.norm(c_axis)
        center_list.append(c_center)
#        cvel_list.append(c_vel)
        axis_list.append(c_axis)
        zz = np.sum((cells.points-c_center[np.newaxis,:]) * c_axis[np.newaxis,:], axis=1)
        rr = np.sqrt( np.sum((cells.points-c_center[np.newaxis,:])**2,axis=1) - zz**2)
        ax[0].scatter(rr*len_AU,zz*len_AU,s=0.5)
        ncell = len(rr) 
        r_cyl = np.sort(rr)[-ncell/100]*len_AU
        z_cyl = np.sort(abs(zz))[-ncell/100]*len_AU
        print np.max(rr)*len_AU, r_cyl, np.max(abs(zz))*len_AU, z_cyl

        ax[0].plot([0,r_cyl,r_cyl,0],[z_cyl,z_cyl,-z_cyl,-z_cyl])
        time_list.append(time*1.e6)
        rr_list.append(r_cyl)
        zz_list.append(z_cyl)

        flux_r, flux_z_out, flux_z_in, jflux_r, jflux_z_out, jflux_z_in, r_shell, M_shell, face_surf, rad_surf = cylinder_flux(ro,c_center=c_center,c_axis=c_axis,c_vel=c_vel,c_radius=r_cyl/len_AU,c_height=z_cyl/len_AU*2,n_shell=nr2,n_h=8,random_sample=False,scl=scl,add_j=True)
        Rflux_list.append(flux_r*unit_flux)
        Vfluxout_list.append(flux_z_out*unit_flux)
        Vfluxin_list.append(flux_z_in*unit_flux)
        jRflux_list.append(jflux_r*unit_jflux)
        jVfluxout_list.append(jflux_z_out*unit_jflux)
        jVfluxin_list.append(jflux_z_in*unit_jflux)
        Msink_list.append(Msink)
        Mdisc_list.append(np.sum(M_shell)*mass_sol)
        Sigma_list.append(M_shell/face_surf)
        rad_surf_list.append(rad_surf)
        ax[3].loglog(r_shell[1:]*len_AU,M_shell/face_surf*mass_kg/len_cm**2*1.e4)
        ax[4].loglog(r_shell[1:]*len_AU, -flux_r*unit_flux)
        #ax[4].semilogx(r_shell[1:]*len_AU, flux_z_out*unit_flux,':')
        ax[4].loglog(r_shell[1:]*len_AU, -flux_z_in*unit_flux,'--')
        print flux_r, flux_z_out, flux_z_in 
    ax[2].plot(time_list,rr_list)
    ax[2].plot(time_list,zz_list)
    ax[1].legend()
    ax[5].plot(time_list,Msink_list)
    ax[5].plot(time_list,Mdisc_list)
    ax[5].plot(time_list,[ms+md for ms, md in zip(Msink_list,Mdisc_list)])
    save_pickle(path_out+'/disc_props_csink.pkl',[time_list,Msink_list,Mdisc_list, rr_list, zz_list, Sigma_list, Rflux_list,Vfluxout_list, Vfluxin_list,jRflux_list,jVfluxout_list,jVfluxin_list,center_list,cvel_list,axis_list,ioutput_list,thres_rho_list,center_sink_list,cvel_sink_list,center_disc_list,rad_surf_list])
    save_fig(path_out+'/disc_time_csink',tight=False)
    plt.show()
#############################################################################################
## plot the tracjectories of sink, disc, and sink+disc
def trace_disc_movement(path,path_out=None,overwrite=True):
    if path_out is None: path_out = path
    filename = path_out+'/sink_trajectory'
    #filename = 'sink_trajectory'
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    outputs = search_ro(path)
    if os.path.exists(filename+'.npz'):
        sinkload = np.load(filename+'.npz')
        mass = sinkload['mass'] #Ms
        positions = sinkload['positions'] #AU
        velocities = sinkload['velocities'] #m/s
        ages = sinkload['ages'] #year
        outputs_s = sinkload['outputs_s']
    else:
        nout = len(outputs)
        mass = []
        positions = []
        velocities = []
        ages = []
        outputs_s = []
    if len(outputs_s)==0: last = 0
    else: last = outputs_s[-1]
    if outputs[-1]>last:
        if last>0:
            mass = mass.tolist()
            positions = positions.tolist()
            velocities = velocities.tolist()
            ages = ages.tolist()
            outputs_s = outputs_s.tolist()
            outputs = [out for out in outputs if out>last]
            nout = len(outputs)
        for iout in range(nout):
            ioutput = outputs[iout]
            sinks = read_sink_cvs(ioutput,path,deli=',')
            if sinks.shape[0]>=1:
                ind_star = np.argmax(sinks[:,1])
                mass.append(sinks[ind_star,1])
                positions.append(sinks[ind_star,3:6]/boxlen-0.5)
                velocities.append(sinks[ind_star,6:9])
                ages.append(sinks[ind_star,15])
                outputs_s.append(ioutput)
        mass = np.asarray(mass)*mass_sol; positions = np.asarray(positions)*len_AU; velocities = np.asarray(velocities)*vel_ms; ages=np.asarray(ages)
        np.savez(filename,mass=mass,positions=positions,velocities=velocities,ages=ages,outputs_s=outputs_s)
    istart = -20
    iend = -1
    mass = mass[istart:iend];positions=positions[istart:iend,:];velocities=velocities[istart:iend,:];ages=ages[istart:iend]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(positions[:,0],positions[:,1], positions[:,2])
    ax.quiver(positions[:,0],positions[:,1], positions[:,2],velocities[:,0],velocities[:,1],velocities[:,2],length=0.5,normalize=True,alpha=0.4)
    ax.set_xlabel('x(AU)')
    ax.set_ylabel('y(AU)')
    ax.set_zlabel('z(AU)')
    vmean = np.mean(np.linalg.norm(velocities,axis=1))
    vstd = np.std(np.linalg.norm(velocities,axis=1))
    vdrift = np.linalg.norm(positions[-1,:]-positions[0,:])/(ages[-1]-ages[0]) * 1.5e11 / (365*86400)
    plt.title(r'$v_\mathrm{local} = %i\pm%i$ (m/s), $v_\mathrm{drift} = %i$ (m/s)'%(vmean,vstd,vdrift))
    plt.savefig(path_out+'/sink_trajectory.pdf')
#    plt.show()
#    [time,Msink,Mdisc, rdisc, zdisc, Sigma_list, Rflux_list,Vfluxout_list,Vfluxin_list,jRflux_list,jVfluxout_list,jVfluxin_list,center_list,cvel_list,axis_list,outputs,thres_rho_list,center_sink_list,cvel_sink_list,center_disc_list,rad_surf_list] = read_pickle(path_out+'/disc_props_csink.pkl',21)
#np.savez(path_out+'/disc_basic_parameters', ioutput_list=ioutput_list, time_list=time_list, Msink_list=Msink_list, Mdisc_list=Mdisc_list, thres_rho_list=thres_rho_list, center_sink_list=center_sink_list, center_disc_list=center_disc_list, cvel_sink_list=cvel_sink_list, cvel_disc_list=cvel_disc_list, ax_csink_list=ax_csink_list, ax_ctot_list=ax_ctot_list, r_csink_list=r_csink_list, r_ctot_list=r_ctot_list)
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    Msink=np.asarray(discload['Msink_list'])[istart:iend]; 
    center_sink_list=(np.asarray(discload['center_sink_list'])-0.5)[istart:iend]*len_AU; 
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])[istart:iend]*vel_ms 
    Mdisc=np.asarray(discload['Mdisc_list'])[istart:iend]; 
    center_disc_list=(np.asarray(discload['center_disc_list'])[istart:iend]-0.5)*len_AU; 
    cvel_disc_list=np.asarray(discload['cvel_disc_list'])[istart:iend]*vel_ms
    ax_ctot = np.asarray(discload['ax_ctot_list'])[istart:iend]
    center_tot = (center_sink_list*Msink[:,np.newaxis]+center_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    cvel_tot = (cvel_sink_list*Msink[:,np.newaxis]+cvel_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    print cvel_tot, np.linalg.norm(cvel_tot,axis=1)
    print np.arcsin(np.sum((center_disc_list-center_sink_list)*ax_ctot,axis=1)/np.linalg.norm(center_disc_list-center_sink_list,axis=1))/np.pi*180.
    inds_nosink = np.where(Msink!=0)
    Msink = Msink[inds_nosink]
    center_sink_list = center_sink_list[inds_nosink]
    cvel_sink_list = cvel_sink_list[inds_nosink]
#    istart = -50
    ax.plot(center_sink_list[:,0],center_sink_list[:,1], center_sink_list[:,2])
    ax.quiver(center_sink_list[:,0],center_sink_list[:,1], center_sink_list[:,2],cvel_sink_list[:,0],cvel_sink_list[:,1],cvel_sink_list[:,2],length=0.5,normalize=True,alpha=0.4)
    ax.plot(center_disc_list[:,0],center_disc_list[:,1], center_disc_list[:,2])
    ax.plot(center_tot[:,0],center_tot[:,1], center_tot[:,2])
    ax.quiver(center_disc_list[:,0],center_disc_list[:,1], center_disc_list[:,2],cvel_disc_list[:,0],cvel_disc_list[:,1],cvel_disc_list[:,2],length=0.5,normalize=True,alpha=0.4)
    #ax.quiver(center_tot[:,0],center_tot[:,1], center_tot[:,2], cvel_tot[:,0], cvel_tot[:,1], cvel_tot[:,2],color='k',length=2,normalize=True,alpha=0.8)    
    ax.quiver(center_tot[:,0],center_tot[:,1], center_tot[:,2], ax_ctot[:,0], ax_ctot[:,1], ax_ctot[:,2],color='k',length=0.5,normalize=True,alpha=0.4) 
    save_fig(filename,ps=True)
    fig, ax = plt.subplots(1,3,subplot_kw=dict(aspect='equal'),figsize=(8,3))
    extend = np.amax([np.ptp(positions[:,0]),np.ptp(positions[:,1]),np.ptp(positions[:,2])])*0.54
    xc, yc, zc = map(np.median, [positions[:,0],positions[:,1],positions[:,2]])
    cen = [xc,yc,zc]
    for axis, ind, lab in zip(ax,[[0,1],[1,2],[0,2]], [['x','y'],['y','z'],['x','z']]):
        points = np.array([positions[:,ind[0]],positions[:,ind[1]]]).T.reshape(-1,1,2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=plt.get_cmap('spring'), norm=plt.Normalize(ages[0]/1.e3, ages[-1]/1.e3))
        lc.set_array(ages/1.e3)
        axis.add_collection(lc)
        axis.axis([cen[ind[0]]-extend, cen[ind[0]]+extend, cen[ind[1]]-extend, cen[ind[1]]+extend])
        axis.set_xlabel(lab[0]+'(AU)');axis.set_ylabel(lab[1]+'(AU)')
    cbar = fig.colorbar(lc,ax=ax.ravel().tolist())
    cbar.set_label('sink age (kyr)', rotation=270)
    save_fig(filename+'_projections',ps=True)
    plt.show()
#############################################################################################
def visulize_disc_flows(path,path_out=None,overwrite=True,order='<',scl=20):
    if path_out == None: path_out = path
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    istart = 3; iend = -1
    outputs = discload['ioutput_list'][istart:iend]
    Msink=np.asarray(discload['Msink_list'])[istart:iend];
    thres_rho_list = discload['thres_rho_list'][istart:iend]
    center_sink_list=(np.asarray(discload['center_sink_list']))[istart:iend];
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])[istart:iend]
    Mdisc=np.asarray(discload['Mdisc_list'])[istart:iend];
    center_disc_list=(np.asarray(discload['center_disc_list'])[istart:iend]);
    cvel_disc_list=np.asarray(discload['cvel_disc_list'])[istart:iend]
    ax_ctot = np.asarray(discload['ax_ctot_list'])[istart:iend]
    r_ctot = np.asarray(discload['r_ctot_list'])[istart:iend]
    center_tot = (center_sink_list*Msink[:,np.newaxis]+center_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    cvel_tot = (cvel_sink_list*Msink[:,np.newaxis]+cvel_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    print cvel_tot, np.linalg.norm(cvel_tot,axis=1)
    print np.arcsin(np.sum((center_disc_list-center_sink_list)*ax_ctot,axis=1)/np.linalg.norm(center_disc_list-center_sink_list,axis=1))/np.pi*180.
#    inds_nosink = np.where(Msink!=0)
#    Msink = Msink[inds_nosink]
#    center_sink_list = center_sink_list[inds_nosink]
#    cvel_sink_list = cvel_sink_list[inds_nosink]
    vec_vtot = cvel_tot/np.linalg.norm(cvel_tot,axis=1)[:,np.newaxis]
#    yvisu = np.cross(ax_ctot, vec_vtot)
    yvisu = np.cross(ax_ctot, (center_sink_list-center_disc_list))
    yvisu = yvisu/np.linalg.norm(yvisu,axis=1)[:,np.newaxis]
    xvisu = np.cross(yvisu, ax_ctot)
    xvisu = xvisu/np.linalg.norm(xvisu,axis=1)[:,np.newaxis] ## v_center projected onto disc plane as x-axis of visulation
    nout = len(outputs)
    nr = 10
    nphi = 33
    pp = np.linspace(-np.pi,np.pi, nphi)
    for iout in range(nout)[:]:
        ioutput=outputs[iout]
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        amr = ro.amr_source(["rho","vel","P","phi","g","Br","Bl"])
        center = center_sink_list[iout]
        cvel = cvel_sink_list[iout]
        #center = center_tot[iout]
        #cvel = cvel_tot[iout]
        if center[0]==0:
            center = center_disc_list[iout]
            cvel = cvel_disc_list[iout]
        vec_z = ax_ctot[iout]
        rdisc = r_ctot[iout]
        rdisc = 20/len_AU
        rr = np.linspace(0, rdisc*2,2*nr+1)
        thres_rho = thres_rho_list[iout]*5
        dd = [0, thres_rho, 1.e50]
        xp = xvisu[iout,:]; yp = yvisu[iout,:]
        disc_region = Cylinder(center,vec_z,rdisc*2,rdisc*0.6)
        amr_disc = RegionFilter(disc_region,amr)
        cell_disc = CellsToPoints(amr_disc, smallest_cell_level=scl).flatten()
        #filt_dense = lambda dset: (dset["rho"] >= thres_rho)
        #cell_dense = PointFunctionFilter(filt_dense,cell_disc)    
        #filt_diffu = lambda dset: (dset["rho"] < thres_rho)
        #cell_diffu = PointFunctionFilter(filt_diffu,cell_disc)
        #inds_disc = np.where(cell_disc["rho"]>=thres_rho)
        vel_disc = cell_disc["vel"]-cvel[np.newaxis,:]
        dV_disc = cell_disc.get_sizes()**3
        dm_disc = cell_disc["rho"]*dV_disc

        posi_disc = cell_disc.points-center[np.newaxis,:]
        dist_disc = np.linalg.norm(posi_disc,axis=1) ##3D distance
        vec_r = posi_disc / dist_disc[:,np.newaxis]
        vec_phi = np.cross(vec_z,vec_r)
        vec_phi = vec_phi / np.linalg.norm(vec_phi,axis=1)[:,np.newaxis]
        vec_rc = np.cross(vec_phi,vec_z)
        vec_rc = vec_rc / np.linalg.norm(vec_rc,axis=1)[:,np.newaxis]
        up_low = np.sign(np.sum(vec_r*vec_z[np.newaxis,:],axis=1))
        distc_disc = np.sum(posi_disc*vec_rc,axis=1) ##horizontal distance
        z_disc = np.sum(posi_disc*vec_z[np.newaxis,:],axis=1)
        vel_r = np.sum(vel_disc * vec_rc, axis=1)
        vel_phi = np.sum(vel_disc * vec_phi, axis=1)
        vel_z = np.sum(vel_disc * vec_z, axis=1)
        vel_x = np.sum(vel_disc * xp[np.newaxis,:],axis=1)
        vel_y = np.sum(vel_disc * yp[np.newaxis,:],axis=1)
        vel_zin = np.ma.masked_where(vel_z*up_low>0, vel_z*up_low)
        B_disc = (cell_disc["Br"]+cell_disc["Bl"])*0.5
        B_r = np.sum(B_disc * vec_rc, axis=1)
        B_phi = np.sum(B_disc * vec_phi, axis=1)
        B_z = np.sum(B_disc * vec_z, axis=1)
        g_r = np.sum(cell_disc["g"] * vec_rc, axis=1)
        g_phi = np.sum(cell_disc["g"] * vec_phi, axis=1)

        x_disc = np.sum(posi_disc*xp[np.newaxis,:],axis=1)
        y_disc = np.sum(posi_disc*yp[np.newaxis,:],axis=1)
        phi = np.arctan2(y_disc,x_disc)
        grid_x = lambda x: np.histogramdd((distc_disc, phi,cell_disc["rho"]), bins=[rr,pp,dd], weights=x)[0]
        grid_m = grid_x(dm_disc)
        grid_V = grid_x(dV_disc)
        grid_mrho = grid_x(cell_disc["rho"]*dm_disc)
        grid_rho_v = grid_m/grid_V
        grid_rho_m = grid_mrho/grid_m
        grid_rho = grid_rho_m
        grid_rho = np.ma.masked_where(np.isnan(grid_rho), grid_rho)
#        sig_rho = np.sqrt(grid_mrho/grid_V - grid_rho**2)
        grid_mz2 = grid_x(dm_disc*z_disc**2)
        grid_mvr = grid_x(dm_disc*vel_r)
        grid_mvphi = grid_x(dm_disc*vel_phi)
        grid_mvz = grid_x(dm_disc*vel_z*up_low)
        grid_mvzin = grid_x(dm_disc*vel_zin)
        grid_mvx = grid_x(dm_disc*vel_x)
        grid_mvy = grid_x(dm_disc*vel_y)
        grid_mvphivzin = grid_x(dm_disc*vel_phi*vel_zin)
        grid_mvr2 = grid_x(dm_disc*vel_r**2)
        grid_mvphi2 = grid_x(dm_disc*vel_phi**2)
        grid_mvz2 = grid_x(dm_disc*vel_z**2)
        grid_mvrvphi = grid_x(dm_disc*vel_r*vel_phi)
        #sig_vr = np.sqrt(grid_mvr2/grid_m - (grid_mvr/grid_m)**2)
        #sig_vphi = np.sqrt(grid_mvphi2/grid_m - (grid_mvphi/grid_m)**2)
        #sig_vz = np.sqrt(grid_mvz2/grid_m - (grid_mvz/grid_m)**2)
        grid_VBr = grid_x(dV_disc*B_r)
        grid_VBr2 = grid_x(dV_disc*B_r**2)
        grid_VBphi = grid_x(dV_disc*B_phi)
        grid_VBphi2 = grid_x(dV_disc*B_phi**2)
        grid_VBz = grid_x(dV_disc*B_z)
        grid_VBz2 = grid_x(dV_disc*B_z**2)
        grid_VBrBphi = grid_x(dV_disc*B_r*B_phi)
        grid_VBzBphi = grid_x(dV_disc*B_z*B_phi)
        grid_Vgr = grid_x(dV_disc*g_r)
        grid_Vgphi = grid_x(dV_disc*g_phi)
        grid_Vgrgphi = grid_x(dV_disc*g_r*g_phi)
        grid_VP = grid_x(dV_disc*cell_disc["P"])
        grid_VPmag = grid_x(dV_disc * np.sum(B_disc**2,axis=1))

        grid_vr = grid_mvr/grid_m
        grid_vr = np.ma.masked_where(np.isnan(grid_vr), grid_vr)
        grid_vphi = grid_mvphi/grid_m
        grid_vphi = np.ma.masked_where(np.isnan(grid_vphi), grid_vphi)
        grid_vx = grid_mvx/grid_m
        grid_vx = np.ma.masked_where(np.isnan(grid_vx), grid_vx)
        grid_vy = grid_mvy/grid_m
        grid_vy = np.ma.masked_where(np.isnan(grid_vy), grid_vy)
        grid_vz = grid_mvz/grid_m
        grid_vz = np.ma.masked_where(np.isnan(grid_vz), grid_vz)
        grid_vzin = grid_mvzin/grid_m
        grid_vzin = np.ma.masked_where(np.isnan(grid_vzin), grid_vzin)
        #sig_Br = np.sqrt(grid_mBr2/grid_m - (grid_mBr/grid_m)**2)
        #sig_Bphi = np.sqrt(grid_mBphi2/grid_m - (grid_mBphi/grid_m)**2)
        np.savez(path_out+'/grids_'+str(ioutput).zfill(5),m=grid_m, V=grid_V, mrho=grid_mrho, rho_v=grid_rho_v, rho_m=grid_rho_m, mvr=grid_mvr, mvphi=grid_mvphi, mvrvphi=grid_mvrvphi, mvz=grid_mvz, mvzin=grid_mvzin, mvx=grid_mvx, mvy=grid_mvy, mvphivzin=grid_mvphivzin, VBr=grid_VBr, VBr2=grid_VBr2, VBphi=grid_VBphi, VBphi2=grid_VBphi2, VBz=grid_VBz, VBz2=grid_VBz2, VBrBphi=grid_VBrBphi, VBzBphi=grid_VBzBphi, Vgr=grid_Vgr, Vgphi=grid_Vgphi, Vgrgphi=grid_Vgrgphi,VP=grid_VP, VPmag=grid_VPmag,threshold=thres_rho,mz2=grid_mz2,pp=pp,rr=rr)

        fig, axes = plt.subplots(3,3, subplot_kw=dict(polar=True),figsize=(8,6))
        P_grid, R_grid = np.meshgrid(pp,rr)
        pp_c = 0.5*(pp[:-1]+pp[1:]); rr_c = 0.5*(rr[:-1]+rr[1:])
        #Pc_grid, Rc_grid = np.meshgrid(pp_c,rr_c)
        dphi = (pp[1]-pp[0])*0.5
        Pc_grid, Rc_grid = np.meshgrid(pp+dphi,rr_c)
        vlim = np.maximum(abs(grid_vr[:,:,1]).max(), abs(grid_vphi[:,:,1]).max())/2.
        axes[0,0].pcolormesh(P_grid,R_grid,np.log10(grid_rho[:,:,1]))
        axes[0,0].contour(Pc_grid+dphi,Rc_grid,np.log10(np.hstack((grid_rho[:,:,1],grid_rho[:,0,np.newaxis,1]))),colors='k')
        axes[0,0].set_title(r'$\rho$')
        axes[0,1].pcolormesh(P_grid,R_grid,grid_vr[:,:,1],cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[0,1].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vr[:,:,1],grid_vr[:,0,np.newaxis,1])),colors='k')
        axes[0,1].set_title(r'$v_r$')
        axes[0,2].pcolormesh(P_grid,R_grid,grid_vphi[:,:,1],cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[0,2].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vphi[:,:,1],grid_vphi[:,0,np.newaxis,1])),colors='k')
        axes[0,2].set_title(r'$v_\phi$')
        axes[1,0].pcolormesh(P_grid,R_grid,grid_vz[:,:,0],cmap='seismic_r',vmin=-vlim*0.1,vmax=vlim*0.1)
        axes[1,0].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vz[:,:,0],grid_vz[:,0,np.newaxis,0])),colors='k')
        axes[1,0].set_title(r'$v_z$')
        axes[1,1].pcolormesh(P_grid,R_grid,grid_vx[:,:,0],cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[1,1].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vx[:,:,0],grid_vx[:,0,np.newaxis,0])),colors='k')
        axes[1,1].set_title(r'$v_x$')
        axes[1,2].pcolormesh(P_grid,R_grid,grid_vy[:,:,0],cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[1,2].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vy[:,:,0],grid_vy[:,0,np.newaxis,0])),colors='k')
        axes[1,2].set_title(r'$v_y$')
        axes[2,0].pcolormesh(P_grid,R_grid,grid_vzin[:,:,0],cmap='seismic_r',vmin=-vlim*0.1,vmax=vlim*0.1)
        axes[2,0].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vzin[:,:,0],grid_vzin[:,0,np.newaxis,0])),colors='k')
        axes[2,0].set_title(r'$v_{z,in}$')
        axes[2,1].pcolormesh(P_grid,R_grid,grid_vr[:,:,0],cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[2,1].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vr[:,:,0],grid_vr[:,0,np.newaxis,0])),colors='k')
        axes[2,1].set_title(r'$v_{r,env}$')
        axes[2,2].pcolormesh(P_grid,R_grid,grid_vphi[:,:,0],cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[2,2].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vphi[:,:,0],grid_vphi[:,0,np.newaxis,0])),colors='k')
        axes[2,2].set_title(r'$v_{\phi,env}$')
        for ax in axes.flatten(): ax.set_yticklabels([]); ax.set_xticklabels([])
#        save_fig(path_out+'/disc_visu_faceon_vdisc_'+str(ioutput).zfill(5),tight=False,ps=True)
        save_fig(path_out+'/disc_visu_faceon_'+str(ioutput).zfill(5),tight=False,ps=True)
#    plt.show()
#############################################################################################
def measure_disc_flux(path,path_out=None,overwrite=True,order='<',scl=20):
    if path_out == None: path_out = path
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    unit_flux = dens_gcc * vel_ms * 1.e3
    unit_jflux = unit_flux * vel_ms * len_cm * 1.e-2
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
#np.savez(path_out+'/disc_basic_parameters', ioutput_list=ioutput_list, time_list=time_list, Msink_list=Msink_list, Mdisc_list=Mdisc_list, thres_rho_list=thres_rho_list, center_sink_list=center_sink_list, center_disc_list=center_disc_list, cvel_sink_list=cvel_sink_list, cvel_disc_list=cvel_disc_list, ax_csink_list=ax_csink_list, ax_ctot_list=ax_ctot_list, r_csink_list=r_csink_list, r_ctot_list=r_ctot_list)
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    istart = 3; iend = -1
    outputs = discload['ioutput_list'][istart:iend]
    Msink=np.asarray(discload['Msink_list'])[istart:iend];
    thres_rho_list = discload['thres_rho_list'][istart:iend]
    center_sink_list=(np.asarray(discload['center_sink_list']))[istart:iend];
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])[istart:iend]
    Mdisc=np.asarray(discload['Mdisc_list'])[istart:iend];
    center_disc_list=(np.asarray(discload['center_disc_list'])[istart:iend]);
    cvel_disc_list=np.asarray(discload['cvel_disc_list'])[istart:iend]
    ax_csink = np.asarray(discload['ax_csink_list'])[istart:iend]
    r_csink = np.asarray(discload['r_csink_list'])[istart:iend]
    nout = len(outputs)
    nr = 20
    nphi = 33
    pp = np.linspace(-np.pi,np.pi, nphi)
    #rr_list = []
    #zz_list = []
    Vfluxout_list = []
    Vfluxin_list = []
    Rflux_list = []
    jVfluxout_list = []
    jVfluxin_list = []
    jRflux_list = []
    Sigma_list = []
    rad_surf_list = []
    if path[-8:]=='DC_1_res': isam=7
    elif  path[-8:]=='DC_2_res':isam =3
    elif  path[-8:]=='DC_6_res':isam =3 
    #fig, ax = plt.subplot(3,1)
    for iout in range(nout)[isam:isam+1]:
        fig, ax = plt.subplots(5,1,figsize=(4,6))
        ioutput=outputs[iout]
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        amr = ro.amr_source(["rho","vel","P","phi","g","Br","Bl"])
        c_center = center_sink_list[iout]
        c_vel = cvel_sink_list[iout]
        c_axis = ax_csink[iout]
        c_radius = r_csink[iout]
        c_radius = 12./len_AU
        Mstar = Msink[iout]
        aspects = [0.06,0.18,0.3,0.42,0.54]#[0.01,0.033,0.067,0.08,0.1,0.2,0.3,0.4,0.5]
        FLUXS_ZIN=[]
        M_CYN=[]
        for aspect in aspects:
            filename = path_out+'/flux_res_'+str(iout)+'_'+str(aspect)
            if os.path.exists(filename+'.npz'):
                 fileload = np.load(filename+'.npz')
                 flux_r = fileload['flux_r']
                 flux_z_in = fileload['flux_z_in']
                 flux_z_out = fileload['flux_z_out']
                 M_shell = fileload['M_shell']
                 face_surf = fileload['face_surf']
                 rad_surf = fileload['rad_surf']
                 r_shell = fileload['r_shell']
            else:
                flux_r, flux_z_out, flux_z_in, jflux_r, jflux_z_out, jflux_z_in, r_shell, M_shell, face_surf, rad_surf = cylinder_flux(ro,c_center=c_center,c_axis=c_axis,c_vel=c_vel,c_radius=c_radius*2,c_height=c_radius*aspect,n_shell=nr,n_h=8,random_sample=False,scl=scl,add_j=True)
            #ax[0].semilogy(r_shell[1:]*len_AU,flux_z_out*unit_flux,':')
            ax[0].semilogy(r_shell[1:]*len_AU, -flux_z_in*unit_flux,'--')
            ax[2].semilogy(r_shell[1:]*len_AU,-np.cumsum((flux_z_in-flux_z_out)*face_surf,axis=0)*mass_sol/time_Myr/1.e6/boxlen)
            ax[3].semilogy(r_shell[1:]*len_AU, -flux_r*unit_flux)
            #ax[2].semilogy(r_shell[1:]*len_AU, np.cumsum((-flux_z_out-flux_z_in)*face_surf)*mass_sol/time_Myr/1.e6/boxlen)
            ax[4].semilogy(r_shell[1:]*len_AU, -flux_r*rad_surf*mass_sol/time_Myr/1.e6/boxlen)
            #ax[4].semilogy(r_shell[1:]*len_AU,(np.cumsum((-flux_z_out-flux_z_in)*face_surf)-flux_r*rad_surf)*mass_sol/time_Myr/1.e6/boxlen)
            #ax[5].plot(r_shell[1:]*len_AU,Mstar+np.cumsum(M_shell)*mass_sol)
            np.savez(filename,flux_z_out=flux_z_out,flux_z_in=flux_z_in,flux_r=flux_r,M_shell=M_shell,face_surf=face_surf,rad_surf=rad_surf,r_shell=r_shell)
            FLUXS_ZIN.append(flux_z_in-flux_z_out)
            M_CYN.append(M_shell)
        FLUXS_ZIN = np.asarray(FLUXS_ZIN)
        M_CYN = np.asarray(M_CYN)
        print np.asarray(aspects).shape, r_shell[1:].shape, FLUXS_ZIN.shape
        func_flux = interp2d(np.asarray(aspects)*c_radius,r_shell[1:],FLUXS_ZIN.T)
        gridfile = path_out+'/grids_'+str(ioutput).zfill(5)+'.npz'
        grids = np.load(gridfile)
        rr = grids['rr']
        rc = 0.5*(rr[:-1]+rr[1:])
        M_shell = np.sum(grids['m'],axis=1)
        h_r = np.sqrt(np.sum(grids['mz2'][:,:,1],axis=1)/M_shell[:,1])
        #h_r = h_r*5
        FLUX_surf = func_flux(h_r,rc).diagonal()
#        FLUX_surf = func_flux(np.asarray(aspects)*c_radius,r_shell[1:])
        print np.asarray(aspects)*c_radius*len_AU
        print r_shell[1:]*len_AU
        print h_r.shape, rc.shape, FLUX_surf.shape
        print h_r*len_AU
        print rc*len_AU
        ax[1].semilogy(np.asarray(aspects)*c_radius*len_AU,-FLUXS_ZIN[:,:]*unit_flux)
        ax[0].semilogy(rc*len_AU,-FLUX_surf*unit_flux,'-o',c='k',lw=1)
        ax[1].semilogy(h_r*len_AU,-FLUX_surf*unit_flux,'-o',c='k',lw=1)
        print FLUXS_ZIN.shape, face_surf.shape, (FLUXS_ZIN*face_surf[np.newaxis,:]).shape, np.cumsum(FLUXS_ZIN*face_surf[np.newaxis,:],axis=0).shape
        #ax[2].semilogy(rc*len_AU,-np.cumsum(FLUXS_ZIN*face_surf[np.newaxis,:],axis=0)*mass_sol/time_Myr/1.e6/boxlen)
        h_r = h_r*5
        FLUX_surf = func_flux(h_r,rc).diagonal()
        ax[0].semilogy(rc*len_AU,-FLUX_surf*unit_flux,c='k',lw=1)
        ax[1].semilogy(h_r*len_AU,-FLUX_surf*unit_flux,c='k',lw=1)
#        ax[0].semilogy(r_shell[1:]*len_AU,-FLUX_surf*unit_flux,c='k',lw=1)
#        ax[1].semilogy(np.asarray(aspects)*c_radius*len_AU,-FLUX_surf.T*unit_flux,c='k',lw=1)
        plt.show()
#############################################################################################
def write_grids(path,path_out=None,overwrite=True,order='<',scl=20):
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    unit_flux = dens_gcc * vel_ms * 1.e3
    unit_jflux = unit_flux * vel_ms * len_cm * 1.e-2
    unit_j = vel_ms * len_cm * 1.e-2
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    istart = 0; iend = -1
    outputs = discload['ioutput_list'][istart:iend]
    Msink=np.asarray(discload['Msink_list'])[istart:iend];
    thres_rho_list = discload['thres_rho_list'][istart:iend]
    center_sink_list=(np.asarray(discload['center_sink_list']))[istart:iend];
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])[istart:iend]
    Mdisc=np.asarray(discload['Mdisc_list'])[istart:iend];
    ax_csink = np.asarray(discload['ax_csink_list'])[istart:iend]
    nout = len(outputs)
    time_list = discload['time_list'][istart:iend]
    dx = 1./2**lmax
    rdisc = 40./len_AU
    nr = int(rdisc/dx/4)
    rr = np.linspace(0, rdisc, 2*nr+1)
    redges = rr[::2]
    rcenters = rr[1::2]
    zdisc = 20/len_AU
    nz = int(zdisc/dx/4)
    zz = np.linspace(0,zdisc,2*nz+2)
    zedges = np.append(0,zz[1::2])
    zcenters = zz[::2]
    print 'nr = ', nr, 'nz = ', nz
    #overwrite=False
    #overwrite=True
    isam = 1
    print Mdisc
    if path[-8:]=='DC_1_res': isam=7
    elif  path[-8:]=='DC_2_res':isam =7
    elif  path[-8:]=='DC_6_res':isam =7
    for iout in range(nout)[isam:]:
#        np.savez(path_out+'/grids_'+str(ioutput).zfill(5),m=grid_m, V=grid_V, mrho=grid_mrho, rho_v=grid_rho_v, rho_m=grid_rho_m, mvr=grid_mvr, mvphi=grid_mvphi, mvz=grid_mvz, mvzin=grid_mvzin, mvx=grid_mvx, mvy=grid_mvy, mvphivzin=grid_mvphivzin, VBr=grid_VBr, VBr2=grid_VBr2, VBphi=grid_VBphi, VBphi2=grid_VBphi2, VBz=grid_VBz, VBz2=grid_VBz2, VBrBphi=grid_VBrBphi, Vgr=grid_Vgr, Vgphi=grid_Vgphi, Vgrgphi=grid_Vgrgphi,VP=grid_VP, VPmag=grid_VPmag,threshold=thres_rho,mz2=grid_mz2,pp=pp,rr=rr)
        ioutput = outputs[iout]
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        amr = ro.amr_source(["rho","vel","P","phi","g","Br","Bl"])
        center = center_sink_list[iout]
        cvel = cvel_sink_list[iout]
        vec_z = ax_csink[iout]
        filedir = path_out+'/grid_rz_'+str(ioutput).zfill(5)+'/'
        fileh = 'grid_rz_'+str(ioutput).zfill(5)
        filef = '.dat'
        if overwrite:
            if not os.path.exists(filedir): os.mkdir(filedir)
            np.savetxt(filedir+'/timeyear_MsinkMs.dat',np.array([time_list[iout]*1.e6, Msink[iout]]))
            disc_region = Cylinder(center,vec_z,rdisc,zdisc*2)
            amr_disc = RegionFilter(disc_region,amr)
            cell_disc = CellsToPoints(amr_disc, smallest_cell_level=scl).flatten()
            vel_disc = cell_disc["vel"]-cvel[np.newaxis,:]
            dV_disc = cell_disc.get_sizes()**3
            dm_disc = cell_disc["rho"]*dV_disc
            posi_disc = cell_disc.points-center[np.newaxis,:]
            dist_disc = np.linalg.norm(posi_disc,axis=1) ##3D distance
            vec_r = posi_disc / dist_disc[:,np.newaxis]
            vec_phi = np.cross(vec_z,vec_r)
            vec_phi = vec_phi / np.linalg.norm(vec_phi,axis=1)[:,np.newaxis]
            vec_rc = np.cross(vec_phi,vec_z)
            vec_rc = vec_rc / np.linalg.norm(vec_rc,axis=1)[:,np.newaxis]
            up_low = np.sign(np.sum(vec_r*vec_z[np.newaxis,:],axis=1))
            distc_disc = np.sum(posi_disc*vec_rc,axis=1) ##horizontal distance
            z_disc = np.sum(posi_disc*vec_z[np.newaxis,:],axis=1) ##vertical distance
            zabs_disc = abs(z_disc)
            vel_r = np.sum(vel_disc * vec_rc, axis=1)
            vel_phi = np.sum(vel_disc * vec_phi, axis=1)
            vel_z = np.sum(vel_disc * vec_z, axis=1)
            vel_zsym = vel_z * up_low
            cs2 = cell_disc["P"]/cell_disc["rho"]
            jz = distc_disc * vel_phi
            #B_disc = (cell_disc["Br"]+cell_disc["Bl"])*0.5
            #B_r = np.sum(B_disc * vec_rc, axis=1)
            #B_phi = np.sum(B_disc * vec_phi, axis=1)
            #B_z = np.sum(B_disc * vec_z, axis=1)

            grid_x = lambda x: np.histogramdd((distc_disc, zabs_disc), bins=[redges,zedges], weights=x)[0]
            grid_m = grid_x(dm_disc)
            grid_V = grid_x(dV_disc)
            print grid_V
            grid_rho = grid_m/grid_V
            #grid_rho = np.ma.masked_where(np.isnan(grid_rho), grid_rho)
            grid_mvr = grid_x(dm_disc*vel_r)
            grid_mvphi = grid_x(dm_disc*vel_phi)
            grid_mvz = grid_x(dm_disc*vel_zsym)
            grid_PV = grid_x(cell_disc["P"]*dV_disc)
            grid_cs = np.sqrt(grid_PV/grid_m)
            grid_Jz = grid_x(dm_disc*jz)
            grid_jz = grid_Jz/grid_m
            grid_Jzvz = grid_x(dm_disc*jz*vel_zsym)
            grid_jzvz = grid_Jzvz/grid_m

            #np.savetxt(filedir+'/timeyear_Mskg.dat',np.array(time_list[iout]*time_Myr*1.e3, Msink[iout]*mass_kg))
            np.savetxt(filedir+'/rs.dat', rcenters*len_AU)
            np.savetxt(filedir+'/zs.dat', zcenters*len_AU)
            np.savetxt(filedir+'/rho.dat',grid_rho*dens_gcc*1.e3)
            np.savetxt(filedir+'/vr.dat',grid_mvr/grid_m)
            np.savetxt(filedir+'/vphi.dat',grid_mvphi/grid_m)
            np.savetxt(filedir+'/vz.dat',grid_mvz/grid_m*vel_ms)
            np.savetxt(filedir+'/cs.dat',grid_cs*vel_ms)
            np.savetxt(filedir+'/jz.dat',grid_jz*unit_j)
            np.savetxt(filedir+'/jzvz.dat',grid_jzvz*unit_j*vel_ms)

            grid_x = lambda x: np.histogramdd((distc_disc, zabs_disc), bins=[redges,zcenters], weights=x)[0]
            grid_m = grid_x(dm_disc)
            grid_m = np.cumsum(grid_m,axis=1)
            grid_m = np.pad(grid_m, ((0, 0), (1, 0)), 'constant', constant_values=(0.))
            np.savetxt(filedir+'/mz.dat',grid_m*mass_kg)

            
#############################################################################################
# calculate flux at all position and evaluate at disc surface
def gaussian(x, a, b, c):
    val = a * np.exp(-(x - b)**2 / c**2)
    return val
def erfunction(n0,z,sig):
    val = n0 * sp.special.erf(z/sig/np.sqrt(2))
def measure_flux(path,path_out=None,overwrite=True,order='<',scl=20):
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    unit_flux = dens_gcc * vel_ms * 1.e3
    unit_jflux = unit_flux * vel_ms * len_cm * 1.e-2
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    istart = 3; iend = -1
    outputs = discload['ioutput_list'][istart:iend]
    Msink=np.asarray(discload['Msink_list'])[istart:iend];
    thres_rho_list = discload['thres_rho_list'][istart:iend]
    center_sink_list=(np.asarray(discload['center_sink_list']))[istart:iend];
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])[istart:iend]
    Mdisc=np.asarray(discload['Mdisc_list'])[istart:iend];
    center_disc_list=(np.asarray(discload['center_disc_list'])[istart:iend]);
    cvel_disc_list=np.asarray(discload['cvel_disc_list'])[istart:iend]
    ax_csink = np.asarray(discload['ax_csink_list'])[istart:iend]
    r_ctot = np.asarray(discload['r_ctot_list'])[istart:iend]
    center_tot = (center_sink_list*Msink[:,np.newaxis]+center_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    cvel_tot = (cvel_sink_list*Msink[:,np.newaxis]+cvel_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    nout = len(outputs)
    time_list = discload['time_list'][istart:iend]
    Sigma_t = []
    Mdot_r_t = []
    Mdot_v_t = []
    Source_t = []
    Sourcein_t = []
    Omega_t = []
    Omega_k_t = []
    Qtoom_t = []
    beta_t = []
    h_r_t = []
    alpha_Rey_t = []
    alpha_Max_t = []
    alpha_Maxz_t = []
    alpha_Poi_t = []
    time_t = []
    SigmaR2W = []; R2Omega = []
    nr = 6
    nz = 12
    isam=7
    if path[-8:]=='DC_1_res': isam=7
    elif  path[-8:]=='DC_2_res':isam =7
    elif  path[-8:]=='DC_6_res':isam =7
    for iout in range(nout)[isam:isam+1]:
#        np.savez(path_out+'/grids_'+str(ioutput).zfill(5),m=grid_m, V=grid_V, mrho=grid_mrho, rho_v=grid_rho_v, rho_m=grid_rho_m, mvr=grid_mvr, mvphi=grid_mvphi, mvz=grid_mvz, mvzin=grid_mvzin, mvx=grid_mvx, mvy=grid_mvy, mvphivzin=grid_mvphivzin, VBr=grid_VBr, VBr2=grid_VBr2, VBphi=grid_VBphi, VBphi2=grid_VBphi2, VBz=grid_VBz, VBz2=grid_VBz2, VBrBphi=grid_VBrBphi, Vgr=grid_Vgr, Vgphi=grid_Vgphi, Vgrgphi=grid_Vgrgphi,VP=grid_VP, VPmag=grid_VPmag,threshold=thres_rho,mz2=grid_mz2,pp=pp,rr=rr)
        ioutput = outputs[iout]
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        dx = 1./2**ro.info["levelmax"]
        amr = ro.amr_source(["rho","vel","P","phi","g","Br","Bl"])
        center = center_sink_list[iout]
        cvel = cvel_sink_list[iout]
        vec_z = ax_csink[iout]
        rdisc = 8./len_AU
        rr = np.linspace(0, rdisc*2,2*nr+1)
        zdisc = rdisc/2.
        zz = np.linspace(-zdisc,zdisc,2*nz+1)
        thres_rho = thres_rho_list[iout]*8
        dd = [0, thres_rho, 1.e50]
        overwrite=False
        #overwrite=True
        gridfile = path_out+'/grids_rz_'+str(ioutput).zfill(5)+'.npz'
        if overwrite and os.path.exists(gridfile): os.remove(gridfile)
        if (not os.path.exists(gridfile)):
            disc_region = Cylinder(center,vec_z,rdisc*2,zdisc*2)
            amr_disc = RegionFilter(disc_region,amr)
            cell_disc = CellsToPoints(amr_disc, smallest_cell_level=scl).flatten()
            vel_disc = cell_disc["vel"]-cvel[np.newaxis,:]
            dV_disc = cell_disc.get_sizes()**3
            dm_disc = cell_disc["rho"]*dV_disc
    
            posi_disc = cell_disc.points-center[np.newaxis,:]
            dist_disc = np.linalg.norm(posi_disc,axis=1) ##3D distance
            vec_r = posi_disc / dist_disc[:,np.newaxis]
            vec_phi = np.cross(vec_z,vec_r)
            vec_phi = vec_phi / np.linalg.norm(vec_phi,axis=1)[:,np.newaxis]
            vec_rc = np.cross(vec_phi,vec_z)
            vec_rc = vec_rc / np.linalg.norm(vec_rc,axis=1)[:,np.newaxis]
            up_low = np.sign(np.sum(vec_r*vec_z[np.newaxis,:],axis=1))
            distc_disc = np.sum(posi_disc*vec_rc,axis=1) ##horizontal distance
            z_disc = np.sum(posi_disc*vec_z[np.newaxis,:],axis=1)
            vel_r = np.sum(vel_disc * vec_rc, axis=1)
            #vel_phi = np.sum(vel_disc * vec_phi, axis=1)
            vel_z = np.sum(vel_disc * vec_z, axis=1)
            #vel_x = np.sum(vel_disc * xp[np.newaxis,:],axis=1)
            #vel_y = np.sum(vel_disc * yp[np.newaxis,:],axis=1)
            vel_zin = np.ma.masked_where(vel_z*up_low>0, vel_z*up_low)
            #B_disc = (cell_disc["Br"]+cell_disc["Bl"])*0.5
            #B_r = np.sum(B_disc * vec_rc, axis=1)
            #B_phi = np.sum(B_disc * vec_phi, axis=1)
            #B_z = np.sum(B_disc * vec_z, axis=1)
            #g_r = np.sum(cell_disc["g"] * vec_rc, axis=1)
            #g_phi = np.sum(cell_disc["g"] * vec_phi, axis=1)

            #x_disc = np.sum(posi_disc*xp[np.newaxis,:],axis=1)
            #y_disc = np.sum(posi_disc*yp[np.newaxis,:],axis=1)
            #phi = np.arctan2(y_disc,x_disc)
            grid_x = lambda x: np.histogramdd((distc_disc, z_disc,cell_disc["rho"]), bins=[rr,zz,dd], weights=x)[0]
            grid_m = grid_x(dm_disc)
            grid_mz2 = grid_x(dm_disc*z_disc**2)
            grid_V = grid_x(dV_disc)
            grid_mrho = grid_x(cell_disc["rho"]*dm_disc)
            grid_rho_v = grid_m/grid_V
            grid_rho_m = grid_mrho/grid_m
            grid_rho = grid_rho_m
            grid_rho = np.ma.masked_where(np.isnan(grid_rho), grid_rho)
    #        sig_rho = np.sqrt(grid_mrho/grid_V - grid_rho**2)
            #grid_mz2 = grid_x(dm_disc*z_disc**2)
            grid_mvr = grid_x(dm_disc*vel_r)
            #grid_mvphi = grid_x(dm_disc*vel_phi)
            grid_mvz = grid_x(dm_disc*vel_z*up_low)
            grid_mvzin = grid_x(dm_disc*vel_zin)
            grid_mz = grid_x(dm_disc*abs(z_disc))
            grid_z10 = grid_x(z_disc**40)
            np.savez(gridfile,rr=rr,zz=zz,grid_m=grid_m,grid_V=grid_V,grid_mrho=grid_mrho,grid_rho_v=grid_rho_v,grid_rho_m=grid_rho_m,grid_mvr=grid_mvr,grid_mvz=grid_mvz,grid_mvzin=grid_mvzin,grid_mz2=grid_mz2,grid_mz=grid_mz,grid_z10=grid_z10)
        else:
            grids = np.load(gridfile)
            rr = grids['rr']; zz = grids['zz']
            grid_m=grids['grid_m']
            grid_V=grids['grid_V']
            grid_mrho=grids['grid_mrho']
            grid_rho_v=grids['grid_rho_v']
            grid_rho_m=grids['grid_rho_m']
            grid_mvr=grids['grid_mvr']
            grid_mvz=grids['grid_mvz']
            grid_mvzin=grids['grid_mvzin'] 
            grid_mz2=grids['grid_mz2']
            grid_mz=grids['grid_mz']
            grid_z10=grids['grid_z10']
        dr = np.diff(rr); rc = 0.5*(rr[:-1]+rr[1:])
        dz = np.diff(zz); zc = 0.5*(zz[:-1]+zz[1:])
        Surf_z = np.pi*np.diff(rr**2)
        Surf_r = 2*np.pi*rc[:,np.newaxis]*dz[np.newaxis,:]
        M_shell = np.sum(grid_m,axis=1)
        V_shell = np.sum(grid_V,axis=1)
        grid_rhovr = grid_mvr/grid_V
        grid_rhovr = np.ma.masked_where(np.isnan(grid_rhovr), grid_rhovr)
        #grid_vphi = grids['mvphi']/grids['m']
        #grid_vphi = np.ma.masked_where(np.isnan(grid_vphi), grid_vphi)
        #grid_vx = grid_mvx/grid_m
        #grid_vx = np.ma.masked_where(np.isnan(grid_vx), grid_vx)
        #grid_vy = grid_mvy/grid_m
        #grid_vy = np.ma.masked_where(np.isnan(grid_vy), grid_vy)
        grid_rhovz = grid_mvz/grid_V
        #grid_rhovz = np.ma.masked_where(np.isnan(grid_rhovz), grid_rhovz)
        grid_rhovzin = grid_mvzin/grid_V
        #grid_rhovzin = np.ma.masked_where(np.isnan(grid_rhovzin), grid_rhovzin)
        mirrow = lambda x: x[:,:nz,:] + x[:,-1:-nz-1:-1,:]
        grid_rhovz2 = mirrow(grid_mvz)/mirrow(grid_V)
        #grid_rhovz2 = np.ma.masked_where(np.isnan(grid_rhovz2), grid_rhovz2)
        grid_rhovzin2 = mirrow(grid_mvzin)/mirrow(grid_V)
        #grid_rhovzin2 = np.ma.masked_where(np.isnan(grid_rhovzin2), grid_rhovzin2)
        grid_rhovr2 = mirrow(grid_mvr)/mirrow(grid_V)
        #grid_rhovr2 = np.ma.masked_where(np.isnan(grid_rhovr2), grid_rhovr2)
        grid_rhovzin_all = np.sum(mirrow(grid_mvzin),axis=2)/np.sum(mirrow(grid_V),axis=2)
        #grid_rhovzin_all = np.ma.masked_where(np.isnan(grid_rhovzin_all), grid_rhovzin_all)
        grid_rhovr_all = np.sum(mirrow(grid_mvr),axis=2)/np.sum(mirrow(grid_V),axis=2)
        #grid_rhovr_all = np.ma.masked_where(np.isnan(grid_rhovr_all), grid_rhovr_all)
        grid_mz22 = mirrow(grid_mz2)
        grid_mzm = mirrow(grid_mz)
        grid_m2 = mirrow(grid_m)
        sum_mz2 = np.cumsum(grid_mz22[:,::-1,:],axis=1)
        sum_mzm = np.cumsum(grid_mzm[:,::-1,:],axis=1)
        sum_m2 = np.cumsum(grid_m2[:,::-1,:],axis=1)
        #h_ex = np.sqrt((sum_mz2[:,:,1] - zc[np.newaxis,nz:]*sum_mzm[:,:,1])/sum_m2[:,:,1])*np.sqrt(np.pi/2.)
        #print h_ex*len_AU
        #h_all_ex = np.sqrt((np.sum(sum_mz2,axis=2) - zc[np.newaxis,nz:]*np.sum(sum_mzm,axis=2))/np.sum(sum_m2,axis=2))*np.sqrt(np.pi/2.)
        #print h_all_ex*len_AU
        h2_dense = np.sqrt(sum_mz2[:,:,1]/sum_m2[:,:,1])*np.sqrt(np.pi/2.)
        h_dense = (sum_mzm[:,:,1]/sum_m2[:,:,1])*np.pi/2.
        h_all = (np.sum(sum_mzm,axis=2)/np.sum(sum_m2,axis=2))*np.pi/2.
        h2_all = np.sqrt(np.sum(sum_mz2,axis=2)/np.sum(sum_m2,axis=2))*np.sqrt(np.pi/2.)
        h_r = np.amax(h_dense,axis=1)
        h_max_r = np.sum(grid_z10[:,:,1],axis=1)**(1./40)
        h_max_r[h_max_r>0.]+=dx*0.5
        print h_r*len_AU
        print h_max_r * len_AU
        print zc[nz:]* len_AU
        #func_flux = interp2d(zc[nz:],rc,grid_rhovzin_all[:,::-1],kind='linear')
        #FLUX_surf = func_flux(h_r,rc).diagonal()
        #print FLUX_surf* unit_flux
        FLUX_surf = np.zeros_like(rc)
        FLUX_surf2 = np.zeros_like(rc)
        for ir in range(nr*2):
            func_flux = interp1d(zc[nz:],grid_rhovzin_all[ir,::-1])
            FLUX_surf[ir] = func_flux(h_r[ir])
            if h_max_r[ir]>zc[nz]: FLUX_surf2[ir] = func_flux(h_max_r[ir])
        print FLUX_surf* unit_flux
        print rc*len_AU, zc[nz:]*len_AU
        stream_line = np.zeros_like(grid_rhovzin_all)
        stream_line[:,0] = rc[:]
        for iz in range(1,nz):
            rhovz = 0.5*(grid_rhovzin_all[:,-iz]+grid_rhovzin_all[:,-iz-1])
            rhovr = 0.5*(grid_rhovr_all[:,-iz]+grid_rhovr_all[:,-iz-1])
            #print 
            stream_line[:,iz] = stream_line[:,iz-1] + (zc[iz]-zc[iz-1])*rhovr/rhovz
            stream_line[np.where(stream_line[:,iz]>rc[-1]),iz] = rc[-1]
        #print grid_rhovzin_all[:,::-1]* unit_flux
        fig, ax = plt.subplots(3,1)
        ax[0].plot(rc*len_AU,np.sqrt(grid_mz22[:,::-1,1]/grid_m2[:,::-1,1])*len_AU)
        ax[0].set_prop_cycle(None)
        ax[0].plot(rc*len_AU,np.sqrt(grid_mz22[:,::-1,0]/grid_m2[:,::-1,0])*len_AU,':')
        #ax[0].set_prop_cycle(None)
        #ax[0].plot(rc*len_AU,h_dense*len_AU)
        #ax[0].set_prop_cycle(None)
        #ax[0].plot(rc*len_AU,h10_dense*len_AU)
        ax[1].plot(rc*len_AU,h2_dense*len_AU)
        ax[1].set_prop_cycle(None)
        ax[1].plot(rc*len_AU,h2_all*len_AU,':')
        ax[2].plot(rc*len_AU,h_dense*len_AU)
        ax[2].set_prop_cycle(None)
        ax[2].plot(rc*len_AU,h_all*len_AU,':')
        #ax[1].set_prop_cycle(None)
        #ax[1].plot(rc*len_AU,h_ex*len_AU,'-*')
        #ax[2].plot(zc[nz:]*len_AU,h_all.T*len_AU)
        for axi in ax[:]: axi.set_yscale('log')
        #plt.show()
        #ax[0,0].plot(rc*len_AU,-grid_rhovz[:,:,0]*unit_flux)
        #ax[0,0].set_prop_cycle(None)
        #ax[0,0].plot(rc*len_AU,grid_rhovz[:,:,0]*unit_flux,':')
        #ax[1].plot(rc*len_AU,-grid_rhovz[:,:,1]*unit_flux)
        #ax[1].set_color_cycle(None)
        #ax[1].plot(rc*len_AU,grid_rhovz[:,:,1]*unit_flux,':')
        fig, ax = plt.subplots(4,2,figsize=(5,8))
        ax[1,0].plot(rc*len_AU,-grid_rhovzin2[:,::-1,0]*unit_flux)
        ax[1,0].plot(rc*len_AU,-FLUX_surf*unit_flux,'k')
        ax[1,0].plot(rc*len_AU,-FLUX_surf2*unit_flux,'k:')
        ax[1,0].set_prop_cycle(None)
        ax[1,0].plot(rc*len_AU,-grid_rhovzin_all[:,::-1]*unit_flux,':')
        ax[2,0].plot(rc*len_AU,-grid_rhovzin2[:,::-1,1]*unit_flux)
        ax[1,1].plot(rc*len_AU,-grid_rhovr2[:,::-1,0]*unit_flux)
        ax[2,1].plot(rc*len_AU,-grid_rhovr2[:,::-1,1]*unit_flux)
        ax[0,1].plot(rc*len_AU,grid_m[:,:,1]*mass_sol)
        ax[0,1].set_prop_cycle(None)
        ax[0,1].plot(rc*len_AU,grid_m[:,:,0]*mass_sol,':')
        for axi in ax[:,0]: axi.set_yscale('log')
        ax[0,0].plot(zc[-1:-nz-1:-1]*len_AU,-grid_rhovzin2[:nr,:,0].T*unit_flux)
        ax[2,0].legend((zc[nz:]*len_AU).tolist())
        ax[3,0].plot(rc*len_AU,-np.cumsum(grid_rhovzin_all[:,::-1]*Surf_z[:,np.newaxis],axis=0)*mass_sol/time_Myr/1.e6/boxlen)
        ax[3,0].plot(rc*len_AU,-np.cumsum(FLUX_surf2*Surf_z)*mass_sol/time_Myr/1.e6/boxlen,'k:')
        ax[3,1].plot(rc*len_AU,-np.cumsum(grid_rhovr_all[:,::-1]*Surf_r[:,nz:],axis=1)*mass_sol/time_Myr/1.e6/boxlen)
        #ax[2,0].plot(rc*len_AU,-grid_rhovzin_all[:,::-1]*Surf_z[:,np.newaxis]*mass_sol/time_Myr/1.e6/boxlen)
        #ax[2,1].plot(rc*len_AU,-grid_rhovr_all[:,::-1]*Surf_r[:,nz:]*mass_sol/time_Myr/1.e6/boxlen)
        for axi in ax[:,1]: axi.set_yscale('log')
        fig2,ax2 = plt.subplots(3,1,figsize=(5,6))
        r_grid, z_grid = np.meshgrid(rc*len_AU,zc[nz:]*len_AU)
        #logs = lambda x: np.sign(x) * np.log10(abs(x))
        #grid_rhovr2=logs(grid_rhovr2);grid_rhovz2=logs(grid_rhovz2)
        h_max_r = np.ma.masked_where(h_max_r==0., h_max_r)
        ax2[0].quiver(r_grid,z_grid,grid_rhovr2[:,::-1,0].T,grid_rhovz2[:,::-1,0].T,pivot='mid')
        ax2[1].quiver(r_grid,z_grid,grid_rhovr2[:,::-1,1].T,grid_rhovz2[:,::-1,1].T,pivot='mid')
        ax2[1].plot(rc*len_AU,h_max_r*len_AU)
        ax2[1].plot(rc*len_AU,h_r*len_AU)
        ax2[2].quiver(r_grid,z_grid,grid_rhovr_all[:,::-1].T,grid_rhovzin_all[:,::-1].T,pivot='mid')
        ax2[2].plot(rc*len_AU,h_max_r*len_AU)
        ax2[2].plot(rc*len_AU,h_r*len_AU)
        print stream_line.T*len_AU
        ax2[2].plot(stream_line.T*len_AU, zc[nz:]*len_AU)
        plt.show()
#############################################################################################
#measure flux around in particle in across a cubic surface
def flux_cube(path,path_out=None,outexist=None,overwrite=True,order='<',scl=20,Msink_sam=None,time_sam=None,out_sam=None):
    fileparas = path_out+'/disc_basic_parameters.npz'
    if os.path.exists(fileparas):
        discload = np.load(fileparas)
        outputs = discload['ioutput_list'].tolist()
        time_list = discload['time_list'].tolist()
        Msink_list = discload['Msink_list'].tolist()
        center_sink_list = discload['center_sink_list'].tolist()
        cvel_sink_list = discload['cvel_sink_list'].tolist()
    nout = len(outputs)
    ir_samples = np.array([1,3,5,7,9,11,13])
    nr = len(ir_samples)
    flux_time = []
    print outputs
    nsam = max(1,nout/30)
    for isam in range(nsam,nout,nsam):
        ioutput = outputs[isam]
        if outexist is not None: 
            if ioutput not in outexist: continue
        time = time_list[isam]
        units = norm_units(path,ioutput=ioutput)
	len_AU = units["len_AU"]
        boxlen = units["boxlen"]
        Ms_Myr = units["mass_sol"] / units["time_Myr"]
        Ms_yr = Ms_Myr * 1.e-6
        kg_m2 = units["mass_kg"] / units["len_m"]**2
        vel_ms = units["vel_ms"]
        unit_dens = kg_m2 / units["len_m"]
        unit_mflux  = unit_dens * vel_ms
        dx = 1./2**units["lmax"]
        dS = dx**2
        c_center = np.asarray(center_sink_list[isam])
        cvel = np.asarray(cvel_sink_list[isam])
        print c_center, cvel
        mass_sink = Msink_list[isam] / units["mass_sol"]
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        amr = ro.amr_source(["rho","vel"])
        psamp = PointSamplingProcessor(amr)
        flux_samples = np.zeros((2,nr))
        for i in range(nr):
            ir = ir_samples[i]
            rsur = ir*dx
            ps = (np.arange(2.*ir)-ir+0.5)*dx
            ys, zs = np.meshgrid(ps,ps)
            surf = np.vstack((ir*dx*np.ones(4*ir**2),ys.flatten(),zs.flatten())).T
            print surf.shape, 'surf.shape'
            mflux = 0.
            mflux_rel = 0.
            for idir in range(3):
              surfdir = np.roll(surf,idir,axis=1)
              for rl in [1,-1]:
    #            print surfdir*rl
                rsurf = c_center[np.newaxis,:]+surfdir*rl
                sample_dset = psamp.process(rsurf,add_level=True)
                rflux = np.sum( sample_dset["rho"]*sample_dset["vel"][:,idir]*rl )
                rflux_rel = np.sum( sample_dset["rho"]*(sample_dset["vel"][:,idir]-cvel[idir]) *rl)
                print 'rflux', rflux*unit_mflux
                mflux -= rflux
                mflux_rel -= rflux_rel
            flux_samples[0,i] = mflux; flux_samples[1,i] = mflux_rel
        plt.plot(ir_samples,flux_samples[0,:]*dS*Ms_yr,'C'+str(isam%7),label='%.1f kyr'%(time* units["time_Myr"]*1.e3))
        plt.plot(ir_samples,flux_samples[1,:]*dS*Ms_yr,'C'+str(isam%7)+':')
        flux_time.append([time* units["time_Myr"],flux_samples*dS*Ms_yr])       
        #print 'time', time* units["time_Myr"], 'flux', mflux * unit_mflux, mflux_rel * unit_mflux
    for itime in flux_time:
        print itime[0]
        print itime[1]
    plt.yscale('symlog',linthreshy=1.e-9)
    plt.legend()
    plt.savefig(path_out+'/flux_cube_r.pdf')
    plt.show()


#############################################################################################
# measure and plot disc properties as function of radius
def some_plots(path,path_out=None,overwrite=True,order='<',scl=20,Msink_sam=None,time_sam=None,out_sam=None):
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    outputs = discload['ioutput_list']
    time_list = discload['time_list']
    Msink=np.asarray(discload['Msink_list'])
    thres_rho_list = discload['thres_rho_list']
    center_sink_list=(np.asarray(discload['center_sink_list']))
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])
    ax_list = np.asarray(discload['ax_csink_list'])
    r_list = np.asarray(discload['r_csink_list'])
    nout = len(outputs)
    outstr = [str(o) for o in outputs]
    print ",".join(outstr)
    for i in range(nout):
        print i, outputs[i], time_list[i]*1.e3, Msink[i]
    if time_sam: isam = np.argmin(abs(time_list-time_sam))
    elif Msink_sam: isam = np.argmin(abs(Msink-Msink_sam))
    elif out_sam: isam = np.where(outputs==out_sam)[0][0]
    else: isam = np.argmin(abs(time_list-0.10098)) #0.10098
    print isam, outputs[isam], time_list[isam]*1.e3, Msink[isam]


    ioutput = outputs[isam]
#    filename = path_out+'/disk_data_'+str(ioutput).zfill(5)+'.npz'
#    if os.path.exists(filename): return

    units = norm_units(path,ioutput=ioutput)
    len_AU = units["len_AU"]
    boxlen = units["boxlen"]
    kg_m2 = units["mass_kg"] / units["len_m"]**2
    vel_ms = units["vel_ms"]
    kg_m2_s = kg_m2/units["time_s"]
    kg_s = units["mass_kg"] / units["time_s"]
    Ms_Myr = units["mass_sol"] / units["time_Myr"]
    Ms_yr = Ms_Myr * 1.e-6
    unit_dens = kg_m2 / units["len_m"]
    unit_mflux  = unit_dens * vel_ms
    unit_J = kg_m2 * vel_ms
    unit_jflux = unit_J * vel_ms
    unit_jspe = vel_ms * units["len_m"]
    unit_B = units["mag_gauss"]
    dx = 1./2**units["lmax"]
    xmin = 0.1; xmax = 100.
    if ioutput > 1800: xmax = 400.
    if ioutput > 2500: xmax = 1000.
    if ioutput > 2900: xmax = 1000.
    npad = 64
    r_sel = xmax / len_AU 
    r_cyl = r_sel + npad*dx
    z_sel = 100 / len_AU
    print 'rmax', r_cyl *len_AU
    #if ioutput > 2900: z_sel = 200 / len_AU
    c_center = center_sink_list[isam,:]
    vec_z = ax_list[isam,:]
    #vec_z = np.array([0,0,-1]); 
    #r_sel*=0.9; r_cyl *= 0.9; z_sel *=0.9
    vec_z = vec_z /np.linalg.norm(vec_z)
    print 'vec_z', vec_z
    cvel = cvel_sink_list[isam,:]
    filecell = path_out+'/selected_cells_'+str(ioutput).zfill(5)+'.npz'
    if os.path.exists(filecell):
        cell_disc = np.load(filecell)
        dV_disc = cell_disc["dx"]**3
        dS_disc = cell_disc["dx"]**2
        dx_disc = cell_disc["dx"]
        posi_disc = cell_disc["x"]-c_center[np.newaxis,:]
        B_disc = cell_disc["B"]
        T_disc = cell_disc["T"]
    else:
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        amr = ro.amr_source(["rho","vel","P","phi","g","Br","Bl","T"])
        disc_region = Cylinder(c_center,vec_z,r_cyl,z_sel*2)
        amr_disc = RegionFilter(disc_region,amr)
        cell_disc = CellsToPoints(amr_disc, smallest_cell_level=scl).flatten()
        B_disc = 0.5 * (cell_disc["Br"] + cell_disc["Bl"])
        T_disc = cell_disc["T"]
        np.savez(filecell,rho=cell_disc["rho"],vel=cell_disc["vel"],P=cell_disc["P"],g=cell_disc["g"],B=B_disc,x=cell_disc.points,dx=cell_disc.get_sizes(),T=T_disc)
        dx_disc = cell_disc.get_sizes()
        dV_disc = dx_disc**3
        dS_disc = dx_disc**2
        posi_disc = cell_disc.points-c_center[np.newaxis,:]

    time = time_list[isam]
    thres_rho = thres_rho_list[isam]
    mass_sink = Msink[isam] / units["mass_sol"]   
    dm = cell_disc["rho"] * dV_disc

#### define grids and sums ####
    zd = np.sum(posi_disc * vec_z[np.newaxis,:],axis=1)
    rsph = np.linalg.norm(posi_disc,axis=1)
    rd = np.sqrt(rsph**2-zd**2)
    print 'lmax', units["lmax"]
    ndx0 = max(1,2**(units["lmax"]-17)) ## smoothing factor
    ndx1 = max(1,2**(units["lmax"]-16)) ## binning factor
    ndx = ndx1 / ndx0
    print 'smoothing=',ndx0, 'binninb=',ndx1, 'ratio=',ndx
    dx0 = 1./2**units["lmax"]
    dx = dx * ndx0 #/  abs(vec_z[2]) * 2
    rr_cyl = np.arange(0,r_cyl,dx)
#    rr = np.arange(0,r_sel,dx)
    rr = rr_cyl[:-npad]
    zz = np.arange(0,z_sel,dx)    
    dd = [0, thres_rho*2, 1.e40]
    rc = 0.5 * ( rr[:-1]+rr[1:] )
    zc = 0.5 * ( zz[:-1]+zz[1:] )
    zcs = np.append(-zc[-1::-1],zc)
    nr = len(rc); nz=len(zc)
    surf  = np.pi * np.diff(rr**2)
    hist_x = lambda x: np.histogramdd((rd,cell_disc["rho"]), bins=(rr,dd), weights=x)[0]
    histr_x = lambda x: np.histogram(rsph,bins=rr,weights=x)[0]
    grid_x = lambda x: np.histogramdd((rd,abs(zd)), bins=(rr_cyl,zz), weights=x)[0]
    grid_x_c = lambda x: np.cumsum(np.histogramdd((rd,abs(zd)), bins=(rr,zz), weights=x)[0],axis=1)
    cumz = lambda x: np.cumsum(x,axis=1)
    zpm = np.append(-zz[-1:0:-1],zz)
    zcpm = np.append(-zc[-1::-1],zc)
    grid_xpm = lambda x: np.histogramdd((rd,zd), bins=(rr_cyl,zpm), weights=x)[0]
    def grid_lev_x(x,func_grid=grid_x,dims=[nr,nz],wdim=3,symm=[1,1]): #wdim: dimension of the weight (3 or 2)
        level = np.log(dx_disc/dx)/np.log(2)
        lm = int(np.max(level))
        #lmin = -int(np.log(ndx)/np.log(2))
        print 'level max',lm
#        inds = np.where(abs(level-0)>0.5)[0]
        inds = np.where(level>0)[0]
        y = x.copy(); y[inds]=0.
        hists = func_grid(y)
        ndim = len(dims)
        r2 = np.ones([2]*ndim)/2**(ndim)
        for il in range(1,lm+1):
            inds = np.where(abs(level-il)>0.5)[0]
            y = x.copy(); y[inds]=0.
            hist_i = func_grid(y)
            conv = np.ones([2**il]*ndim)/(2**il)**(ndim-(3-wdim))
            conv = sp.signal.fftconvolve(conv,r2)
            hist_i = sp.signal.fftconvolve(hist_i,conv)
            nmar = 2**(il-1)
            if ndim == 2:
                 hist_i[nmar:nmar*2,:] += symm[0]*hist_i[nmar-1::-1,:]
                 hist_i[:,nmar:nmar*2] += symm[1]*hist_i[:,nmar-1::-1]
                 hists = hists + hist_i[nmar:-nmar,nmar:-nmar]
            if ndim ==1:
                 hist_i[nmar:nmar*2] += symm[0]*hist_i[nmar-1::-1]
                 hists = hists + hist_i[nmar:-nmar]
        return hists[:-npad]
    ## usage for full plane: grid_lev_x(x, func_grid=grid_xpm dims=[nr,2*nz],symm=[1,0])

#### variables inside disk (and envelope) for each cell ####
    vel_disc = cell_disc["vel"]-cvel[np.newaxis,:]
    vecrsph = posi_disc / np.linalg.norm(posi_disc,axis=1)[:,np.newaxis]
    vecphi = np.cross(vec_z, vecrsph) # phi vector in cylindrical coordinate 
    vecphi = vecphi / np.linalg.norm(vecphi,axis=1)[:,np.newaxis]
    vecr = np.cross(vecphi, vec_z)
    vecr = vecr / np.linalg.norm(vecr,axis=1)[:,np.newaxis] # r vector in cylindrical coordinate 

    vrsph = np.sum( vel_disc * vecrsph, axis=1)
    vr = np.sum( vel_disc * vecr, axis=1)
    vphi = np.sum( vel_disc * vecphi, axis=1)
    vz = np.sum( vel_disc * vec_z, axis=1)
    B_r = np.sum(B_disc * vecr, axis=1)
    B_phi = np.sum(B_disc * vecphi, axis=1)
    B_z = np.sum(B_disc * vec_z, axis=1)
    g_r = np.sum(cell_disc["g"] * vecr, axis=1)
    g_phi = np.sum(cell_disc["g"] * vecphi, axis=1)
    g_z = np.sum(cell_disc["g"] * vec_z, axis=1)

    tick_size = 16
    fontlabel_size = 16
    legendlabel_size = 12
    axislabel_size = 16
    params = {'backend': 'wxAgg', 'lines.markersize' : 3, 'axes.labelsize': axislabel_size, 'font.size': fontlabel_size, 'legend.fontsize': legendlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'axes.linewidth' : 2., 'legend.labelspacing': 0.2, 'legend.handletextpad':0.1}
    plt.rcParams.update(params)

#### properties along r ####
    mass = hist_x(dm)
    mzz = hist_x(dm*zd**2)
    mz = hist_x(dm*abs(zd))
    mcc = hist_x(dV_disc * cell_disc["P"]) # * 1.5 
    mass_z = np.sum(mass,axis=1) ## mass integrated in z
    mass_r = np.cumsum(mass_z) ## cumulated mass in r
    hd_1 = np.sqrt(np.pi/2.) * np.sum(mz,axis=1) / mass_z
    hd_2 = np.sqrt(np.sum(mzz,axis=1)/mass_z)
    hd_1[mass_z==0]=0.; hd_2[mass_z==0]=0.
    hd_diff = abs(hd_1/hd_2-1.)  
    hd = hd_1.copy() ## disc scale height from vertical mass distribution
    ind_1AU = np.where(hd[5:]*len_AU>1.)[0][0]+5 # index where h>1AU
#    cs = np.sqrt(mcc[:,1]/mass[:,1]) ## sound speed selecting only dense cells
#    hp = rc**1.5 * cs / boxlen / np.sqrt(mass_sink)  ## disc scale height from hydrostatic equilibrium
    cs_all = np.sqrt(np.sum(mcc,axis=1)/mass_z)
    hp_all = rc**1.5 * cs_all / boxlen / np.sqrt(mass_sink)
    hp_all2 = rc**1.5 * cs_all / boxlen / np.sqrt(mass_sink+mass_r) #take into account disk mass
    hs = np.vstack((hd_1,hd_2,hp_all2))
    hd_var = np.std(hs[:2,:],axis=0) / np.mean(hs[:2,:],axis=0)
    h_var = np.std(hs,axis=0) / np.mean(hs,axis=0)
    ind_rmax = np.where(h_var[ind_1AU:]>2.e-2)[0][0]+ind_1AU
    print 'Rd=',rc[ind_rmax]*len_AU, ind_rmax

    def plot_reso(f,ylim,dx0=dx0,a=0.1): ## add resolution limit to plot
        plt.figure(f.number)
        plt.plot(np.ones(2)*dx0*len_AU, ylim,'k',zorder=0,alpha=a)
        plt.plot(np.ones(2)*4*dx0*len_AU, ylim,'k',zorder=0,alpha=a)

    f_h = False
    f_Sigma = False
    f_Md = False
    f_vel = False
    f_beta = False
    f_Q = False
    f_alpha = False
    f_alpha_bins = False
    f_T = True
    compare=True
    f_jflux_imshow = False

    f_Mdot = False
    f_S = False
    f_flux = False
    f_Mdot_imshow = False
    f_Mdot_layer = False
    f_vertical = False
    f_alpha_imshow = False
    f_pm = False
    f_Mdot_alpha = False
    overplot_S = False

##    f_pm = True
#    f_h = True
#    f_Sigma = True
#    f_Md = True
#    f_vel = True
#    f_beta = True
#    f_Q = True
#    f_T = True
#    f_S = True
#    overplot_S = True
#    f_Mdot = True
#    f_Mdot_alpha = True
#    f_alpha_imshow = True
#    f_jflux_imshow = True
#    f_alpha_bins = True
###    f_flux = True
#    f_Mdot_imshow = True
#    f_Mdot_layer = True
##    f_vertical = True

    if f_h: ## plot disk scale height from several definitions 
        ylim_h = [1.e-1,1.e2]
        f_h = plt.figure(figsize=(7,3))
        plt.plot(rc*len_AU,hp_all*len_AU, label=r'$h_p$')
        plt.plot(rc*len_AU,hp_all2*len_AU,c='C0',alpha=0.5)
        plt.plot(rc*len_AU,hd_1*len_AU, label=r'$h_1$')
        plt.plot(rc*len_AU,hd_2*len_AU, label=r'$h_2$')
        plot_reso(f_h,ylim_h)
        plt.plot([1.e-2,1.e3], np.ones(2)*dx0*len_AU,'k',zorder=0,alpha=0.1)
        plt.xscale('log'); plt.yscale('log')
        plt.xlabel('r(AU)');plt.ylabel('H(AU)')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_h)
        l_g = plt.legend()
        plt.savefig(path_out+'/disk_height_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

#### properties on r-z grid ####
#    mvrvphi_d = grid_x( dm*vr*vphi ); mvrvphi_c = cumz( mvrvphi_d )
#    mvr_d = grid_x( dm*vr );  mvr_c = cumz( mvr_d)
#    mvphi_d = grid_x( dm*vphi); mvphi_c = cumz( mvphi_d )
#    mvz_d = grid_x( dm*vz*np.sign(zd))
##    rhoSvr_d = grid_x( cell_disc["rho"]*dS_disc*vr)
##    rhoSvz_d = grid_x( cell_disc["rho"]*dS_disc*vz*np.sign(zd))
#    mass_d = grid_x( dm ); mass_c = cumz( mass_d )
#    vol_d = grid_x( dV_disc )
##    sur_d = grid_x( dS_disc )
#    mcc_d = grid_x( dV_disc * cell_disc["P"]) ; mcc_c = cumz( mcc_d )
#    mdvrdvphi_c = mvrvphi_c - mvr_c * mvphi_c / mass_c
#    VBrBphi_d = grid_x(dV_disc * B_r * B_phi); VBrBphi_c = cumz( VBrBphi_d )
#    VBzBphi_d = grid_x(dV_disc * B_z * B_phi * np.sign(zd))
#    Vgrgphi_d = grid_x(dV_disc * g_r * g_phi)/ (4.*np.pi); Vgrgphi_c = cumz( Vgrgphi_d )
#    mT_d = grid_x(dm * T_disc); mT_c = cumz(mT_d)
#    T_d = mT_d / mass_d; T_c = mT_c/mass_c
##    mz_c = grid_x_c( dm * abs(zd) ) 
#    VBB_d = grid_x(dV_disc * np.sum(B_disc**2,axis=1)); VBB_c = cumz(VBB_d)

    signz = np.sign(zd)
    signsym = 1.
#    signz = np.ones_like(zd)
    #signz = (signz-1)*0.5
    #signsym = -signz
    mvrvphi_l = grid_lev_x( dm*vr*vphi*signsym); mvrvphi_c = cumz( mvrvphi_l )
    mvzvphi_l = grid_lev_x( dm*vz*vphi * signz, symm=[-1,-1])#; mvzvphi_c = cumz( mvzvphi_l )
    mvr_l = grid_lev_x( dm*vr*signsym , symm=[-1,1]); mvr_c = cumz( mvr_l)
    mvphi_l = grid_lev_x( dm*vphi*signsym, symm=[-1,1] ); mvphi_c = cumz( mvphi_l )
    mvz_l = grid_lev_x( dm*vz*signz, symm=[1,-1])
    mass_l = grid_lev_x( dm*signsym );  mass_c = cumz( mass_l )
    vol_l = grid_lev_x(dV_disc*signsym); vol_c = cumz( vol_l ) 
    mcc_l = grid_lev_x( dV_disc * cell_disc["P"]*signsym) ; mcc_c = cumz( mcc_l )
    mdvrdvphi_l = mvrvphi_l - mvr_l * mvphi_l / mass_l
    mdvrdvphi_c = mvrvphi_c - mvr_c * mvphi_c / mass_c
    mdvzdvphi_l = mvzvphi_l - mvz_l * mvphi_l / mass_l

    VBrBphi_l = -grid_lev_x(dV_disc * B_r * B_phi*signsym); VBrBphi_c = cumz( VBrBphi_l )
    VBzBphi_l = -grid_lev_x(dV_disc * B_z * B_phi * signz, symm=[-1,-1])
    Vgrgphi_l = grid_lev_x(dV_disc * g_r * g_phi*signsym)/ (4.*np.pi); Vgrgphi_c = cumz( Vgrgphi_l )
    Vgzgphi_l = grid_lev_x(dV_disc * g_z * g_phi * signz, symm=[-1,-1])/ (4.*np.pi)
    mT_l = grid_lev_x(dm * T_disc*signsym); mT_c = cumz(mT_l)
    T_l = mT_l / mass_l; T_c = mT_c/mass_c
    VBB_l = grid_lev_x(dV_disc * np.sum(B_disc**2,axis=1)*signsym); VBB_c = cumz(VBB_l)
    mass_r = np.cumsum(mass_c[:,-1])
    VBr_l = grid_lev_x(dV_disc * B_r * signz, symm=[-1,-1]) 
    VBphi_l = grid_lev_x(dV_disc * B_phi * signz, symm=[-1,-1])
    VBz_l = grid_lev_x(dV_disc * B_z*signsym)

#### surface density ####
    Sigma_sel = mass_c[:,-1] / surf
    Sigma_cor = Sigma_sel / sp.special.erf(zc[-1]/hd/np.sqrt(2.))
    if f_Sigma: ## plot disk surface density (measured and corrected)
        ylim_Sig = [1.e0,5.e4]
        ylim_Sig = [1.e-2,5.e4]
        f_Sigma = plt.figure(figsize=(7,3))
        plt.loglog(rc*len_AU,Sigma_sel*kg_m2,label=r'$\Sigma_{cyl}$')
#        plt.loglog(rc*len_AU,Sigma_cor*kg_m2,c='C0',alpha=0.5) #label=r'$\Sigma_{tot}$') 
        plot_reso(f_Sigma,ylim_Sig)
        plt.legend()
        plt.xlim([xmin,xmax]); plt.ylim(ylim_Sig)
        plt.xlabel('r(AU)');plt.ylabel(r'$\Sigma(kg/m^2)$')
        plt.savefig(path_out+'/surface_density_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#    S_crit = 100. / kg_m2 #kg/m^2 disk height at given critical density of ionization
#    z_crit = np.sqrt(2) * hd * sp.special.erfinv(np.maximum(1.-2*S_crit/Sigma_cor,0.))
#    if f_h and False:  # plot z_crit onto disk scale height
#        plt.figure(f_h.number)
#        plt.plot(rc*len_AU, z_crit*len_AU,label='hcrit')
#        plt.legend()
#        plt.savefig(path_out+'/disk_height_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

#### values at selected altitudes ####
    ind_end = nr
    #h_disk = hp_all2
    fac2 = 3; nfloor = 4
    def sel_height(x,h=hp_all2,cum=True,cen=False,flo=False): #cum: cumulated, cen: centered, flo: floored
        #h = h_disk
        xh = np.zeros(ind_end)
        x2 = np.zeros(ind_end)
        for ir in range(ind_end):
          xloc = x[ir,:].copy()
          if flo: xloc[:flo] = xloc[flo] #floored
          if cum: func_x_z = interp1d(zz,np.append(0,xloc)) #if is an accumulated function
          if cen: func_x_z = interp1d(np.append(0,zc),np.append(xloc[0],xloc)) # if is centered on cell
          if cum and flo: func_x_z = interp1d(zz,np.append(xloc[0],xloc))
          if h[ir]>0.:
            xh[ir] = func_x_z(np.minimum(h[ir],zz[-2]))
            x2[ir] = func_x_z(np.minimum(h[ir]*fac2,zz[-2]))
        return [xh, x2]
#### values at selected altitudes ####
#    ind_end = nr
#    h_disk = hp_all2
##    fac2 = 3; nfloor = 4
#    def sel_height(x,h=hp_all2,fac2=3,cum=True,cen=False,flo=False): #cum: cumulated, cen: centered, flo: floored
##        h = h_disk
#        xh = np.zeros(ind_end)
#        x2 = np.zeros(ind_end)
#        xc = np.zeros(ind_end)
#        xd = np.zeros(ind_end)
#        for ir in range(ind_end):
#          xloc = x[ir,:].copy()
##          if ir > 17 and ir<21: print ir, xloc[5:10]
#          if flo: xloc[:flo] = xloc[flo] #floored
#          if cum: func_x_z = interp1d(zz,np.append(0,xloc)) #if is an accumulated function
#          if cen: func_x_z = interp1d(np.append(0,zc),np.append(xloc[0],xloc)) # if is centered on cell
#          if cum and flo: func_x_z = interp1d(zz,np.append(xloc[0],xloc))
#          if h[ir]>0.:
#            xh[ir] = func_x_z(np.minimum(h[ir],zz[-2]))
#            x2[ir] = func_x_z(np.minimum(h[ir]*fac2,zz[-2]))
##          xc[ir] = func_x_z(z_crit[ir])
##          if z_crit[ir]==0.: xc[ir]=0.
##          if h_disk[ir]*fac2 > z_crit[ir]: xd[ir] = x2[ir]-xc[ir]
##          else: xd[ir] = 0.
#        return [xh, x2, xc, xd]

#### deviation from keplerian rotation as first radius estimation ####
#    MM_l = sel_height(mass_l)
    MM_c = sel_height(mass_c)
    MVp_c = sel_height(mvphi_c)
    MVr_c = sel_height(mvr_c)
    Vp_c = map(truediv, MVp_c, MM_c)
    Vr_c = map(truediv, MVr_c, MM_c)
    #vphi_kep = np.sqrt( mass_sink / rc ) * boxlen 
    vphi_kep2 = np.sqrt( (mass_sink+mass_r)/rc ) * boxlen
    vphi_dev = abs(1.-Vp_c[1]/MM_c[1]/vphi_kep2[:ind_end])
    ind_rrot = np.where(vphi_dev[int(ind_rmax*0.6):]>0.02)[0][0]+int(ind_rmax*0.6)

#### keplerian disk and pseudodisk radius ####    
    Scrit = Sigma_cor[ind_rrot]
    Rcrit = rc[ind_rrot]
    sa = -np.sum(np.log(Sigma_cor[:ind_rrot]/Scrit)*np.log(rc[:ind_rrot]/Rcrit)) / np.sum(np.log(rc[:ind_rrot]/Rcrit)**2)
    sb = -np.sum(np.log(Sigma_cor[ind_rrot+1:-npad]/Scrit)*np.log(rc[ind_rrot+1:-npad]/Rcrit)) / np.sum(np.log(rc[ind_rrot+1:-npad]/Rcrit)**2)
    func_Sig = lambda rr,s,r,a,b,c,d: np.log(s * (1./((rr/r)**a+(rr/r)**b) + d*(r/rr)**c))
    #Sig_params = curve_fit(func_Sig,rc[4:],np.log(Sigma_sel[4:]),p0=(Scrit,Rcrit*1.2,sa,sb,1,5.e-3))[0]
#    plt.figure()
#    plt.loglog(rc[4:-npad],Sigma_sel[4:-npad])
#    plt.show()
    print Sigma_sel[4:-npad]
    Sig_params = curve_fit(func_Sig,rc[4:-npad],np.log(Sigma_sel[4:-npad]),p0=(Scrit,Rcrit*1.,sa,sb,1,1.e-2))[0]
#    Sig_params = curve_fit(func_Sig,rc[4:-npad],np.log(Sigma_sel[4:-npad]),p0=(Scrit,Rcrit*1.,sa,sb,1,5.e-3))[0]
    print Scrit*kg_m2,Rcrit*len_AU,sa,sb
    print Sig_params[0]*kg_m2,Sig_params[1]*len_AU,Sig_params[1]*len_AU*Sig_params[5]**(1./(Sig_params[4]-Sig_params[3])),Sig_params[2],Sig_params[3],Sig_params[4],Sig_params[5]
    Sigma_par = np.exp(func_Sig(rc,Sig_params[0],Sig_params[1],Sig_params[2],Sig_params[3],Sig_params[4],Sig_params[5]))
    Rk = Sig_params[1]; Rp = Rk*Sig_params[5]**(1./(Sig_params[4]-Sig_params[3]))
    func_mr = interp1d(rc,mass_r)
    print 'radius', Rk*len_AU, Rp*len_AU
    mass_Rk = func_mr(Rk); mass_Rp = func_mr(Rp)
    def plot_Rk_Rp(f, ylim,Rk=Rk,Rp=Rp,a=1):
        plt.figure(f.number)
        plt.plot(np.ones(2)*Rk*len_AU, ylim,ls=':',c='C8',zorder=1,alpha=a)
        plt.plot(np.ones(2)*Rp*len_AU, ylim,ls=':',c='C8',zorder=1,alpha=a)
    if f_Sigma:
        plt.figure(f_Sigma.number)
        plt.plot(rc[4:]*len_AU,Sigma_par[4:]*kg_m2,label=r'$\Sigma \propto r^{%.2f,%.2f,%.2f}$'%(-Sig_params[2],-Sig_params[3],-Sig_params[4]))
        plot_Rk_Rp(f_Sigma,ylim_Sig)
        plt.legend()
        plt.text(xmin*1.1,ylim_Sig[1]*2.e-3,r'$M_\ast=%.3f (M_\odot)$''\n'r'$M_d=%.4f+%.4f (M_\odot)$'%(Msink[isam],mass_Rk*units["mass_sol"], (mass_Rp-mass_Rk)*units["mass_sol"]))
        plt.savefig(path_out+'/surface_density_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    if f_h: 
        plot_Rk_Rp(f_h,ylim_h)
        plt.savefig(path_out+'/disk_height_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    ind_kep = np.where(rr>Rk)[0][0]-1
    ind_pse = np.where(rr>Rp)[0][0]-1
    #plt.show()
 
#### disk mass #### 
    if f_Md: ## plot disk surface density (measured and corrected)
        ylim_Md = [1.e-4,1.e-1]
        f_Md = plt.figure(figsize=(7,3))
        #f_Md, (ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(7,3))
        plt.loglog(rc*len_AU,mass_r*units["mass_sol"])
        #ax1.loglog([xmin,xmax],np.ones(2)*Msink[isam],'k--')
        plot_reso(f_Md,ylim_Md)
        plot_Rk_Rp(f_Md,ylim_Md)
        #plt.legend()
        plt.xlim([xmin,xmax]); plt.ylim(ylim_Md)
        #ax2.set_ylim(ylim_Md[0],2.e-2); ax1.set_ylim(0.1,ylim_Md[1])
        #ax1.spines['bottom'].set_visible(False)
        #ax2.spines['top'].set_visible(False)
        #ax1.xaxis.tick_top()
        #ax1.tick_params(labeltop='off')  # don't put tick labels at the top
        #ax2.xaxis.tick_bottom()
        #d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        #kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        #ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        #ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
        #kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        #ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        #ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  
        plt.xlabel('r(AU)');plt.ylabel(r'$M_{d}(M_\odot)$')
        plt.savefig(path_out+'/disk_m_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    #plt.show()
### smooth fit of hp
    wght = rc**0
    ind0 = 4
    sum_hr = np.sum(np.log(hp_all2[ind0:ind_pse])*np.log(rc[ind0:ind_pse])*wght[ind0:ind_pse])
    sum_1 = np.sum(wght[ind0:ind_pse])
    sum_h = np.sum(np.log(hp_all2[ind0:ind_pse])*wght[ind0:ind_pse])
    sum_r = np.sum(np.log(rc[ind0:ind_pse])*wght[ind0:ind_pse])
    sum_r2 = np.sum(np.log(rc[ind0:ind_pse])**2*wght[ind0:ind_pse])
    deno = sum_r2*sum_1 - sum_r**2
    exp_h = (sum_hr*sum_1-sum_h*sum_r)/deno
    fac_h = (sum_h*sum_r2-sum_hr*sum_r)/deno
    fac_h1 = fac_h + np.log(len_AU)*(1-exp_h)
    h_disk = np.exp(fac_h)*rc**exp_h
    print h_disk[ind0:ind_pse] / hp_all2[ind0:ind_pse], fac_h1
    if f_h:
        plt.figure(f_h.number) 
        plt.plot(rc[ind0:ind_pse]*len_AU,h_disk[ind0:ind_pse]*len_AU,'k--',alpha=0.4)
        l_g.get_texts()[0].set_text(r'$h_p \approx %.2f \/ r^{%.2f}$'%(np.exp(fac_h1),exp_h))
        plt.savefig(path_out+'/disk_height_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#    plt.show()


#### values at selected altitudes ####
    ind_end = nr
#    fac2 = 3; nfloor = 4
    def sel_height(x,h=h_disk,cum=True,cen=False,flo=False): #cum: cumulated, cen: centered, flo: floored
        xh = np.zeros(ind_end)
        x2 = np.zeros(ind_end)
#        xc = np.zeros(ind_end)
#        xd = np.zeros(ind_end)
        for ir in range(ind_end):
          xloc = x[ir,:].copy()
          if flo: xloc[:flo] = xloc[flo] #floored
          if cum: func_x_z = interp1d(zz,np.append(0,xloc)) #if is an accumulated function
          if cen: func_x_z = interp1d(np.append(0,zc),np.append(xloc[0],xloc)) # if is centered on cell
          if cum and flo: func_x_z = interp1d(zz,np.append(xloc[0],xloc))
          if h[ir]>0.:
            xh[ir] = func_x_z(np.minimum(h[ir],zz[-2]))
            x2[ir] = func_x_z(np.minimum(h[ir]*fac2,zz[-2]))
#          xc[ir] = func_x_z(z_crit[ir])
#          if z_crit[ir]==0.: xc[ir]=0.
#          if h_disk[ir]*fac2 > z_crit[ir]: xd[ir] = x2[ir]-xc[ir]
#          else: xd[ir] = 0.
        return [xh, x2] #, xc, xd]

#### velocity properties ####
    MM_l = sel_height(mass_l,cum=False,cen=True)
    MM_c = sel_height(mass_c)
    MVp_c = sel_height(mvphi_c)
    MVr_c = sel_height(mvr_c)
    Vp_c = map(truediv, MVp_c, MM_c)
    Vr_c = map(truediv, MVr_c, MM_c)
    vphi_kep = np.sqrt( mass_sink / rc ) * boxlen
    vphi_kep2 = np.sqrt( (mass_sink+mass_r)/rc ) * boxlen
    vphi_dev = abs(1.-Vp_c[1]/MM_c[1]/vphi_kep[:ind_end])
#    ind_rrot = np.where(vphi_dev[int(ind_rmax*0.6):]>0.02)[0][0]+int(ind_rmax*0.6)
    if f_vel: ## plot velocity profiles
        ylim_vel = [-1.e3,5.e4]
        f_vel = plt.figure(figsize=(7,3))
        plt.loglog(rc*len_AU,vphi_kep * vel_ms, label=r'$u_{kep}$')
        plt.loglog(rc*len_AU,vphi_kep2 * vel_ms, c='C0',alpha=0.5)
        plt.plot(rc[:ind_end]*len_AU,Vp_c[0] * vel_ms, label=r'$u_\phi(H)$')
        plt.plot(rc[:ind_end]*len_AU,Vp_c[1] * vel_ms, label=r'$u_\phi(%i H)$'%fac2)
        plt.plot(rc[:ind_end]*len_AU,-Vr_c[0] * vel_ms,c='C3', label=r'$-u_r(H)$')
        plt.plot(rc[:ind_end]*len_AU,-Vr_c[1] * vel_ms,c='C4', label=r'$-u_r(%i H)$'%fac2)
    #    plt.plot(rc[:ind_end]*len_AU,Vr_c[0] * vel_ms,c='C3', ls=':')
    #    plt.plot(rc[:ind_end]*len_AU,Vr_c[1] * vel_ms,c='C4', ls=':')
        plot_reso(f_vel,ylim_vel)
        plot_Rk_Rp(f_vel,ylim_vel)
        plt.xscale('log'); #plt.yscale('log')
        plt.yscale('symlog',linthreshy = 1.e2)
        plt.xlabel('r(AU)');plt.ylabel(r'$u$(m/s)')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_vel)
        plt.legend()
        plt.savefig(path_out+'/disk_velocity_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

#### mass accretion in layers ####
    r_img = int(Rp*1.2*len_AU/10)*10; z_img = int(r_img*0.4/5)*5
    if f_Mdot_imshow:
        #r_img = 40; z_img = 15
        vr_l = mvr_l/mass_l
        vz_l = mvz_l/mass_l
        rhovr_l = mvr_l/vol_l / boxlen
        rhovz_l = mvz_l/vol_l / boxlen
        colornorm = lambda lim1,lim2: SymLogNorm(10**lim1, vmin = -10**lim2, vmax = 10**lim2)
        rlast=int(r_img/len_AU/dx); zlast = int(z_img/len_AU/dx)
    #    rlast=int(20/len_AU/dx); zlast = int(5/len_AU/dx)
        imgext = np.asarray([rr[0],rr[rlast],zz[0],zz[zlast]])*len_AU
        rgrid, zgrid = np.meshgrid(rr[:rlast]*len_AU,zz[:zlast]*len_AU)
        vres = max(int(Rp/40./dx),1)
    
        vnorm  = np.sqrt(vr_l**2+vz_l**2)
        vr_dir = vr_l[vres/2:rlast:vres,vres/2:zlast:vres]/vnorm[vres/2:rlast:vres,vres/2:zlast:vres]
        vz_dir = vz_l[vres/2:rlast:vres,vres/2:zlast:vres]/vnorm[vres/2:rlast:vres,vres/2:zlast:vres]
        rhovnorm = vnorm * mass_l / vol_l / boxlen  
     
        ind_end_h = ind_kep
        plt.figure(figsize=(7,3))
        plt.imshow((vr_l[:rlast,:zlast]*vel_ms).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm(0,5),aspect='equal')
        plt.colorbar(extend='both')
        plt.contour(rgrid,zgrid,vr_l[:rlast,:zlast].T,[0.])
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.quiver(rgrid[vres/2::vres,vres/2::vres],zgrid[vres/2::vres,vres/2::vres],vr_dir.T,vz_dir.T,pivot='middle',headlength=2,headwidth=2,color='w')
        plt.title(r'$v_r (m/s, t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/vr_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        plt.figure(figsize=(7,3))
#        plt.imshow((mvr_l[:rlast,:zlast]/(dx*boxlen)*Ms_Myr*1.e-6).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm(-8,-5),aspect='equal')
        plt.imshow((rhovr_l[:rlast,:zlast]*kg_m2_s).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm(-10,-6),aspect='equal')
        plt.colorbar(extend='both')
        plt.contour(rgrid,zgrid,mvr_l[:rlast,:zlast].T,[0.])
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.quiver(rgrid[vres/2::vres,vres/2::vres],zgrid[vres/2::vres,vres/2::vres],vr_dir.T,vz_dir.T,pivot='middle',headlength=2,headwidth=2,color='w')
        plt.title(r'$Flux_r (kg/m^2/s, t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.title(r'$\dot{M}_r (M_\odot/yr, t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/Fluxr_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        plt.figure(figsize=(7,3))
        plt.imshow((vz_l[:rlast,:zlast]*vel_ms).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm(0,5),aspect='equal')
        plt.colorbar(extend='both')
        plt.contour(rgrid,zgrid,vz_l[:rlast,:zlast].T,[0.])
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.quiver(rgrid[vres/2::vres,vres/2::vres],zgrid[vres/2::vres,vres/2::vres],vr_dir.T,vz_dir.T,pivot='middle',headlength=2,headwidth=2,color='w')
        plt.title(r'$v_z (m/s, t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/vz_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        plt.figure(figsize=(7,3))
#        plt.imshow((mvz_l[:rlast,:zlast]/(dx*boxlen)*Ms_Myr*1.e-6).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm(-8,-5),aspect='equal')
        plt.imshow((rhovz_l[:rlast,:zlast]*kg_m2_s).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm(-10,-6),aspect='equal')
        plt.colorbar(extend='both')
        plt.contour(rgrid,zgrid,mvz_l[:rlast,:zlast].T,[0.])
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.quiver(rgrid[vres/2::vres,vres/2::vres],zgrid[vres/2::vres,vres/2::vres],vr_dir.T,vz_dir.T,pivot='middle',headlength=2,headwidth=2,color='w')
        plt.title(r'$Flux_z (kg/m^2/s, t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/Fluxz_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#    plt.show()

#### some disk radial internal properties #### 
    PP_l = sel_height(mcc_l,cum=False,cen=True)
    PP_c = sel_height(mcc_c)
    CC_l = map(truediv, PP_l, MM_l); CC_c = map(truediv, PP_c, MM_c)
    BB_c = sel_height(VBB_c/2.)
    MT_c = sel_height(mT_c)

    if f_beta: ## plot plasma beta
        ylim_beta = [1.e-1,1.e4]
        f_beta = plt.figure(figsize=(7,3))
        be = PP_c[0]/BB_c[0]; be2 = PP_c[1]/BB_c[1] 
        plt.loglog(rc*len_AU, be,'C1',label=r'$H$')
        plt.plot(rc*len_AU, be2,'C2',label=r'$%i H$'%fac2)
        plot_reso(f_beta,ylim_beta)
        plot_Rk_Rp(f_beta,ylim_beta)
        plt.xlabel('r(AU)');plt.ylabel(r'$\beta$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_beta)
        plt.legend()
        plt.savefig(path_out+'/disk_beta_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    if f_Q: ## plot Toomre Q parameter
        ylim_Q = [1.e0,1.e3]
        f_Q = plt.figure(figsize=(7,3))
        Q_n = np.sqrt(mcc_c) * mvphi_c 
        Q_d = mass_c**(2.5)
        QQ_n = sel_height(Q_n)  
        QQ_d = sel_height(Q_d)
        QQ = map(truediv, QQ_n, QQ_d)
        QQ = [Q/(rc/surf*np.pi*boxlen**2) for Q in QQ]
        plt.loglog(rc*len_AU, QQ[0],'C1',label=r'$H$')
        plt.plot(rc*len_AU, QQ[1],'C2',label=r'$%i H$'%fac2)
        plot_reso(f_Q,ylim_Q)
        plot_Rk_Rp(f_Q,ylim_Q)
        plt.xlabel('r(AU)');plt.ylabel(r'$Q$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_Q)
        plt.legend()
        plt.savefig(path_out+'/disk_Q_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    def calc_snowlines(func_r_T, Ts):
        Rsnow = []
        for T in Ts[:-1]:
            if T>Ts[-1]:
                Rsnow.append(0.); continue
            r = func_r_T(T); Rsnow.append(r)
        return Rsnow
    def plot_snowlines(f, ylim, Ts,Rsnow,a=1):
        plt.figure(f.number)
        for T, r in zip(Ts,Rsnow):
            if r==0.: continue
            plt.plot(np.ones(2)*r*len_AU, ylim,ls='--',c='C9',lw=1,zorder=1,alpha=a)
            plt.plot(np.array([xmin,xmax]),np.ones(2)*T,ls='--',c='C9',lw=1,zorder=1,alpha=a)
    TT_c = map(truediv, MT_c, MM_c)
    func_r_T = interp1d(TT_c[1],rc)
    Tcond = [150,1500,np.max(TT_c[1])]
    Rsnow = calc_snowlines(func_r_T,Tcond)
    if f_T: ## plot temperature
        ylim_T = [1.e1,1.e4]
        f_T = plt.figure(figsize=(7,3))
        plt.loglog(rc*len_AU, TT_c[0],'C1',label=r'$H$')
        plt.plot(rc*len_AU, TT_c[1],'C2',label=r'$%i H$'%fac2)
        plot_reso(f_T,ylim_T)
        plot_Rk_Rp(f_T,ylim_T)
        plot_snowlines(f_T,ylim_T,Tcond,Rsnow)
        plt.xlabel('r(AU)');plt.ylabel(r'$T(K)$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_T)
        plt.legend()
        np.savez(path_out+'/disk_T_'+str(ioutput).zfill(5)+'.npz',TT_c=TT_c, rc=rc, Rk=Rk,Rp=Rp,Rsnow=Rsnow,dx0=dx0)
        if compare:
            tc = np.load(path_out+'_nofeed/disk_T_'+str(ioutput).zfill(5)+'.npz')
            plt.loglog(rc*len_AU, tc["TT_c"][0],'C1',alpha=0.4)
            plt.plot(rc*len_AU, tc["TT_c"][1],'C2',alpha=0.4)
            plot_Rk_Rp(f_T,ylim_T,Rk=tc["Rk"],Rp=tc["Rp"],a=0.4)
            plot_snowlines(f_T,ylim_T,Tcond,tc["Rsnow"],a=0.4)
            plt.savefig(path_out+'/disk_T_comp_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
       # else: plt.savefig(path_out+'/disk_T_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    if f_vertical: ## plot vertical profiles
        R_sel = [1,5,10]
        R_ind = [int(R/dx) for R in R_sel]
        z_lim = [0,5]
        f_vert_T = plt.figure(figsize=(7,3))
        f_vert_rho = plt.figure(figsize=(7,3))
        f_vert_P = plt.figure(figsize=(7,3))

### measure source function and mass transport ###
    nfloor=0
#    Sr_d = mvr_d / vol_d / boxlen # mean rho * vr
#    Sz_d = mvz_d / vol_d / boxlen # mean rho * vz
#    extr = 240; extz = 80
#    Sshow = Sz_d[:extr,:extz]
#    Sshow = (np.sign(Sshow)*np.log(abs(Sshow))).T
#    plt.figure()
#    plt.imshow(Sshow,extent=(0,extr*dx*len_AU,0,extz*dx*len_AU),origin='lower',cmap='seismic')
#
    Sr_d = mvr_l / vol_l / boxlen
    Sz_d = mvz_l / vol_l / boxlen
    Sr_d[vol_l==0.]=0.
    Sz_d[vol_l==0.]=0.
#    extr = 240; extz = 80
#    Sshow = Sz_d[:extr,:extz]
#    Sshow = (np.sign(Sshow)*np.log(abs(Sshow))).T
#    plt.figure()
#    plt.imshow(Sshow,extent=(0,extr*dx*len_AU,0,extz*dx*len_AU),origin='lower',cmap='seismic')

#    H_slope = np.zeros_like(h_disk)
#    H_slope[1:-1] = (h_disk[2:]-h_disk[:-2])/(2*dx)
#    H_slope[0] = (h_disk[1]-h_disk[0])/dx; H_slope[-1] = (h_disk[-1]-h_disk[-2])/dx
    H_slope = h_disk * exp_h / rc
    H_slope2 = H_slope.copy() * fac2
#    f = plt.figure()
#    plt.plot(rc[:ind_pse],hp_all2[:ind_pse])
#    plt.plot(rc[:ind_pse],h_disk[:ind_pse])
#    plt.show()

    if f_Mdot_layer:
        f_Mdot_layer = plt.figure(figsize=(12,3)) ## S measured at several surfaces
        f_Mdot_floored = plt.figure(figsize=(12,8))
        SSz = sel_height(Sz_d,cen=True,cum=False)
        SSr = sel_height(Sr_d,cen=True,cum=False)
        ifloor = [4,5,6,7,8,10,15]
        
        Source = -SSz[0] + SSr[0] * H_slope
        Source2 = -SSz[1] + SSr[1] * H_slope2
        plt.figure(f_Mdot_layer.number)
        plt.plot(rc*len_AU, Source * kg_m2_s, label=r'$H$')
        plt.plot(rc*len_AU, Source2 * kg_m2_s, label=r'$%i H$'%fac2)
#        plt.plot(rc*len_AU, -Sz_d[:,nfloor]* kg_m2_s, label='%i.5dx'%nfloor) 
#        plt.figure(f_Mdot_floored.number)
#        plt.plot(rc*len_AU, np.cumsum(Source*surf)*Ms_yr, label=r'$H$')
#        plt.plot(rc*len_AU, np.cumsum(Source2*surf)*Ms_yr, label=r'$%i H$'%fac2)
        for fl in ifloor:
            indf = np.where(h_disk*fac2> (fl+0.5)*dx)[0][0]
            Sfl = np.append(-Sz_d[:indf,fl],Source2[indf:])
            plt.figure(f_Mdot_layer.number)
            plt.plot(rc*len_AU,Sfl* kg_m2_s, label='%.1f dx'%(fl+0.5))
            Mdot_fl = np.cumsum(Sfl*surf)
            plt.figure(f_Mdot_floored.number)
            plt.plot(rc*len_AU,Mdot_fl*Ms_yr, label='%.1f dx'%(fl+0.5))
#        plt.plot(rc*len_AU, -Sz_d[:,5]* kg_m2_s, label='5.5dx')
#        plt.plot(rc*len_AU, -Sz_d[:,7]* kg_m2_s, label='7.5dx')
#        plt.plot(rc*len_AU, -Sz_d[:,10]* kg_m2_s, label='10.5dx')
#        plt.plot(rc*len_AU, -Sz_d[:,20]* kg_m2_s, label='20.5dx')
#        plt.plot(rc*len_AU, -Sz_d[:,30]* kg_m2_s, label='30.5dx')
        plt.figure(f_Mdot_layer.number)
        plt.xscale('log');plt.yscale('symlog',linthreshy=1.e-10)
        plt.legend()
        plt.xlabel('r(AU)'); plt.ylabel('flux(kg/m$^2$/s)')
        plt.savefig(path_out+'/fluxes_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
        plt.figure(f_Mdot_floored.number)
        plt.xscale('log');plt.yscale('log')
        plt.legend()
        plt.xlabel('r(AU)'); plt.ylabel('Mdot(Ms/yr)')
        plt.savefig(path_out+'/mdots_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    ### flooring S ###
#    floor=False
    floor=nfloor
    SSr = sel_height(Sr_d,cen=True,cum=False,flo=floor)
    SSz = sel_height(Sz_d,cen=True,cum=False,flo=floor)
    H_slope[h_disk<(nfloor+0.5)*dx] = 0.
    H_slope2[h_disk*fac2<(nfloor+0.5)*dx] = 0.

    Source2 = -SSz[1] + SSr[1] * H_slope2
    Source = -SSz[0] + SSr[0] * H_slope ## rho v_perp / cos flux perpendicular to disk surface, corrected to projected disk surface
    Mdot_surf2 = Source2 * surf#[:,np.newaxis]
    Mdot_surf = Source * surf#[:,np.newaxis] # bin-wise (r) mass accretion rate onto surface
    Mdot_r_d = -Sr_d * surf[:,np.newaxis] # bin-wise (z) radial mass accretion 
    Mdot_r_c = cumz(Mdot_r_d) ## radial mass accretion rate (accumulated to z) at rs
    Mdtr = sel_height(Mdot_r_c,cum=True,flo=floor)

#    if f_Mdot_layer:
#       plt.plot(rc*len_AU, Source * kg_m2_s, ':',label=r'$H$ floored')
#       plt.plot(rc*len_AU, Source2 * kg_m2_s, ':',label=r'$%i H$ floored'%fac2)
#       plt.legend()
#       plt.savefig(path_out+'/fluxes_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    if f_flux:
        Mdotz = -np.cumsum(Sz_d * surf[:,np.newaxis], axis=0)
        Mdotr = -np.cumsum(Sr_d,axis=1)* surf[:,np.newaxis] 

        ylim_Mdotr = [-1.e2,1.e2]
        f_Mdotr = plt.figure(figsize=(7,3))
        plt.loglog(rc*len_AU, -mvrsph/Vsph*4*np.pi*rc**2*Ms_Myr,label=r'$\dot{M}_r$')
        for iz in range(20,80,8):
            plt.loglog(rc*len_AU, Mdotz[:,iz]*Ms_Myr,label=zc[iz]*len_AU)
        plot_reso(f_Mdotr,ylim_Mdotr)
        plot_Rk_Rp(f_Mdotr,ylim_Mdotr)
        plt.yscale('symlog',linthreshy=1.e-2)
        plt.xlabel('r(AU)');plt.ylabel(r'$\dot{M}(M_\odot/Myr)$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_Mdotr)
        plt.legend()
        plt.savefig(path_out+'/cylinder_Mdot_z_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        f_Mdotr2 = plt.figure(figsize=(7,3))
        plt.loglog(rc*len_AU, -mvrsph/Vsph*4*np.pi*rc**2*Ms_Myr,label=r'$\dot{M}_r$')
        for iz in range(20,80,8):
            plt.loglog(rc*len_AU, Mdotr[:,iz] *Ms_Myr,label=zc[iz]*len_AU)
        plot_reso(f_Mdotr2,ylim_Mdotr)
        plot_Rk_Rp(f_Mdotr2,ylim_Mdotr)
        plt.yscale('symlog',linthreshy=1.e-2)
        plt.xlabel('r(AU)');plt.ylabel(r'$\dot{M}(M_\odot/Myr)$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_Mdotr)
        plt.legend()
        plt.savefig(path_out+'/cylinder_Mdot_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')        

        f_Mdotr3 = plt.figure(figsize=(7,3))
        plt.loglog(rc*len_AU, -mvrsph/Vsph*4*np.pi*rc**2*Ms_Myr,label=r'$\dot{M}_r(V)$')
        plt.loglog(rc*len_AU, -rhoSvrsph/Ssph*4*np.pi*rc**2*Ms_Myr,label=r'$\dot{M}_r(S)$')
        for iz in range(20,80,8):
            plt.loglog(rc*len_AU, (Mdotz[:,iz]+Mdotr[:,iz])*Ms_Myr,label=zc[iz]*len_AU)
        plot_reso(f_Mdotr3,ylim_Mdotr)
        plot_Rk_Rp(f_Mdotr3,ylim_Mdotr)
        plt.yscale('symlog',linthreshy=1.e-2)
        plt.xlabel('r(AU)');plt.ylabel(r'$\dot{M}(M_\odot/Myr)$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_Mdotr)
        plt.legend()
        plt.savefig(path_out+'/cylinder_Mdot_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    Mdot_cum = np.cumsum(Mdot_surf); Mdot_cum = Mdot_cum[ind_kep]-Mdot_cum
    Mdot_cum2 = np.cumsum(Mdot_surf2); Mdot_cum2 = Mdot_cum2[ind_kep]-Mdot_cum2
    ifloor = [4] #,6,7,8,10,15]

    for fl in ifloor:
        indf = np.where(h_disk> (fl+0.5)*dx)[0][0]
        indf2 = np.where(h_disk*fac2> (fl+0.5)*dx)[0][0]
        Sfl = np.append(-Sz_d[:indf,fl],Source[indf:])
        Sfl2 = np.append(-Sz_d[:indf2,fl],Source2[indf2:])
        Mdot_fl = np.cumsum(Sfl*surf)
        Mdot_fl2 = np.cumsum(Sfl2*surf)
        Mdotr_fl = np.append(Mdot_r_c[:indf,fl],Mdtr[0][indf:])
        Mdotr_fl2 = np.append(Mdot_r_c[:indf2,fl],Mdtr[1][indf2:])

    if f_S: ## plot source function and Mdot
#    Source2 = -SSz[1] + SSr[1] * H_slope2
#    Source = -SSz[0] + SSr[0] * H_slope ## rho v_perp / cos flux perpendicular to disk surface, corrected to projected disk surface
#    Mdot_surf2 = Source2 * surf#[:,np.newaxis]
#    Mdot_surf = Source * surf#[:,np.newaxis] # bin-wise (r) mass accretion rate onto surface
#    Mdot_r_d = -Sr_d * surf[:,np.newaxis] # bin-wise (z) radial mass accretion 
#    Mdot_r_c = cumz(Mdot_r_d) ## radial mass accretion rate (accumulated to z) at rs
#    Mdtr = sel_height(Mdot_r_c,cum=True,flo=floor)
#        SSz = sel_height(Sz_d,cen=True,cum=False)
#        SSr = sel_height(Sr_d,cen=True,cum=False)
        ifloor = [4] #,6,7,8,10,15]

        for fl in ifloor:
            indf = np.where(h_disk> (fl+0.5)*dx)[0][0]
            indf2 = np.where(h_disk*fac2> (fl+0.5)*dx)[0][0]
            Sfl = np.append(-Sz_d[:indf,fl],Source[indf:])
            Sfl2 = np.append(-Sz_d[:indf2,fl],Source2[indf2:])
            Mdot_fl = np.cumsum(Sfl*surf)
            Mdot_fl2 = np.cumsum(Sfl2*surf)
            Mdotr_fl = np.append(Mdot_r_c[:indf,fl],Mdtr[0][indf:])
            Mdotr_fl2 = np.append(Mdot_r_c[:indf2,fl],Mdtr[1][indf2:])

        ylim_S = [1.e-12,1.e-4]
        f_S = plt.figure(figsize=(7,3))
        plt.loglog(rc[:ind_pse]*len_AU, Source[:ind_pse]*kg_m2_s,'C1',label=r'$H$')
        plt.plot(rc[:ind_pse]*len_AU, Source2[:ind_pse]*kg_m2_s,'C2',label=r'$%i H$'%fac2)
        plt.plot(rc[:indf]*len_AU, -Sz_d[:indf,fl]*kg_m2_s,'gray',zorder=0,label=r'$4dx$')
        plot_reso(f_S,ylim_S)
        plot_Rk_Rp(f_S,ylim_S)
        plt.xlabel('r(AU)');plt.ylabel(r'$S(kg/m^2/s)$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_S)
        plt.legend()
#        plt.savefig(path_out+'/disk_S_p_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
        plt.savefig(path_out+'/disk_S_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

#    if f_Mdot:    
        ylim_Mdot = [-1.e-4,1.e-4]
        f_Mdot = plt.figure(figsize=(7,3))
#        plt.loglog(rc*len_AU, Mdot_cum*Ms_yr,'C1',label=r'$\dot{M}_{s}(H)$')
#        plt.plot(rc*len_AU, Mdot_cum2*Ms_yr,'C2',label=r'$\dot{M}_{s}(%i H)$'%fac2)
#        plt.loglog(rc[:ind_kep]*len_AU, Mdtr[0][:ind_kep]*Ms_yr,'C1--',label=r'$\dot{M}_{r}(H)$')
#        plt.plot(rc[:ind_kep]*len_AU, Mdtr[1][:ind_kep]*Ms_yr,'C2--',label=r'$\dot{M}_{r}(%i H)$'%fac2)

        plt.loglog(rc[:ind_kep]*len_AU, Mdot_fl[:ind_kep]*Ms_yr,'C1',label=r'$\dot{M}_{s}(H)$')
        plt.plot(rc[:ind_kep]*len_AU, Mdot_fl2[:ind_kep]*Ms_yr,'C2',label=r'$\dot{M}_{s}(%i H)$'%fac2)
        plt.loglog(rc[:ind_kep]*len_AU, Mdotr_fl[:ind_kep]*Ms_yr,'C1--',label=r'$\dot{M}_{r}(H)$')
        plt.plot(rc[:ind_kep]*len_AU, Mdotr_fl2[:ind_kep]*Ms_yr,'C2--',label=r'$\dot{M}_{r}(%i H)$'%fac2)

#        plt.loglog(rc*len_AU, np.cumsum(Mdot_surf)*Ms_Myr,'C1',label=r'$\dot{M}_{surface}(3 au,<r)$')
#        plt.plot(rc*len_AU, np.cumsum(Mdot_surf2)*Ms_Myr,'C2',label=r'$\dot{M}_{surface}(4 au,<r)$')
#        plt.loglog(rc[:ind_pse]*len_AU, Mdtr[:ind_pse]*Ms_Myr,'C1--',label=r'$\dot{M}_{radial}(3 au,r)$')
#        plt.plot(rc[:ind_pse]*len_AU, Mdtr2[:ind_pse]*Ms_Myr,'C2--',label=r'$\dot{M}_{radial}(4 au,r)$')
        plot_reso(f_Mdot,ylim_Mdot)
        plot_Rk_Rp(f_Mdot,ylim_Mdot)
        plt.yscale('symlog',linthreshy=1.e-7)
        plt.xlabel('r(AU)');plt.ylabel(r'$\dot{M}(M_\odot/yr)$')
        plt.xlim([xmin,xmax]); plt.ylim(ylim_Mdot)
        plt.legend(fontsize=13)
#        plt.savefig(path_out+'/disk_Mdot_p_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
        plt.savefig(path_out+'/disk_Mdot_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    if f_S and overplot_S: 
        S_p = lambda r: 1./(4.* r * np.sqrt(1-r))
        S_k = lambda r: 1./(8.*r**1.5 * np.sqrt(1-np.sqrt(r)))
        S_scale = Mdot_fl[ind_kep] / ( np.pi * rr[ind_kep]**2 ) *kg_m2_s
        S_scale2 = Mdot_fl2[ind_kep] / ( np.pi * rr[ind_kep]**2 ) *kg_m2_s
        r_norm = rc[:ind_kep+1]/rr[ind_kep]
        plt.figure(f_S.number)
        plt.plot(rc[:ind_kep+1]*len_AU,S_p(r_norm)*S_scale,'C1',alpha=0.5)
        plt.plot(rc[:ind_kep+1]*len_AU,S_k(r_norm)*S_scale,'C1',alpha=0.5)
        plt.plot(rc[:ind_kep+1]*len_AU,S_p(r_norm)*S_scale2,'C2',alpha=0.5)
        plt.plot(rc[:ind_kep+1]*len_AU,S_k(r_norm)*S_scale2,'C2',alpha=0.5)
        plt.savefig(path_out+'/disk_S_c_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    print 'Mdot_disk', np.sum(Mdot_surf[:ind_pse])*Ms_Myr, np.sum(Mdot_surf2[:ind_pse])*Ms_Myr, np.sum(Mdot_surf2[:ind_pse])*Ms_Myr+Mdtr[1][ind_pse]*Ms_Myr
#    plt.show()

### binning function to smooth out fluctuations ###
    def binning(x,n):
        xp = x[::n].copy()
        np = len(x)/n
        nc = len(x)%n
        for i in range(1,nc): xp += x[i::n]
        for i in range(nc,n): xp[:np] += x[i::n]
        return xp
    ratio_bin = lambda W,P,n: binning(W,n)/binning(P,n)

### different alpha terms ###
#    Trv_l = sel_height(mdvrdvphi_l,cum=False,cen=True)
    Trv_c = sel_height(mdvrdvphi_c)
    Tzv_l = sel_height(mdvzdvphi_l,cum=False,cen=True)
#    TrB_l = sel_height(VBrBphi_l,cum=False,cen=True)
    TrB_c = sel_height(VBrBphi_c)
    TzB_l = sel_height(VBzBphi_l,cum=False,cen=True)
#    Trg_l = sel_height(Vgrgphi_l,cum=False,cen=True)
    Trg_c = sel_height(Vgrgphi_c)
    Tzg_l = sel_height(Vgzgphi_l,cum=False,cen=True)
    VV_l = sel_height(vol_l,cum=False,cen=True)
    VV_c = sel_height(vol_c)

    alphaT_h = Trv_c[0]/PP_c[0]
    alphaG_h = Trg_c[0]/PP_c[0]
    freq_alphaT = np.fft.fft(alphaT_h)
    freq_alphaG = np.fft.fft(alphaG_h)
    freq = np.fft.fftfreq(alphaT_h.shape[-1])
    nf = len(freq)/2


### 2D binning function to smooth out fluctuations ###
    def binning2D(x,n,cum=False):
        xdim = x.shape
        np0 = xdim[0]/n; nc0 = xdim[0]%n
        np1 = xdim[1]/n; nc1 = xdim[1]%n
        xp0 = x[:np0*n:n,:].copy()
        for i in range(1,nc0): xp0 += x[i::n,:]
        for i in range(nc0,n): xp0[:np0,:] += x[i::n,:]
        if cum:
          xp1 = xp0[:,n/2:np1*n:n].copy()
        else:
          xp1 = xp0[:,:np1*n:n].copy()
          for i in range(1,nc1): xp1 += xp0[:,i::n]
          for i in range(nc1,n): xp1[:,:np1] += xp0[:,i::n]
        return xp1

### 2D plots of various alpha terms ###
    if f_alpha_imshow:
        
#        r_img = int(Rk*len_AU/10)*10; z_img = int(max(r_img*0.4/5,1))*5
        r_img = int(np.ceil(Rk*len_AU/10))*10; z_img = r_img/2
        colornorm = lambda lim1,lim2: SymLogNorm(10**lim1, vmin = -10**lim2, vmax = 10**lim2)
        norm_a = colornorm(-3,-1)
        rlast=int(r_img/len_AU/dx); zlast = int(z_img/len_AU/dx)
        print r_img, z_img, len(rr),len(zz)
        print r_img, z_img, rlast, zlast, rr[rlast]*len_AU, zz[zlast]*len_AU
        imgext = np.asarray([rr[0],rr[rlast],zz[0],zz[zlast]])*len_AU
        #rgrid, zgrid = np.meshgrid(rr[:rlast]*len_AU,zz[:zlast]*len_AU)
        ind_end_h = ind_kep

        if ndx>1:
            rlast = (rlast / ndx)*ndx; zlast = (zlast / ndx)*ndx
            imgext = np.asarray([rr[0],rr[rlast],zz[0],zz[zlast]])*len_AU

# cumulated tur_r
#        alphaTr = mdvrdvphi_c[:rlast,:zlast]/mcc_c[:rlast,:zlast]
        TVtur_bin = binning2D(mdvrdvphi_c[:rlast,:zlast],ndx,cum=True)
        PV_bin = binning2D(mcc_c[:rlast,:zlast],ndx,cum=True)
        alphaTr = TVtur_bin/PV_bin
        plt.figure(figsize=(7,3))
        plt.imshow((alphaTr).T,origin='lower',cmap='seismic',extent=imgext,norm=norm_a,aspect='equal')
        plt.colorbar(extend='both')
        #plt.contour(rgrid,zgrid,alphaTr.T,[0.])
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
        plt.title(r'$\alpha_{Rey,r} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/alpha_tur_r_img_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

# cumulated mag_r
#        alphaBr = VBrBphi_c[:rlast,:zlast]/mcc_c[:rlast,:zlast]
        TVmag_bin = binning2D(VBrBphi_c[:rlast,:zlast],ndx,cum=True)
        alphaBr = TVmag_bin/PV_bin
        plt.figure(figsize=(7,3))
        plt.imshow((alphaBr).T,origin='lower',cmap='seismic',extent=imgext,norm=norm_a,aspect='equal')
        plt.colorbar(extend='both')
        #plt.contour(rgrid,zgrid,alphaBr.T,[0.])
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
        plt.title(r'$\alpha_{Max,r} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/alpha_mag_r_img_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

# local tur_z
#        rhoc2 = mcc_c[:rlast,:zlast]/vol_c[:rlast,:zlast] 
#        alphaTz = mdvzdvphi_l[:rlast,:zlast]/vol_l[:rlast,:zlast]/rhoc2
        V_bin = binning2D(vol_c[:rlast,:zlast],ndx,cum=True)
        P_bin = PV_bin/V_bin
        Vl_bin = binning2D(vol_l[:rlast,:zlast],ndx)
        TVturl_bin = binning2D(mdvzdvphi_l[:rlast,:zlast],ndx)
        alphaTz = TVturl_bin / Vl_bin / P_bin
        plt.figure(figsize=(7,3))
        plt.imshow((alphaTz).T,origin='lower',cmap='seismic',extent=imgext,norm=norm_a,aspect='equal')
        plt.colorbar(extend='both')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
        plt.title(r'$\alpha_{Rey,z} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/alpha_tur_z_img_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')


#        alphaBz = VBzBphi_l[:rlast,:zlast]/vol_l[:rlast,:zlast]/rhoc2
        TVmagl_bin = binning2D(VBzBphi_l[:rlast,:zlast],ndx)
        alphaBz = TVmagl_bin / Vl_bin / P_bin
        plt.figure(figsize=(7,3))
        plt.imshow((alphaBz).T,origin='lower',cmap='seismic',extent=imgext,norm=norm_a,aspect='equal')
        plt.colorbar(extend='both')
        #plt.contour(rgrid,zgrid,alphaBz.T,[0.])
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
        plt.title(r'$\alpha_{Max,z} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/alpha_mag_z_img_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

#    plt.show()

    if f_jflux_imshow:
        r_img = int(Rk*len_AU/10*2.)*10; z_img=r_img; #z_img = int(max(r_img*0.4/5,1)*3)*5
#        r_img = 150.; z_img = 90.
#        r_img = 100.; z_img = 100.
        r_img = 30.; z_img = 30.
        colornormsym = lambda lim1,lim2: SymLogNorm(10**lim1, vmin = -10**lim2, vmax = 10**lim2)
        colornorm = lambda lim1,lim2: LogNorm(vmin = 10**lim1,vmax = 10**lim2)
        rlast=int(r_img/len_AU/dx); zlast = int(z_img/len_AU/dx)
        rlast = min(rlast,nr-ndx/2); zlast = min(zlast,nz-ndx/2)
        rlast = (rlast/ndx)*ndx; zlast = (zlast/ndx)*ndx; 
        imgext = np.asarray([rr[0],rr[rlast],zz[0],zz[zlast]])*len_AU
        rc_img = rr[:rlast:ndx]+dx*ndx/2. ; zc_img = zz[:zlast:ndx]+dx*ndx/2.
        rgrid, zgrid = np.meshgrid(rc_img*len_AU,zc_img*len_AU)
        ind_end_h = ind_kep

        V_bin = binning2D(vol_l[:rlast,:zlast],ndx)
        M_bin = binning2D(mass_l[:rlast,:zlast],ndx)
        Rho_bin = M_bin/V_bin
        MVphi_bin = binning2D(mvphi_l[:rlast,:zlast],ndx)
        MVr_bin = binning2D(mvr_l[:rlast,:zlast],ndx)
        MVz_bin = binning2D(mvz_l[:rlast,:zlast],ndx)
        Vphi_bin = MVphi_bin / M_bin
        J_bin = binning2D(mvphi_l[:rlast,:zlast]*rc[:rlast,np.newaxis],ndx)
        Jspe_bin = J_bin / M_bin
        J_bin = J_bin / V_bin
        MVrVphi_bin = binning2D(mvrvphi_l[:rlast,:zlast],ndx)
        MVzVphi_bin = binning2D(mvzvphi_l[:rlast,:zlast],ndx)
        Tvr = MVrVphi_bin / V_bin - MVr_bin * Vphi_bin / V_bin
        Tvz = MVzVphi_bin / V_bin - MVz_bin * Vphi_bin / V_bin
        Tv = np.sqrt(Tvr**2 + Tvz**2)
#        plt.figure(figsize=(7,3)) # Turbulent M transport / density
#        stream=plt.streamplot(rgrid,zgrid,Tvr.T,Tvz.T, color=Tv.T*unit_dens*vel_ms**2,norm=colornorm(-10,-7),cmap='cool')
#        plt.colorbar(stream.lines,extend='both')
#        img=plt.imshow(Rho_bin.T*unit_dens,origin='lower',cmap='gray',extent=imgext,norm=colornorm(-15,-10),aspect='equal')
#        plt.colorbar(img,extend='both')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
#        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
#        plt.title(r'$\rho \delta u_\phi \delta u_p\,(\rho)$')
#        plt.savefig(path_out+'/Tensor_v_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        RhoVr_bin = MVr_bin / V_bin
        RhoVz_bin = MVz_bin / V_bin
        RhoV_bin = np.sqrt(RhoVr_bin**2 + RhoVz_bin**2)
        plt.figure(figsize=(7,3)) # mass flux / specific angular momentum
        stream=plt.streamplot(rgrid,zgrid,RhoVr_bin.T,RhoVz_bin.T, color=RhoV_bin.T*unit_mflux,norm=colornorm(-13,-8),cmap='cool')
        plt.colorbar(stream.lines,extend='both')
        img=plt.imshow(Jspe_bin.T*unit_jspe,origin='lower',cmap='gray',extent=imgext,norm=colornorm(15,17),aspect='equal')
        plt.colorbar(img,extend='both')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
        plt.title(r'$\rho u_p\,(ru_\phi)$')
        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
        plt.savefig(path_out+'/mflux_v_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
        #return
        jf_lVr = Vphi_bin * RhoVr_bin * rc_img[:,np.newaxis]
        jf_lVz = Vphi_bin * RhoVz_bin * rc_img[:,np.newaxis]
        jf_lV = np.sqrt(jf_lVr**2 + jf_lVz**2)
#        plt.figure(figsize=(7,3)) # Laminar J transport / J
#        stream=plt.streamplot(rgrid,zgrid,jf_lVr.T,jf_lVz.T, color=jf_lV.T*unit_jflux,norm=colornorm(2,7),cmap='cool')
#        plt.colorbar(stream.lines,extend='both')
#        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgext,norm=colornorm(0,5),aspect='equal')
#        plt.colorbar(img,extend='both')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
#        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
#        plt.title(r'$r\rho u_\phi u_p\,(\rho ru_\phi)$')
#        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.savefig(path_out+'/jflux_v_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        jf_Vr = MVrVphi_bin / V_bin * rc_img[:,np.newaxis]
        jf_Vz = MVzVphi_bin / V_bin * rc_img[:,np.newaxis]
        jf_dVr = jf_Vr - jf_lVr
        jf_dVz = jf_Vz - jf_lVz
        jf_dV = np.sqrt(jf_dVr**2 + jf_dVz**2)
#        plt.figure(figsize=(7,3)) # Turbulent J transport / J
#        stream=plt.streamplot(rgrid,zgrid,jf_dVr.T,jf_dVz.T, color=jf_dV.T*unit_jflux,norm=colornorm(2,7),cmap='cool')
#        plt.colorbar(stream.lines,extend='both')
#        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgext,norm=colornorm(0,5),aspect='equal')
#        plt.colorbar(img,extend='both')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
#        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
#        plt.title(r'$r\rho \delta u_\phi\delta  u_p\,(\rho ru_\phi)$')
#        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.savefig(path_out+'/jflux_dv_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        jf_Br = binning2D(VBrBphi_l[:rlast,:zlast],ndx) * rc_img[:,np.newaxis] / V_bin
        jf_Bz = binning2D(VBzBphi_l[:rlast,:zlast],ndx) * rc_img[:,np.newaxis] / V_bin
        jf_B = np.sqrt(jf_Br**2 + jf_Bz**2)
#        plt.figure(figsize=(7,3)) # magnetic J flux / J
#        stream=plt.streamplot(rgrid,zgrid,jf_Br.T,jf_Bz.T, color=jf_B.T*unit_jflux,norm=colornorm(2,7),cmap='cool')
#        plt.colorbar(stream.lines,extend='both')
#        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgext,norm=colornorm(0,5),aspect='equal')
#        plt.colorbar(img,extend='both')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
#        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
#        plt.title(r'$ -r B_\phi B_p\,(\rho ru_\phi)$')
#        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.savefig(path_out+'/jflux_B_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        jf_r = jf_Vr + jf_Br
        jf_z = jf_Vz + jf_Bz
        jf_tot = np.sqrt(jf_r**2 + jf_z**2)
#        plt.figure(figsize=(7,3)) # total J flux / J
#        stream=plt.streamplot(rgrid,zgrid,jf_r.T,jf_z.T, color=jf_tot.T*unit_jflux,norm=colornorm(2,7),cmap='cool')
#        plt.colorbar(stream.lines,extend='both')
#        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgext,norm=colornorm(0,5),aspect='equal')
#        plt.colorbar(img,extend='both')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
#        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
#        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
#        plt.title(r'$r\rho u_\phi u_p + r\rho \delta u_\phi \delta u_p - rB_\phi B_p\,(\rho ru_\phi)$')
#        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.savefig(path_out+'/jflux_tot_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        fig = plt.figure(figsize=(8,6)) # total J flux / J
        stream=plt.streamplot(-rgrid,zgrid,-jf_lVr.T,jf_lVz.T, color=jf_lV.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        stream=plt.streamplot(-rgrid,-zgrid,-jf_dVr.T,-jf_dVz.T, color=jf_dV.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        stream=plt.streamplot(rgrid,zgrid,jf_r.T,jf_z.T, color=jf_tot.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        stream=plt.streamplot(rgrid,-zgrid,jf_Br.T,-jf_Bz.T, color=jf_B.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        plt.colorbar(stream.lines,extend='both')
        imgextpm = imgext.copy(); imgextpm[3] *= -1
        imgextmp = imgext.copy(); imgextmp[1] *= -1
        imgextmm = imgext.copy(); imgextmm[1] *= -1; imgextmm[3] *= -1
        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgextmm,norm=colornorm(0,6),aspect='equal',zorder=0)
        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgextpm,norm=colornorm(0,6),aspect='equal',zorder=0)
        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgextmp,norm=colornorm(0,6),aspect='equal',zorder=0)
        img=plt.imshow(J_bin.T*unit_J,origin='lower',cmap='gray',extent=imgext,norm=colornorm(0,6),aspect='equal',zorder=0)
        plt.colorbar(img,extend='both')
        #cbaxes = fig.add_axes([0.1, 0.1, 0.03, 0.8])  # This is the position for the colorbar
        #cb = plt.colorbar(img, cax = cbaxes,extend='both')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='y',zorder=1)
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,-hp_all2[:min(ind_end_h,rlast)]*len_AU,c='y',zorder=1)
        plt.plot(-rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='y',zorder=1)
        plt.plot(-rc[:min(ind_end_h,rlast)]*len_AU,-hp_all2[:min(ind_end_h,rlast)]*len_AU,c='y',zorder=1)
        plt.plot([-imgext[1],imgext[1]],[0,0],'k'); plt.plot([0,0],[-imgext[3],imgext[3]],'k')
        plt.xlim([-r_img,r_img]); plt.ylim([-z_img,z_img])
        plt.title('Angular momentum transport')
        plt.savefig(path_out+'/jflux4_tot_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        fig = plt.figure(figsize=(10,4)) # total J flux / J zoom (not binned)
        rzm = int(rlast*ndx*0.6); zzm = zlast*ndx/4
        rrgrid, zzgrid = np.meshgrid(rc[:rzm]*len_AU,zc[:zzm]*len_AU)
        #rrgrid = rgrid[:zzm,:rzm]; zzgrid = zgrid[:zzm,:rzm]
        #print rrgrid.shape, zzgrid.shape, jf_lVr[:rzm,:zzm].shape, jf_lVz[:rzm,:zzm].shape, jf_lV[:rzm,:zzm].shape
        rweight = rc[:rzm,np.newaxis] / (mass_l[:rzm,:zzm] * vol_l[:rzm,:zzm])
        jf_lr = mvr_l[:rzm,:zzm] * mvphi_l[:rzm,:zzm] * rweight
        jf_lz = mvz_l[:rzm,:zzm] * mvphi_l[:rzm,:zzm] * rweight
        jf_l = np.sqrt(jf_lr**2+jf_lz**2)
        jf_dr = mdvrdvphi_l[:rzm,:zzm] / vol_l[:rzm,:zzm] * rc[:rzm,np.newaxis]
        jf_dz = mdvzdvphi_l[:rzm,:zzm] / vol_l[:rzm,:zzm] * rc[:rzm,np.newaxis]
        jf_d = np.sqrt(jf_dr**2+jf_dz**2)
        jf_dbr = VBrBphi_l[:rzm,:zzm] / vol_l[:rzm,:zzm] * rc[:rzm,np.newaxis]
        jf_dbz = VBzBphi_l[:rzm,:zzm] / vol_l[:rzm,:zzm] * rc[:rzm,np.newaxis]
        jf_db = np.sqrt(jf_dbr**2 + jf_dbz**2)
        jf_zr = jf_lr+jf_dr+jf_dbr
        jf_zz = jf_lz+jf_dz+jf_dbz 
        jf_zoom = np.sqrt(jf_zr**2+jf_zz**2)
        jz = mvphi_l[:rzm,:zzm] * rc[:rzm,np.newaxis] / vol_l[:rzm,:zzm]
        stream=plt.streamplot(-rrgrid,zzgrid,-jf_lr.T,jf_lz.T, color=jf_l.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        stream=plt.streamplot(-rrgrid,-zzgrid,-jf_dr.T,-jf_dz.T, color=jf_d.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        stream=plt.streamplot(rrgrid,zzgrid,jf_zr.T,jf_zz.T, color=jf_zoom.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        stream=plt.streamplot(rrgrid,-zzgrid,jf_dbr.T,-jf_dbz.T, color=jf_db.T*unit_jflux,norm=colornorm(2,8),cmap='autumn_r',zorder=2)
        plt.colorbar(stream.lines,extend='both')
        imgextz = np.array([0,rr[rzm],0,zz[zzm]]) * len_AU
        imgextpm = imgextz.copy(); imgextpm[3] *= -1
        imgextmp = imgextz.copy(); imgextmp[1] *= -1
        imgextmm = imgextz.copy(); imgextmm[1] *= -1; imgextmm[3] *= -1
        img=plt.imshow(jz.T*unit_J,origin='lower',cmap='gray',extent=imgextmm,norm=colornorm(0,6),aspect='equal',zorder=0)
        img=plt.imshow(jz.T*unit_J,origin='lower',cmap='gray',extent=imgextpm,norm=colornorm(0,6),aspect='equal',zorder=0)
        img=plt.imshow(jz.T*unit_J,origin='lower',cmap='gray',extent=imgextmp,norm=colornorm(0,6),aspect='equal',zorder=0)
        img=plt.imshow(jz.T*unit_J,origin='lower',cmap='gray',extent=imgextz,norm=colornorm(0,6),aspect='equal',zorder=0)
        plt.colorbar(img,extend='both')
        #cbaxes = fig.add_axes([0.1, 0.1, 0.03, 0.8])  # This is the position for the colorbar
        #cb = plt.colorbar(img, cax = cbaxes,extend='both')
        plt.plot(rc[:min(ind_end_h,rzm)]*len_AU,hp_all2[:min(ind_end_h,rzm)]*len_AU,c='y',zorder=1)
        plt.plot(rc[:min(ind_end_h,rzm)]*len_AU,-hp_all2[:min(ind_end_h,rzm)]*len_AU,c='y',zorder=1)
        plt.plot(-rc[:min(ind_end_h,rzm)]*len_AU,hp_all2[:min(ind_end_h,rzm)]*len_AU,c='y',zorder=1)
        plt.plot(-rc[:min(ind_end_h,rzm)]*len_AU,-hp_all2[:min(ind_end_h,rzm)]*len_AU,c='y',zorder=1)
        plt.plot([-imgextz[1],imgextz[1]],[0,0],'k'); plt.plot([0,0],[-imgextz[3],imgextz[3]],'k')
        plt.title('Angular momentum transport')
        plt.savefig(path_out+'/jflux4zoom_tot_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        Br = binning2D(VBr_l[:rlast,:zlast],ndx) / V_bin
        Bz = binning2D(VBz_l[:rlast,:zlast],ndx) / V_bin
        #Br = VBr_l[:rlast,:zlast]/vol_l[:rlast,:zlast]
        #Bz = VBz_l[:rlast,:zlast]/vol_l[:rlast,:zlast]
        B = np.sqrt(Br**2 + Bz**2)
        Bphi = binning2D(VBphi_l[:rlast,:zlast],ndx) / V_bin
        #Bphi = VBphi_l[:rlast,:zlast]/vol_l[:rlast,:zlast]
        plt.figure(figsize=(7,3))
        stream=plt.streamplot(rgrid,zgrid,Br.T,Bz.T, color=B.T*unit_B,norm=colornorm(-5,-1),cmap='cool')
        plt.colorbar(stream.lines,extend='both')
        img=plt.imshow(Bphi.T*unit_B,origin='lower',cmap='PuOr',extent=imgext,norm=colornormsym(-5,-2),aspect='equal')
        plt.colorbar(img,extend='both')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([zz[0],zz[zlast]])*len_AU)
        plt.title(r'$B_p\,(B_\phi)$')
        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.savefig(path_out+'/B_lines_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    if f_pm:
        r_img = int(Rk*len_AU/10*2.)*10; z_img=r_img; #z_img = int(max(r_img*0.4/5,1)*3)*5
        r_img = 150.; z_img = 90.
#        r_img = 40.; z_img = 40.
        colornormsym = lambda lim1,lim2: SymLogNorm(10**lim1, vmin = -10**lim2, vmax = 10**lim2)
        colornorm = lambda lim1,lim2: LogNorm(vmin = 10**lim1,vmax = 10**lim2)
        rlast=int(r_img/len_AU/dx); zlast = int(z_img/len_AU/dx)
        rlast = min(rlast,nr-ndx/2); zlast = min(zlast,nz-ndx/2)
        rlast = (rlast/ndx)*ndx; zlast = (zlast/ndx)*ndx;
        zi = nz-zlast; zf = nz+zlast
        imgext = np.asarray([rr[0],rr[rlast],zpm[zi],zpm[zf]])*len_AU
        rc_img = rr[:rlast/ndx]+dx*ndx/2. ; zc_img = zz[:zlast/ndx]+dx*ndx/2. ; zc_img = np.append(-zc_img[-1::-1],zc_img)
        rgrid, zgrid = np.meshgrid(rc_img*len_AU,zc_img*len_AU)
        ind_end_h = ind_kep

        vol_pm = grid_lev_x(dV_disc, func_grid=grid_xpm, dims=[nr,2*nz],symm=[1,0])
        mass_pm = grid_lev_x(dm, func_grid=grid_xpm, dims=[nr,2*nz],symm=[1,0])
        mvr_pm = grid_lev_x(dm*vr, func_grid=grid_xpm, dims=[nr,2*nz],symm=[-1,0])
        mvphi_pm = grid_lev_x(dm*vphi, func_grid=grid_xpm, dims=[nr,2*nz],symm=[-1,0])
        mvz_pm = grid_lev_x(dm*vz, func_grid=grid_xpm, dims=[nr,2*nz],symm=[1,0])
        mvrvphi_pm = grid_lev_x(dm*vr*vphi, func_grid=grid_xpm, dims=[nr,2*nz],symm=[1,0])
        mvzvphi_pm = grid_lev_x(dm*vz*vphi, func_grid=grid_xpm, dims=[nr,2*nz],symm=[-1,0])
        VBr_pm =  grid_lev_x(dV_disc*B_r, func_grid=grid_xpm, dims=[nr,2*nz],symm=[-1,0])
        VBphi_pm =  grid_lev_x(dV_disc*B_phi, func_grid=grid_xpm, dims=[nr,2*nz],symm=[-1,0])
        VBz_pm =  grid_lev_x(dV_disc*B_z, func_grid=grid_xpm, dims=[nr,2*nz],symm=[1,0])
        VBrBphi_pm =  grid_lev_x(dV_disc*B_r*B_phi, func_grid=grid_xpm, dims=[nr,2*nz],symm=[1,0])
        VBzBphi_pm =  grid_lev_x(dV_disc*B_z*B_phi, func_grid=grid_xpm, dims=[nr,2*nz],symm=[-1,0])

        V_bin = binning2D(vol_pm[:rlast,zi:zf],ndx)
        M_bin = binning2D(mass_pm[:rlast,zi:zf],ndx)
        Rho_bin = M_bin/V_bin
        MVphi_bin = binning2D(mvphi_pm[:rlast,zi:zf],ndx)
        MVr_bin = binning2D(mvr_pm[:rlast,zi:zf],ndx)
        MVz_bin = binning2D(mvz_pm[:rlast,zi:zf],ndx)
        Vphi_bin = MVphi_bin / M_bin
        J_bin = binning2D(mvphi_pm[:rlast,zi:zf]*rc[:rlast,np.newaxis],ndx)
        Jspe_bin = J_bin / M_bin
#        J_bin = J_bin / V_bin
        MVrVphi_bin = binning2D(mvrvphi_pm[:rlast,zi:zf],ndx)
        MVzVphi_bin = binning2D(mvzvphi_pm[:rlast,zi:zf],ndx)
        Br_bin = binning2D(VBr_pm[:rlast,zi:zf],ndx) / V_bin
        Bphi_bin = binning2D(VBphi_pm[:rlast,zi:zf],ndx) / V_bin
        Bz_bin = binning2D(VBz_pm[:rlast,zi:zf],ndx) / V_bin
      
        RhoVr_bin = MVr_bin / V_bin
        RhoVz_bin = MVz_bin / V_bin
        RhoV_bin = np.sqrt(RhoVr_bin**2 + RhoVz_bin**2)
        plt.figure(figsize=(7,6)) # mass flux / specific angular momentum
        stream=plt.streamplot(rgrid,zgrid,RhoVr_bin.T,RhoVz_bin.T, color=RhoV_bin.T*unit_mflux,norm=colornorm(-13,-8),cmap='cool')
        plt.colorbar(stream.lines,extend='both')
        img=plt.imshow(Jspe_bin.T*unit_jspe,origin='lower',cmap='gray',extent=imgext,norm=colornorm(15,17),aspect='equal')
        plt.colorbar(img,extend='both')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([-zz[zlast],zz[zlast]])*len_AU)
        plt.title(r'$\rho u_p\,(ru_\phi)$')
        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.savefig(path_out+'/mflux2_v_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        B_bin = np.sqrt(Br_bin**2 + Bz_bin**2)
        plt.figure(figsize=(7,6))
        stream=plt.streamplot(rgrid,zgrid,Br_bin.T,Bz_bin.T, color=B_bin.T*unit_B,norm=colornorm(-5,-1),cmap='cool')
        plt.colorbar(stream.lines,extend='both')
        img=plt.imshow(Bphi_bin.T*unit_B,origin='lower',cmap='PuOr',extent=imgext,norm=colornormsym(-5,-2),aspect='equal')
        plt.colorbar(img,extend='both')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*len_AU,c='w')
        plt.plot(rc[:min(ind_end_h,rlast)]*len_AU,hp_all2[:min(ind_end_h,rlast)]*fac2*len_AU,c='w')
        plt.ylim(np.asarray([-zz[zlast],zz[zlast]])*len_AU)
        plt.title(r'$B_p\,(B_\phi)$')
        #plt.title(r'$\alpha_{r,tur} (t=%i kyr, M_\ast = %.2f M_\odot)$'%(time*1.e3, Msink[isam]), fontsize=18)
#        plt.savefig(path_out+'/B_lines_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        plt.show()


    if f_alpha_bins: ## plot fft of alpha
        f_fft = plt.figure(figsize=(7,3))
        print 'dx', dx*len_AU
        plt.plot(freq[:nf]/(dx*len_AU), np.absolute(freq_alphaT[:nf])/np.absolute(freq_alphaT[0]),label=r'$\alpha_{Rey,r}$')
        plt.plot(freq[:nf]/(dx*len_AU), np.absolute(freq_alphaG[:nf])/np.absolute(freq_alphaG[0]),'--',label=r'$\alpha_{Grav,r}$')
        plt.yscale('log')
        plt.legend(fontsize=16,frameon=False)
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.xlabel(r'$1/\Delta r (AU^{-1})$',fontsize=16)
        plt.ylabel(r'$\|\widetilde{\mathcal{F}(\alpha)}\|$',fontsize=16)
        plt.savefig(path_out+'/fft_alpha_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    alphaT_bin = [Trv_c[0]/PP_c[0]]
    alphaB_bin = [TrB_c[0]/PP_c[0]]
    alphaG_bin = [Trg_c[0]/PP_c[0]]
    bins = [2,4,8,16]
    for bin in bins:
        alphaT_bin.append(ratio_bin(Trv_c[0],PP_c[0],bin))
        alphaB_bin.append(ratio_bin(TrB_c[0],PP_c[0],bin))
        alphaG_bin.append(ratio_bin(Trg_c[0],PP_c[0],bin))
    ind_rmax=nr
    if f_alpha_bins: ## plot alpha at several binning size 
        ylim_alpha=[-1.e1,1.e1]
        f_alpha_Rey_b = plt.figure(figsize=(7,3))
        plt.plot(rc[:ind_rmax]*len_AU,alphaT_bin[0],label='dx')
        for bin, alp in zip(bins,alphaT_bin[1:]):
            plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alp,label=str(bin)+'dx')
        plot_reso(f_alpha_Rey_b,ylim_alpha)
        plot_Rk_Rp(f_alpha_Rey_b,ylim_alpha)
        plt.xscale('log');plt.yscale('symlog',linthreshy=1.e-4)
        plt.xlim([xmin,xmax]); plt.ylim(ylim_alpha)
        plt.legend()
        plt.xlabel(r'$r(AU)$')#,fontsize=16)
        plt.ylabel(r'$\alpha_{Rey,r}$')#,fontsize=16)
        plt.savefig(path_out+'/alpha_Rey_bins_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        f_alpha_Max_b = plt.figure(figsize=(7,3))
        plt.plot(rc[:ind_rmax]*len_AU,alphaB_bin[0],label='dx')
        for bin, alp in zip(bins,alphaB_bin[1:]):
            plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alp,label=str(bin)+'dx')
        plot_reso(f_alpha_Max_b,ylim_alpha)
        plot_Rk_Rp(f_alpha_Max_b,ylim_alpha)
        plt.xscale('log');plt.yscale('symlog',linthreshy=1.e-4)
        plt.xlim([xmin,xmax]); plt.ylim(ylim_alpha)
        plt.legend()
        plt.xlabel(r'$r(AU)$')
        plt.ylabel(r'$\alpha_{Max,r}$')
        plt.savefig(path_out+'/alpha_Max_bins_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        f_alpha_Gra_b = plt.figure(figsize=(7,3))
        plt.plot(rc[:ind_rmax]*len_AU,alphaG_bin[0],label='dx')
        for bin, alp in zip(bins,alphaG_bin[1:]):
            plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alp,label=str(bin)+'dx')
        plot_reso(f_alpha_Gra_b,ylim_alpha)
        plot_Rk_Rp(f_alpha_Gra_b,ylim_alpha)
        plt.xscale('log');plt.yscale('symlog',linthreshy=1.e-4)
        plt.xlim([xmin,xmax]); plt.ylim(ylim_alpha)
        plt.legend()
        plt.xlabel(r'$r(AU)$')
        plt.ylabel(r'$\alpha_{Grav,r}$')
        plt.savefig(path_out+'/alpha_Grav_bins_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#        print 'alpha_Rey', np.sum(WW)/np.sum(PP)
#        print 'alpha_Max', np.sum(BB)/np.sum(PP)
#        print 'alpha_Grav', np.sum(GG)/np.sum(PP)
#        plt.show()

    if f_alpha: ## plot different alpha terms
        ylim_alpha=[-1.e1,1.e1]
        bin = 4
#        alphaN = 78 *  (mass_sink *  boxlen**2 ) / (CC_c[0] * rc) * (rc/dx)**(-2.85)
#        alphaN2 = 78 *  (mass_sink *  boxlen**2 ) / (CC_c[1] * rc) * (rc/dx)**(-2.85)
#        rBondi = (mass_sink *  boxlen**3 ) / (CC_c[1]) 
#        print 'r_Bondi', rBondi*len_AU
        alphaTr = ratio_bin(Trv_c[0],PP_c[0],bin)
        alphaTr2 = ratio_bin(Trv_c[1],PP_c[1],bin)
        alphaBr = ratio_bin(TrB_c[0],PP_c[0],bin)
        alphaBr2 = ratio_bin(TrB_c[1],PP_c[1],bin)
#        alphaGr = ratio_bin(Trg_c[0],PP_c[0],bin)
#        alphaGr2 = ratio_bin(Trg_c[1],PP_c[1],bin)
        f_alpha_r = plt.figure(figsize=(7,3))
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaTr,label=r'$\alpha_{Rey,r}(H)$',zorder=2)
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaTr2,'C0--',label=r'$\alpha_{Rey,r}(%i H)$'%fac2,zorder=2)
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaBr,label=r'$\alpha_{Max,r}(H)$',zorder=1)
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaBr2,'C1--',label=r'$\alpha_{Max,r}(%i H)$'%fac2,zorder=1)
#        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaGr,label=r'$\alpha_{Grav}(H)$',zorder=0,alpha=0.6)
#        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaGr2,'C2--',label=r'$\alpha_{Grav}(%i H)$'%fac2,zorder=0,alpha=0.6)
        plot_reso(f_alpha_r,ylim_alpha)
        plot_Rk_Rp(f_alpha_r,ylim_alpha)
        plt.xscale('log');plt.yscale('symlog',linthreshy=1.e-4)
        plt.xlim([xmin,xmax]); plt.ylim(ylim_alpha)
        plt.legend(handletextpad=0.1)
        plt.xlabel(r'$r(AU)$'); plt.ylabel(r'$\alpha_{r,%i dx}$'%(bin))
#        plt.show()
        plt.savefig(path_out+'/alpha_r_'+str(bin)+'dx_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

        rcc0 = ratio_bin(PP_c[0],VV_c[0],bin)
        rcc1 = ratio_bin(PP_c[1],VV_c[1],bin)
        alphaTz = ratio_bin(Tzv_l[0],VV_l[0],bin)/rcc0
        alphaTz2 = ratio_bin(Tzv_l[1],VV_l[1],bin)/rcc1
        alphaBz = ratio_bin(TzB_l[0],VV_l[0],bin)/rcc0
        alphaBz2 = ratio_bin(TzB_l[1],VV_l[1],bin)/rcc1
#        alphaGz = ratio_bin(Tzg_l[0],VV_l[0],bin)/rcc0
#        alphaGz2 = ratio_bin(Tzg_l[1],VV_l[1],bin)/rcc1
        f_alpha_z = plt.figure(figsize=(7,3))
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaTz,label=r'$\alpha_{Rey,z}(H)$',zorder=2)
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaTz2,'C0--',label=r'$\alpha_{Rey,z}(%i H)$'%fac2,zorder=2)
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaBz,label=r'$\alpha_{Max,z}(H)$',zorder=1)
        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaBz2,'C1--',label=r'$\alpha_{Max,z}(%i H)$'%fac2,zorder=1)
#        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaGz,label=r'$\alpha_{Grav}(H)$',zorder=0,alpha=0.6)
#        plt.plot((rc[:ind_rmax:bin]+dx*(bin/2))*len_AU,alphaGz2,'C2--',label=r'$\alpha_{Grav}(%i H)$'%fac2,zorder=0,alpha=0.6)
        plot_reso(f_alpha_z,ylim_alpha)
        plot_Rk_Rp(f_alpha_z,ylim_alpha)
        plt.xscale('log');plt.yscale('symlog',linthreshy=1.e-4)
        plt.xlim([xmin,xmax]); plt.ylim(ylim_alpha)
        plt.legend(handletextpad=0.1)
        plt.xlabel(r'$r(AU)$'); plt.ylabel(r'$\alpha_{z,%i dx}$'%(bin))
        plt.savefig(path_out+'/alpha_z_'+str(bin)+'dx_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#    plt.show()

### Mass accretion rate from alpha ###
    nbin = 4
#    alphar = ratio_bin(Trv_c[0]+TrB_c[0]+Trg_c[0],PP_c[0],nbin)
#    alphar2 = ratio_bin(Trv_c[1]+TrB_c[1]+Trg_c[1],PP_c[1],nbin)
    alphar = ratio_bin(Trv_c[0]+TrB_c[0],PP_c[0],nbin)
    alphar2 = ratio_bin(Trv_c[1]+TrB_c[1],PP_c[1],nbin)
    rcc0 = ratio_bin(PP_c[0],VV_c[0],nbin)
    rcc1 = ratio_bin(PP_c[1],VV_c[1],nbin)
#    alphaz = ratio_bin(Tzv_l[0]+TzB_l[0]+Tzg_l[0],VV_l[0],nbin)/rcc0
#    alphaz2 = ratio_bin(Tzv_l[1]+TzB_l[1]+Tzg_l[1],VV_l[1],nbin)/rcc1
    alphaz = ratio_bin(Tzv_l[0]+TzB_l[0],VV_l[0],nbin)/rcc0
    alphaz2 = ratio_bin(Tzv_l[1]+TzB_l[1],VV_l[1],nbin)/rcc1
    rbin = rc[::nbin]+dx*nbin/2
    rbinc = 0.5*(rbin[1:]+rbin[:-1])
    Vp_binned = ratio_bin(MVp_c[0],MM_c[0],nbin)
    Vp_binned2 = ratio_bin(MVp_c[1],MM_c[1],nbin)
    SigC2_binned = ratio_bin(PP_c[0],surf,nbin)
    SigC2_binned2 = ratio_bin(PP_c[1],surf,nbin)
#    hp_binned = ratio_bin(h_disk,np.ones_like(h_disk),nbin)
    hp_binned = np.exp(fac_h)*rbin**exp_h

    R2Omep = np.diff(rbin * Vp_binned / boxlen)
    R2SigC2Alprp = np.diff(rbin**2 * SigC2_binned * alphar / boxlen**2)
    R2SigC2Alpz = rbin**2 * SigC2_binned * alphaz / boxlen**2
    diff_z = (R2SigC2Alpz[1:] + R2SigC2Alpz[:-1]) / (hp_binned[1:]+hp_binned[:-1]) *dx*nbin 

    R2Ome2p = np.diff(rbin * Vp_binned2 / boxlen)
    R2SigC2Alpr2p = np.diff(rbin**2 * SigC2_binned2 * alphar2 / boxlen**2)
    R2SigC2Alpz2 = rbin**2 * SigC2_binned2 * alphaz2 / boxlen**2
    diff_z2 = (R2SigC2Alpz2[1:] + R2SigC2Alpz2[:-1]) / (hp_binned[1:]+hp_binned[:-1]) *dx*nbin


    Mdot_alphar = 2*np.pi * R2SigC2Alprp/R2Omep
    Mdot_alphaz = 2*np.pi * diff_z/R2Omep
    Mdot_alphar2 = 2*np.pi * R2SigC2Alpr2p/R2Ome2p
    Mdot_alphaz2 = 2*np.pi * diff_z2/R2Ome2p
 
    Mdot_alpha = Mdot_alphar + Mdot_alphaz
    Mdot_alpha2 = Mdot_alphar2 + Mdot_alphaz2
    if f_Mdot and f_Mdot_alpha:
        plt.figure(f_Mdot.number)
        #plt.plot(rbinc[:ind_pse/nbin]*len_AU,Mdot_alphar[:ind_pse/nbin]*Ms_Myr,'C1',label=r'$\dot{M}_\alpha(H)$',alpha=0.5,zorder=1)
        #plt.plot(rbinc[:ind_pse/nbin]*len_AU,Mdot_alphar2[:ind_pse/nbin]*Ms_Myr,'C2',label=r'$\dot{M}_\alpha(%i H)$'%fac2,alpha=0.5,zorder=1)
        #plt.plot(rbinc[:ind_pse/nbin]*len_AU,Mdot_alphaz[:ind_pse/nbin]*Ms_Myr,'C1--',label=r'$\dot{M}_\alpha(H)$',alpha=0.5,zorder=1)
        #plt.plot(rbinc[:ind_pse/nbin]*len_AU,Mdot_alphaz2[:ind_pse/nbin]*Ms_Myr,'C2--',label=r'$\dot{M}_\alpha(%i H)$'%fac2,alpha=0.5,zorder=1)
        plt.plot(rbinc[:ind_kep/nbin]*len_AU,Mdot_alpha[:ind_kep/nbin]*Ms_yr,'C1:',label=r'$\dot{M}_\alpha(H)$',alpha=0.7,zorder=1)
        plt.plot(rbinc[:ind_kep/nbin]*len_AU,Mdot_alpha2[:ind_kep/nbin]*Ms_yr,'C2:',label=r'$\dot{M}_\alpha(%i H)$'%fac2,alpha=0.7,zorder=1)
        plt.legend(frameon=False,labelspacing=0.05)
        plt.savefig(path_out+'/disk_Mdot_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#    plt.show() 


    plt.close('all')
    filename = path_out+'/disk_data_f_'+str(ioutput).zfill(5)+'.npz'
    np.savez( filename, time=time_list[isam], Rk=Rk*len_AU, Rp=Rp*len_AU, Ms=Msink[isam], MRk=mass_Rk*units["mass_sol"], MRp=mass_Rp*units["mass_sol"], Rsnow=np.asarray(Rsnow)*len_AU, Mdotsurf=[np.sum(Mdot_surf2[:ind_kep])*Ms_Myr, np.sum(Mdot_surf2[:ind_pse])*Ms_Myr,Mdot_cum2[4/ndx0]*Ms_Myr],fac_h=fac_h,exp_h=exp_h,dx=dx,len_AU=len_AU,Mdotrad=[Mdtr[1][ind_kep]*Ms_Myr,Mdtr[1][ind_pse]*Ms_Myr], Mdotsurffl=[Mdot_fl[ind_kep]*Ms_Myr,Mdot_fl2[ind_kep]*Ms_Myr] )
    return [time_list[isam], Rk*len_AU, Rp*len_AU, Msink[isam], mass_Rk*units["mass_sol"], mass_Rp*units["mass_sol"], np.asarray(Rsnow)*len_AU , np.sum(Mdot_surf2[:ind_kep])*Ms_Myr, np.sum(Mdot_surf2[:ind_pse])*Ms_Myr]

    if False: ## plot temperature
        f_Th = plt.figure(figsize=(7,3))
        plt.plot(rc[:ind_rmax]*len_AU, TT/MM, label='$T(h)$')
        plt.plot(rc[:ind_rmax]*len_AU, TT2/MM2, label='$T(2h)$')
        plt.plot(rc[:ind_rmax]*len_AU, TTc/MMc, label='$T(z_{crit})$')
        plt.loglog(rc[:ind_rmax]*len_AU, TTd/MMd, label="$T'(z_{crit})$")
        plt.loglog(rc[:ind_rmax]*len_AU, T_c[:ind_rmax,0], label="$T(0)$")
        plt.legend(frameon=True)
        plt.xlabel(r'$r(AU)$')
        plt.ylabel(r'$T (K)$')
        plt.savefig(path_out+'/T_h_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
        
        f_Tz = plt.figure(figsize=(7,3))
        for ir in range(0,ind_rmax,8):
            plt.loglog(zc*len_AU,T_d[ir,:], label='%.1f'%(rc[ir]*len_AU))
        plt.legend(frameon=True)
        plt.xlabel(r'$z(AU)$')
        plt.ylabel(r'$T (K)$')
        plt.savefig(path_out+'/T_rz_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#    return
#    plt.show()


    #symlog = lambda x: np.sign(x) * np.log10(abs(x))
    symlog = lambda x:x
    lim1 = -3; lim2 = 2
    colornorm = SymLogNorm(10**lim1, vmin = -10**lim2, vmax = 10**lim2)
    imgext = np.asarray([rr[0],rr[-1],zz[0],zz[-1]])*len_AU

#    sigma_c = mass_c / surf[:,np.newaxis] * kg_m2
#    alpha = mdvrdvphi_c/mcc_c
#    dens = mass_d / surf[:,np.newaxis] / dx
#    hght = np.tile(zz[np.newaxis,1:],(len(rr)-1,1))

#    zmin = 0; zmax = int(7/len_AU/dx) 
#    plt.figure(figsize=(7,6))
#    plt.scatter(sigma_c[:,zmin:zmax].flatten(),(alpha[:,zmin:zmax]).flatten(),c=np.log10(dens[:,zmin:zmax]),s=5.e-4/hght[:,zmin:zmax].flatten(),cmap='YlGn',alpha=0.5)
#    plt.xscale('log')
#    plt.yscale('symlog',linthreshy=1.e-3)
#    plt.xlim([1,np.amax(sigma_c)]); #plt.ylim([1.e-3,20])
#    plt.ylim([-10,20])
#    plt.colorbar()
#    for i in range(zmin,zmax,2):
#        plt.plot(sigma_c[:,i],alpha[:,i])
#    for i in range(0,len(rr)/2,5):
#        plt.plot(sigma_c[i,zmin:zmax],alpha[i,zmin:zmax],'k')#'grey',zorder=0)
#    plt.xlabel(r'$\Sigma$')
#    plt.ylabel(r'$\alpha$')
#
#    print zmax,len(rr)
#    print hght[:,zmin:zmax]

#    zmin = 0; zmax=int(7/len_AU/dx)
#    plt.figure(figsize=(7,6))
#    plt.scatter(dens[:,zmin:zmax].flatten(),(alpha[:,zmin:zmax]).flatten(),c=np.log10(dens[:,zmin:zmax]),s=5.e-4/hght[:,zmin:zmax].flatten(),cmap='YlGn',alpha=0.5)
#    plt.xscale('log')
#    plt.yscale('symlog',linthreshy=1.e-3)
##    plt.xlim([10,2.e3]); #plt.ylim([1.e-3,20])
#    plt.xlim([1.e9,np.amax(dens)])
#    plt.ylim([-10,20])
#    plt.colorbar()
#    for i in range(zmin,zmax,2):
#        plt.plot(dens[:,i],alpha[:,i])
#    for i in range(0,len(rr)/2,5):
#        plt.plot(dens[i,zmin:zmax],alpha[i,zmin:zmax],'k')#,zorder=0)#'grey',zorder=0)
#    plt.xlabel(r'$\rho_{mean}$')
#    plt.ylabel(r'$\alpha$')
    #plt.show()
#    return

#    alpha_Max = -VBrBphi_c/mcc_c
#    alpha_Maxz = -VBzBphi_d/mcc_d
#    alpha_Poi = Vgrgphi_c / mcc_c / (4*np.pi)
#    norm = mcc_c


    plt.figure(figsize=(7,3))
    plt.imshow(symlog(mdvrdvphi_c/mcc_c).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm)
    plt.colorbar(extend='both')
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.savefig(path_out+'/alpha(<z)_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
   #plt.show()
#    return
    plt.figure(figsize=(7,3))
    plt.imshow(symlog(-VBrBphi_c/mcc_c).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm)
    plt.colorbar(extend='both')
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.savefig(path_out+'/alpha_Max(<z)_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
#    return
    plt.figure(figsize=(7,3))
    norm = mcc_d
    plt.imshow(symlog(-VBzBphi_d/mcc_d).T,origin='lower',cmap='Vega20',extent=imgext,norm=colornorm)
    plt.colorbar(extend='both')
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.savefig(path_out+'/alpha_Maxz(<z)_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    plt.figure(figsize=(7,3))
    norm = mcc_c
    plt.imshow(symlog(Vgrgphi_c/mcc_c).T,origin='lower',cmap='seismic',extent=imgext,norm=colornorm)
    plt.colorbar(extend='both')
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.savefig(path_out+'/alpha_grav(<z)_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    #plt.xscale('log');plt.yscale('log')
    plt.show()

    return
    plt.figure(figsize=(7,3))
    alpha_Max = -VBrBphi[:,1] / mcc2[:,1] #/ (4.*np.pi)
    alpha_Maxz = -VBzBphi[:,1] / mcc2[:,1] #/ (4.*np.pi)
    alpha_Poi = Vgrgphi[:,1] / mcc2[:,1] / (4.*np.pi)
    plt.figure(figsize=(7,3))
    plt.plot(rc*len_AU,alpha,label=r'$\alpha_{Rey}$')
    plt.plot(rc*len_AU,alpha_Max,label=r'$\alpha_{Max}$')
#    plt.plot(rr[:-1]*len_AU,alpha_Maxz,label=r'$\alpha_{Max,z}$')
    plt.plot(rc*len_AU,alpha_Poi,label=r'$\alpha_{Grav}$')
#    plt.plot(rr[:-1]*len_AU,alpha+alpha_Max+alpha_Poi)
#    plt.plot(rr[:-1]*len_AU,sig2r,label=r'$\sigma^2_r$')
#    plt.semilogy(rr[:-1]*len_AU,sig2phi,label=r'$\sigma^2_\phi$')
    plt.yscale('symlog',linthreshy=1.e-2)
    plt.xlim([0,35])
    plt.ylim([-1.e1,1.e1])
    plt.xlabel(r'$R$ (AU)',fontsize=20)
    plt.ylabel(r'$\alpha$',fontsize=20)
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.legend(fontsize=16)
    plt.savefig(path_out+'/alpha_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    return

    T_hist = hist_x(dm2 * T_disc)/mass2
    plt.figure(figsize=(7,3))
    plt.loglog(rc*len_AU,T_hist[:,1])
    plt.xlim([1,35])
    plt.ylim([1.e1,2.e3])
    plt.xlabel(r'$R$ (AU)',fontsize=20)
    plt.ylabel(r'$T$ (K)',fontsize=20)
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.savefig(path_out+'/T_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    plt.figure(figsize=(7,3))
    Sigma_f = mass[:,1] / surf * kg_m2
    Sigma_h = mass2[:,1] / surf * kg_m2
    plt.loglog(rc*len_AU, Sigma_f)
    plt.loglog(rc*len_AU, Sigma_h)
    plt.xlim([1,35])
    plt.ylim([1.e0,1.e5])
    plt.xlabel(r'$R$ (AU)',fontsize=20)
    plt.ylabel(r'$\Sigma$ (kg/m$^2$)',fontsize=20)
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    #plt.show()
    plt.savefig(path_out+'/Sigma_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')

    #return

    VB2 = hist_x(dV2 * np.sum(B_disc**2,axis=1))
    plt.figure(figsize=(7,3))
    plt.semilogy(rc*len_AU,2.*mcc2[:,1]/VB2[:,1])
    plt.xlim([0,35])
    plt.ylim([1.e-1,1.e3])
    plt.xlabel(r'$R$ (AU)',fontsize=20)
    plt.ylabel(r'$\beta$',fontsize=20)
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.savefig(path_out+'/beta_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
   
    mvphi = hist_x( (dm2 * vphi) ) 
    Omega = mvphi[:,1] / mass2[:,1] / rc / boxlen
    Csound = np.sqrt(mcc2[:,1] / mass2[:,1] )
    Sigma = mass2[:,1] / (2*np.pi*rc * np.diff(rr))
    Qtoom = Omega * Csound / Sigma / np.pi / boxlen
    plt.figure(figsize=(7,3))
    plt.semilogy(rc*len_AU,Qtoom)
    plt.xlim([0,35])
    plt.ylim([1.e-1,1.e3])
    plt.xlabel(r'$R$ (AU)',fontsize=20)
    plt.ylabel(r'$Q$',fontsize=20)
    plt.title('Time = %i kyr, M* = %.2f Ms'%(time*1.e3, Msink[isam]), fontsize=24)
    plt.savefig(path_out+'/Q_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
    plt.show()
#############################################################################################
# plot disk mass and radius history
def plot_MR_history(path_outs):
    tick_size = 16
    fontlabel_size = 16
    axislabel_size = 16
    params = {'backend': 'wxAgg', 'lines.markersize' : 3, 'axes.labelsize': axislabel_size, 'font.size': fontlabel_size, 'legend.fontsize': fontlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'axes.linewidth' : 2., 'legend.labelspacing': 0.2}
    plt.rcParams.update(params)
    time_lim = [0,1.5e2]

    def read_history(path_out):
        Disk_Hist = read_pickle(path_out+'/disk_mr_histoty.npz',1)[0]
        print Disk_Hist

#        time = np.asarray([ h[0][0] for h in Disk_Hist])
#        Rk = np.asarray([ h[0][1] for h in Disk_Hist])
#        Rp = np.asarray([ h[0][2] for h in Disk_Hist])
#        Msink = np.asarray([ h[0][3] for h in Disk_Hist])
#        Mdisk_k = np.asarray([ h[0][4] for h in Disk_Hist])
#        Mdisk_p = np.asarray([ h[0][5] for h in Disk_Hist])

        time = np.asarray([ h[0] for h in Disk_Hist])
        Rk = np.asarray([ h[1] for h in Disk_Hist])
        Rp = np.asarray([ h[2] for h in Disk_Hist])
        Msink = np.asarray([ h[3] for h in Disk_Hist])
        Mdisk_k = np.asarray([ h[4] for h in Disk_Hist])
        Mdisk_p = np.asarray([ h[5] for h in Disk_Hist])
        Mdot_s = np.zeros_like(Msink)
        Mdot_s[1:-1] = (Msink[2:]-Msink[:-2])/(time[2:]-time[:-2])
        Mdot_s[0] = (Msink[1]-Msink[0])/(time[1]-time[0])
        Mdot_s[-1] = (Msink[-1]-Msink[-2])/(time[-1]-time[-2])
        time_kyr = time*1.e3
        Rsnow = np.asarray([h[6] for h in Disk_Hist])
        Mdot_d = np.asarray([h[7] for h in Disk_Hist])
        return time_kyr, Rk, Rp, Msink, Mdisk_k, Mdisk_p, Mdot_s, Mdot_d, Rsnow

    def read_history_ind(path_out):
        filelist = sorted(glob.glob(path_out+'/disk_data_f_*.npz'))
        time=[]; Rk=[]; Rp=[]; Msink=[]; Mdisk_k=[]; Mdisk_p=[]; Rsnow=[]; Mdot_d=[]; Mdot_r=[]
        Mdot_fl=[]
        nfile = len(filelist)
#        for temp_file in filelist:
        for ifile in range(nfile):
            #if ifile==8 or ifile==10: continue
            temp_file = filelist[ifile]
            temp = np.load(temp_file)
            if 'Mdotrad' not in temp:
                print filelist[ifile], 'outdated'
                continue
            time.append(temp['time'])
            Rk.append(temp['Rk']) 
            Rp.append(temp['Rp'])
            Msink.append(temp['Ms'])
            Mdisk_k.append(temp['MRk'])
            Mdisk_p.append(temp['MRp'])
            Rsnow.append(temp['Rsnow'])
            Mdot_d.append(temp['Mdotsurf'])
            Mdot_r.append(temp['Mdotrad'])
            Mdot_fl.append(temp['Mdotsurffl'])
        Msink = np.asarray(Msink)
        time = np.asarray(time)
        Mdot_s = np.zeros_like(Msink)
        Mdot_s[1:-1] = (Msink[2:]-Msink[:-2])/(time[2:]-time[:-2])
        Mdot_s[0] = (Msink[1]-Msink[0])/(time[1]-time[0])
        Mdot_s[-1] = (Msink[-1]-Msink[-2])/(time[-1]-time[-2])
        Mdot_d = np.asarray(Mdot_d)
        Mdot_r = np.asarray(Mdot_r)
        Mdot_fl = np.asarray(Mdot_fl)
        time_kyr = time*1.e3
        Rsnow = np.asarray(Rsnow)
        return time_kyr, Rk, Rp, Msink, Mdisk_k, Mdisk_p, Mdot_s, Mdot_d, Mdot_r, Rsnow, Mdot_fl        

#    time_kyr, Rk, Rp, Msink, Mdisk_k, Mdisk_p, Mdot_s, Mdot_d, Rsnow = read_history(path_outs[0])
    time_kyr, Rk, Rp, Msink, Mdisk_k, Mdisk_p, Mdot_s, Mdot_d, Mdot_r, Rsnow, Mdot_fl = read_history_ind(path_outs[0])
    Mdot_sc = (Msink[1:]-Msink[:-1])/(time_kyr[1:]-time_kyr[:-1])
    time_kyrc = (time_kyr[1:]+time_kyr[:-1])*0.5
    print time_kyr
    print Rsnow.shape
    print Mdisk_k, Mdisk_p, Mdot_s,Rk,Rp
    print 'Mdot_s', Mdot_sc
    print 'Mdot_d', Mdot_d, Mdot_r
    f_M = plt.figure(figsize=(7,3))
    plt.semilogy(time_kyr,Msink,lw=3,label=r'$M_\ast$')
    plt.plot(time_kyr,Mdisk_k, c='C0',label=r'$M_{Keplerian}$')
    plt.plot(time_kyr,Mdisk_p, '--',c='C0',label=r'$M_{magnetized}$')
    plt.xlim(time_lim); # plt.ylim([1.e-3,2])
    plt.legend()
    plt.xlabel('Time(kyr)')
    plt.ylabel('Mass($M_\odot$)')
#    plt.savefig(path_outs[0]+'/history_M.pdf',bbox_inches='tight')

    f_Mdot = plt.figure(figsize=(7,3))
    plt.semilogy(time_kyrc,Mdot_sc*1.e-3,lw=3,label=r'$\dot{M}_\ast$')
    #plt.plot(time_kyr,Mdot_d[:,0]*1.e-6,'--', c='C0', label=r'$\dot{M}_{Disk,3H}$')
    #plt.plot(time_kyr,(Mdot_d[:,2]+Mdot_r[:,0])*1.e-6,'--', c='C0', label=r'$\dot{M}_{Disk,3H}$')
    plt.plot(time_kyr,Mdot_fl[:,1]*1.e-6,'--', c='C0', label=r'$\dot{M}_{Disk,3H}$')
    plt.xlim(time_lim); # plt.ylim([1.e-4,1.e1])
    plt.legend()
    plt.xlabel('Time(kyr)')
    plt.ylabel(r'$\dot{M}(M_\odot/yr)$')
#    plt.savefig(path_outs[0]+'/history_Mdot.pdf',bbox_inches='tight')
    #plt.show()
    f_R = plt.figure(figsize=(7,3))
    plt.plot(time_kyr,Rk,c='C0',label=r'$R_{Keplerian}$')
    plt.plot(time_kyr,Rp,'--', c='C0',label=r'$R_{magnetized}$')
    plt.xlim(time_lim); # plt.ylim([0,40])
    plt.legend()
    plt.xlabel('Time(kyr)')
    plt.ylabel('Disk radius (AU)')
#    plt.savefig(path_outs[0]+'/history_R.pdf',bbox_inches='tight')

    f_150 = plt.figure(figsize=(7,3))
    plt.plot(time_kyr,Rsnow[:,0],c='C0',label='150 K')
    #plt.plot(time_kyr,Rsnow[:,1],ls='--', c='C0',label='1500 K')
    plt.xlim(time_lim);  plt.ylim([0,15])
    plt.legend()
    plt.xlabel('Time(kyr)')
    plt.ylabel('Snowline (AU)')
#    plt.savefig(path_outs[0]+'/history_snowline_150.pdf',bbox_inches='tight')

    label = ['','R_40ky_l18', 'R_80ky_l18', 'R_l14_nofeed', 'R_40ky_l17', 'R_40ky_l18_lowAD']
   # plt.show()
    for irun in range(1,len(path_outs)):
        print 'irun', irun
#        time_kyr, Rk, Rp, Msink, Mdisk_k, Mdisk_p, Mdot_s, Mdot_d, Rsnow = read_history(path_outs[irun])
        time_kyr, Rk, Rp, Msink, Mdisk_k, Mdisk_p, Mdot_s, Mdot_d, Mdot_r, Rsnow, Mdot_fl = read_history_ind(path_outs[irun])
        Mdot_sc = (Msink[1:]-Msink[:-1])/(time_kyr[1:]-time_kyr[:-1])
        time_kyrc = (time_kyr[:1:]+time_kyr[:-1])*0.5

        print 'time', time_kyr
        print 'Ms', Msink
        print 'Mdot_s', Mdot_sc
        print 'Mdot_d', Mdot_d, Mdot_r, Mdot_fl
        plt.figure(f_150.number)
        plt.plot(time_kyr,Rsnow[:,0],lw=3,c='C'+str(irun),label=label[irun])
#        plt.plot(time_kyr,Rsnow[:,1],c='C'+str(irun),ls='--')
        plt.legend(frameon=False)
        plt.savefig(path_outs[0]+'/history_snowline_150.pdf',bbox_inches='tight')
        if label[irun]=='R_l14_nofeed': continue

        plt.figure(f_M.number)
        plt.plot(time_kyr,Msink,c='C'+str(irun),lw=3,label=label[irun])
        plt.plot(time_kyr,Mdisk_k,lw=3, c='C'+str(irun))
        plt.plot(time_kyr,Mdisk_p,'--',lw=4,c='C'+str(irun))
        plt.legend(frameon=False)
        plt.savefig(path_outs[0]+'/history_M.pdf',bbox_inches='tight')

        plt.figure(f_Mdot.number)
        if irun==1:
            plt.plot(time_kyrc,Mdot_sc*1.e-3,lw=3,c='C'+str(irun),label=label[irun])
            plt.plot(time_kyr[::3],Mdot_fl[::3,1]*1.e-6,'*',c='C'+str(irun),markersize=12)
            print Mdot_fl[::3,1]
        else:
            plt.plot(time_kyrc,Mdot_sc*1.e-3,lw=3,c='C'+str(irun),label=label[irun])
            plt.plot(time_kyr,Mdot_fl[:,1]*1.e-6,'*',c='C'+str(irun),markersize=12)
#        plt.plot(time_kyr,(Mdot_d[:,0]+Mdot_r[:,0])*1.e-6,'--', c='C'+str(irun))
        plt.legend(frameon=False)
        plt.savefig(path_outs[0]+'/history_Mdot.pdf',bbox_inches='tight')

        plt.figure(f_R.number)
        plt.plot(time_kyr,Rk,lw=3,c='C'+str(irun),label=label[irun])
        plt.plot(time_kyr,Rp,'--',lw=4,c='C'+str(irun))
        plt.legend(frameon=False)
        plt.savefig(path_outs[0]+'/history_R.pdf',bbox_inches='tight')


    
#    plt.show() 

 

#############################################################################################
# measure and plot disc properties as function of radius
def trace_radial_properties(path,path_out=None,overwrite=True,order='<',scl=20):
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    istart = 0; iend = -1
    outputs = discload['ioutput_list'][istart:iend]
    Msink=np.asarray(discload['Msink_list'])[istart:iend];
    thres_rho_list = discload['thres_rho_list'][istart:iend]
    center_sink_list=(np.asarray(discload['center_sink_list']))[istart:iend];
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])[istart:iend]
    Mdisc=np.asarray(discload['Mdisc_list'])[istart:iend];
    center_disc_list=(np.asarray(discload['center_disc_list'])[istart:iend]);
    cvel_disc_list=np.asarray(discload['cvel_disc_list'])[istart:iend]
    ax_ctot = np.asarray(discload['ax_ctot_list'])[istart:iend]
    r_ctot = np.asarray(discload['r_ctot_list'])[istart:iend]
    center_tot = (center_sink_list*Msink[:,np.newaxis]+center_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    cvel_tot = (cvel_sink_list*Msink[:,np.newaxis]+cvel_disc_list*Mdisc[:,np.newaxis])/(Msink+Mdisc)[:,np.newaxis]
    print 'outputs', outputs
    nout = len(outputs)
    time_list = discload['time_list'][istart:iend]
    Sigma_t = []
    Mdot_r_t = []
    Mdot_v_t = []
    Source_t = []
    Sourcein_t = []
    Omega_t = []
    Omega_k_t = []
    Qtoom_t = []
    beta_t = []
    h_r_t = []
    alpha_Rey_t = []
    alpha_Max_t = []
    alpha_Maxz_t = []
    alpha_Poi_t = []
    time_t = []
    SigmaR2W = []; R2Omega = []
    for iout in range(nout)[:]:
#        np.savez(path_out+'/grids_'+str(ioutput).zfill(5),m=grid_m, V=grid_V, mrho=grid_mrho, rho_v=grid_rho_v, rho_m=grid_rho_m, mvr=grid_mvr, mvphi=grid_mvphi, mvz=grid_mvz, mvzin=grid_mvzin, mvx=grid_mvx, mvy=grid_mvy, mvphivzin=grid_mvphivzin, VBr=grid_VBr, VBr2=grid_VBr2, VBphi=grid_VBphi, VBphi2=grid_VBphi2, VBz=grid_VBz, VBz2=grid_VBz2, VBrBphi=grid_VBrBphi, Vgr=grid_Vgr, Vgphi=grid_Vgphi, Vgrgphi=grid_Vgrgphi,VP=grid_VP, VPmag=grid_VPmag,threshold=thres_rho,mz2=grid_mz2,pp=pp,rr=rr)
        ioutput = outputs[iout] 
        gridfile = path_out+'/grids_'+str(ioutput).zfill(5)+'.npz'
        if not os.path.exists(gridfile): continue
        time_t.append(time_list[iout])
        grids = np.load(gridfile)
        rr = grids['rr']; dr = np.diff(rr); rc = 0.5*(rr[:-1]+rr[1:])
        Surf = np.pi*np.diff(rr**2)
        pp = grids['pp']
        M_shell = np.sum(grids['m'],axis=1)
        V_shell = np.sum(grids['V'],axis=1)
        Mdot_r = np.sum(grids['mvr'],axis=1)/dr[:,np.newaxis] / boxlen
        print Mdot_r[:,1] * mass_sol /time_Myr
        Sigma = M_shell/Surf[:,np.newaxis]
        Rho_mean = M_shell / V_shell
        Source = np.sum(grids['mvz'][:,:,0],axis=1) / V_shell[:,0] # rho*vz
        Sourcein = np.sum(grids['mvzin'][:,:,0],axis=1) / V_shell[:,0] # rho*vzin
        Mdot_v = np.cumsum( Source * Surf ) / boxlen
        Jdot_z = np.sum(grids['mvphivzin'][:,:,0],axis=1) / V_shell[:,0] # rho*vphi*vzin
        vphi_mean = np.sum(grids['mvphi'],axis=1)/ M_shell
        Omega = vphi_mean / rc[:,np.newaxis] / boxlen
        Mr = np.append(0, np.cumsum(M_shell[:,1]))+Msink[iout]/mass_sol
        Mrc = np.exp(interp1d(np.log(rr),np.log(Mr))(np.log(rc)))
        Mrc[0] = Msink[iout]/mass_sol
        Omega_k = np.sqrt(Mrc/rc**3)
        VPtherm = np.sum(grids['VP'],axis=1)
        Ptherm = VPtherm / V_shell 
        Csound = np.sqrt(Ptherm/Rho_mean)
        Qtoom = Omega * Csound / Sigma / np.pi / boxlen #* (len_cm*1.e-2)**3 / (time_Myr*1.e6*365*86400)**2 / mass_kg / 6.67e-11 
        Mdvrdvphi = np.sum(grids['mvrvphi'][:,:,1],axis=1) - np.sum(grids['mvr'][:,:,1],axis=1)*np.sum(grids['mvphi'][:,:,1],axis=1)/M_shell[:,1]
        alpha_Rey = Mdvrdvphi / VPtherm[:,1]
        alpha_Max = - np.sum(grids['VBrBphi'][:,:,1],axis=1) / VPtherm[:,1]  / (4.*np.pi)
        alpha_Maxz = - np.sum(grids['VBzBphi'][:,:,1],axis=1) / VPtherm[:,1]  / (4.*np.pi)
        alpha_Poi = np.sum(grids['Vgrgphi'][:,:,1],axis=1) / VPtherm[:,1]  / (4.*np.pi)# * boxlen**2
        VPmag = np.sum(grids['VPmag'],axis=1) / (8.*np.pi)
        beta = VPtherm / VPmag        
        print 'Msink', Msink[iout],'Mdisc', Mdisc[iout]
        h_r = np.sqrt(np.sum(grids['mz2'][:,:,1],axis=1)/M_shell[:,1])
       
        Sigma_t.append(Sigma)
        Mdot_r_t.append(Mdot_r)
        Mdot_v_t.append(Mdot_v)
        Source_t.append(Source)
        Sourcein_t.append(Sourcein)
        Omega_t.append(Omega)
        Omega_k_t.append(Omega_k)
        Qtoom_t.append(Qtoom)
        beta_t.append(beta)
        h_r_t.append(h_r*len_AU)
        alpha_Rey_t.append(alpha_Rey)
        alpha_Max_t.append(alpha_Max)
        alpha_Maxz_t.append(alpha_Maxz)
        alpha_Poi_t.append(alpha_Poi)
        indsig = np.where(Sigma[:,1]>0)
        #print 'Sigma', Sigma[:,1], indsig
        Sigfit = np.polyfit(np.log(rc[indsig]),np.log(Sigma[:,1][indsig]),1)
        print Sigfit
        SigmaR2W.append(Sigma[:,1]*rc**2*(alpha_Rey+alpha_Max+alpha_Poi)*Csound[:,1]**2*boxlen)
        R2Omega.append(rc**2*Omega[:,1]*boxlen**2) 
#        f, ax = plt.subplots(2,1)
#        rc_AU = rc*len_AU
#        ax[0].plot(rc_AU,Sigma*mass_kg/len_cm**2*1.e4, label='Sigma (kg/m^2)')
#        ax[0].plot(rc_AU,-Mdot_r * mass_sol/time_Myr, label='Mdot (Ms/Myr)')
#        ax[0].plot(rc_AU,-Mdot_v * mass_sol/time_Myr, label='Mdot_vertical (Ms/Myr)')
#        ax[0].plot(rc_AU,-(Mdot_r+Mdot_v) * mass_sol/time_Myr, label='Mdot_total (Ms/Myr)')
#        ax[0].plot(rc_AU,-Source * dens_gcc * vel_ms * 1.e3, label='Source (kg/s/m^2)')
#        ax[0].plot(rc_AU,Omega/Omega_k, label='Rotation (Omega/Omega_k)')
#        ax[0].plot(rc_AU,Qtoom, label='Qtoom')
#        ax[0].plot(rc_AU,beta[:,1], label='beta_disc')
#        ax[0].plot(rc_AU,beta[:,0], label='beta_env')
#        ax[0].plot(rc_AU,h_r/rc, label='h/r')
#        ax[0].set_yscale('log')
#        ax[0].legend()
#        plt.plot(rc_AU,alpha_Rey, label='Tur')
#        plt.plot(rc_AU,alpha_Max, label='Mag')
#        plt.plot(rc_AU,alpha_Poi, label='Gra')
#        plt.plot(rc_AU,alpha_Rey+alpha_Max+alpha_Poi, label='alpha')
#        ax[1].legend()
#        ax[1].set_yscale('symlog',linthreshy=1.e-3)
#        print Omega/Omega_k
#        print Qtoom
#        plt.show()
    lbs = 18
    lgs = 16 
    fgs = (6,2.5)
    xlim = [0, 40]
    ntime = len(time_t)

    if not os.path.exists(path_out+'/disc_prop'): os.makedirs(path_out+'/disc_prop')
    rc_AU = rc*len_AU
    fig = plt.figure(figsize=fgs)
    for i in range(ntime): plt.loglog(rc_AU, Sigma_t[i][:,1]*mass_kg/len_cm**2*1.e4, label=str(time_t[i]*1.e3)[:5]+' kyr')
    plt.legend(frameon=False, fontsize=lgs, labelspacing=0.2)
    plt.ylabel(r'$\Sigma$ (kg m$^{-2}$)',fontsize=lbs)
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    plt.xlim([rc_AU[0],40]); plt.ylim([1.e1,1.e4])
    plt.gcf().subplots_adjust(bottom=0.25,right=0.95)
    save_fig(path_out+'/disc_prop/Sigma_r')

    fig,ax = plt.subplots(3,1,sharex=True,figsize=(6,6))
    for i in range(ntime): ax[0].plot(rc_AU,-Mdot_r_t[i] * mass_sol/time_Myr, label=time_t[i])
    ax[0].set_ylabel('Mdot_r (Ms/Myr)')
    ax[0].set_yscale('symlog',linthreshy=1.e-1)
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): ax[1].plot(rc_AU,-Mdot_v_t[i] * mass_sol/time_Myr, '--')
    ax[1].set_ylabel('Mdot_v (Ms/Myr)')
    ax[1].set_yscale('symlog',linthreshy=1.e-1)
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): ax[2].plot(rc_AU,-(Mdot_v_t[i]+Mdot_r_t[i][:,1]) * mass_sol/time_Myr, ':')
    plt.ylabel('Mdot_tot (Ms/Myr)')
    ax[2].set_yscale('symlog',linthreshy=1.e-1)

    fig = plt.figure(figsize=fgs)
    plt.plot(xlim, [0,0], c='gray')
    for i in range(ntime): 
        Mdot_r_in = np.ma.masked_where((Sigma_t[i][:,1]==0.),Mdot_r_t[i][:,1])
        Mdot_v_in = np.ma.masked_where((Sigma_t[i][:,1]==0.),Mdot_v_t[i])
        plt.plot(rc_AU,-(Mdot_v_in+Mdot_r_in) * mass_sol/time_Myr, '--')
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): 
        Mdot_r_in = np.ma.masked_where((Sigma_t[i][:,1]==0.),Mdot_r_t[i][:,1])
        plt.plot(rc_AU,-Mdot_r_in * mass_sol/time_Myr, label=time_t[i])
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): 
#        Mdot_r_out = np.ma.masked_where((Sigma_t[i][:,1]>0.),Mdot_r_t[i][:,0])
        Mdot_r_out = Mdot_r_t[i][:,0]
        plt.plot(rc_AU,-Mdot_r_out * mass_sol/time_Myr, ':')
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    plt.ylabel(r'$\dot{M}$ ($M_\odot$ Myr$^{-1}$)',fontsize=lbs)
    #plt.yscale('symlog',linthreshy=1.e0)
    plt.xlim(xlim); plt.ylim([-50,50])
    plt.gcf().subplots_adjust(bottom=0.25,right=0.95)
    save_fig(path_out+'/disc_prop/Mdot_r')

    fig = plt.figure(figsize=fgs)
    for i in range(ntime): plt.semilogy(rc_AU,-Source_t[i] * dens_gcc * vel_ms * 1.e3, label=time_t[i])
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): plt.semilogy(rc_AU,-Sourcein_t[i] * dens_gcc * vel_ms * 1.e3, ':')
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    plt.ylabel('S (kg s$^{-1}$ m$^{-2}$)',fontsize=lbs)
    plt.xlim(xlim); plt.ylim([1.e-10,3.e-8])
    plt.gcf().subplots_adjust(bottom=0.25,right=0.95,left=0.2)
    save_fig(path_out+'/disc_prop/Source_r')

    #fig = plt.figure(figsize=fgs)
    fig,ax = plt.subplots(2,1,sharex=True,figsize=(6,4))
    for i in range(ntime): ax[0].loglog(rc_AU,Omega_t[i][:,1]/(time_Myr*1.e6*365*86400))
    ax[0].set_prop_cycle(None)
    for i in range(ntime): ax[0].loglog(rc_AU,Omega_k_t[i]/(time_Myr*1.e6*365*86400), '--')
    ax[0].set_prop_cycle(None)
    for i in range(ntime): ax[0].loglog(rc_AU,Omega_t[i][:,0]/(time_Myr*1.e6*365*86400),':')
    ax[0].set_ylim([5e-11,5e-7])
    for i in range(ntime): ax[1].semilogx(rc_AU,Omega_t[i][:,1]/Omega_k_t[i], label=time_t[i])
    ax[1].set_prop_cycle(None)
    for i in range(ntime): ax[1].semilogx(rc_AU,Omega_t[i][:,0]/Omega_k_t[i],':')
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    ax[0].set_ylabel(r'$\Omega$ (rad s$^{-1}$)',fontsize=lbs)
    ax[1].set_ylabel(r'$\Omega/\Omega_\mathrm{k}$',fontsize=lbs)
    plt.xlim([rc_AU[0],40]); plt.ylim([0.6,1.])
    plt.subplots_adjust(hspace=0)
    plt.gcf().subplots_adjust(bottom=0.15,right=0.95)
    save_fig(path_out+'/disc_prop/Omega_r')

    fig = plt.figure(figsize=fgs)
    for i in range(ntime): plt.semilogy(rc_AU,Qtoom_t[i][:,1], label=time_t[i])
    plt.gca().set_prop_cycle(None)
    for i in range(ntime):
        Q_out = np.ma.masked_where((Sigma_t[i][:,1]>0.),Qtoom_t[i][:,0])
        plt.semilogy(rc_AU,Q_out, ':')
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    plt.ylabel(r'$Q_\mathrm{Toomre}$',fontsize=lbs)
    plt.xlim(xlim); plt.ylim([1.e0,1.e3])
    plt.gcf().subplots_adjust(bottom=0.15,right=0.95)
    save_fig(path_out+'/disc_prop/Qtoomre_r')

    fig = plt.figure(figsize=fgs)
    for i in range(ntime): plt.semilogy(rc_AU,beta_t[i][:,1], label=time_t[i])
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): plt.semilogy(rc_AU,beta_t[i][:,0], ':')
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    plt.ylabel(r'$\beta$',fontsize=lbs)
    plt.xlim(xlim); plt.ylim([1.e-1,1.e5])
    plt.gcf().subplots_adjust(bottom=0.25,right=0.95)
    save_fig(path_out+'/disc_prop/beta_r')

    fig = plt.figure(figsize=(6,2))
    for i in range(ntime): plt.plot(rc_AU,h_r_t[i], label=time_t[i])
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    plt.ylabel(r'$H$ (AU)',fontsize=lbs)
    plt.xlim(xlim); plt.ylim([0,3])
    plt.gcf().subplots_adjust(bottom=0.25,right=0.95)
    save_fig(path_out+'/disc_prop/H_r')

#    fig = plt.figure()
    fig,ax = plt.subplots(4,1,sharey=True,sharex=True,figsize=(6,12))
    for i in range(ntime): ax[1].plot(rc_AU,alpha_Rey_t[i], '--')
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): ax[2].plot(rc_AU,alpha_Max_t[i], ':')
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): ax[2].plot(rc_AU,alpha_Maxz_t[i], '-')
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): ax[3].plot(rc_AU,alpha_Poi_t[i], '-.')
    plt.gca().set_prop_cycle(None)
    for i in range(ntime): ax[0].plot(rc_AU,(alpha_Rey_t[i]+alpha_Max_t[i]+alpha_Poi_t[i]),  label=time_t[i])
    ax[3].set_xlabel(r'$R$ (AU)',fontsize=lbs)
    ax[0].set_ylabel(r'$\alpha$',fontsize=lbs)
    ax[1].set_ylabel(r'$\alpha_\mathrm{tur}$',fontsize=lbs)
    ax[2].set_ylabel(r'$\alpha_\mathrm{mag}$',fontsize=lbs)
    ax[3].set_ylabel(r'$\alpha_\mathrm{grav}$',fontsize=lbs)
    plt.subplots_adjust(hspace=0)
    plt.yscale('symlog',linthreshy=1.e-2)
    plt.gcf().subplots_adjust(top=0.95,bottom=0.1,left=0.2,right=0.95)
    save_fig(path_out+'/disc_prop/alpha_r')

    fig = plt.figure(figsize=fgs)
    plt.plot(xlim, [0,0], c='gray')
    for i in range(ntime): plt.plot(rr[1:-1]*len_AU,np.diff(SigmaR2W[i])/np.diff(R2Omega[i])*2*np.pi*mass_sol/time_Myr)
    plt.xlabel(r'$R$ (AU)',fontsize=lbs)
    plt.ylabel(r'$\dot{M}$ ($M_\odot$ Myr$^{-1}$)',fontsize=lbs)
    plt.xlim(xlim); plt.ylim([-50,50])
    plt.gcf().subplots_adjust(bottom=0.25,right=0.95)
    save_fig(path_out+'/disc_prop/Mdot_alpha_r')

#    plt.show()
#############################################################################################
def calc_disc_properties(path,path_out=None,overwrite=True,order='<',scl=20):
    if path_out == None: path_out = path
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    #[time,Msink,Mdisc, rdisc, zdisc, Sigma_list, Rflux_list,Vfluxout_list,Vfluxin_list,jRflux_list,jVfluxout_list,jVfluxin_list,center_list,cvel_list,axis_list,outputs,thres_rho_list,center_sink_list,cvel_sink_list,center_disc_list,rad_surf_list] = read_pickle(path_out+'/disc_props_csink.pkl',21)
    #[time,Msink,Mdisc, rdisc, zdisc, Sigma_list, Rflux_list,Vfluxout_list,Vfluxin_list,jRflux_list,jVfluxout_list,jVfluxin_list,center_list,cvel_list,axis_list,outputs,thres_rho_list,center_sink_list,cvel_sink_list,center_disc_list] = read_pickle(path_out+'/disc_props_csink.pkl',20)
    [time,Msink,Mdisc, rdisc, zdisc, Sigma_list, Rflux_list,Vfluxout_list,Vfluxin_list,jRflux_list,jVfluxout_list,jVfluxin_list,center_list,cvel_list,axis_list,outputs,thres_rho_list,center_sink_list,cvel_sink_list,center_disc_list,rad_surf_list] = read_pickle(path_out+'/disc_props_csink.pkl',21)
    nout = len(outputs)
    f_time, ax_t = plt.subplots(6,1,figsize=(6,6))
    for iout in range(nout)[:]:
        ioutput=outputs[iout]
        nr = len(Rflux_list[iout])
        rr = np.linspace(0, rdisc[iout]*2/len_AU,2*nr+1)
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        amr = ro.amr_source(["rho","vel","P","phi","g","Br","Bl"])
#        print "disc parameters", center_sink_list[iout],cvel_sink_list[iout],axis_list[iout],thres_rho_list[iout]
#        return
        center = center_sink_list[iout]
        cvel = cvel_sink_list[iout]
#        center = center_list[iout]
#        center = (Msink*center_sink_list[iout]+Mdisc*centert_disc_list[iout])/(Msink+Mdisc)      
#        cvel = cvel_list[iout]
        vec_z = axis_list[iout] / np.linalg.norm(axis_list[iout])
        vec_z = np.array([-0.08404289, -0.00555191, -0.99644667])
        disc_region = Cylinder(center,vec_z,rdisc[iout]/len_AU*2,zdisc[iout]/len_AU*5.) 
        amr_disc = RegionFilter(disc_region,amr)
        cell_disc = CellsToPoints(amr_disc, smallest_cell_level=scl).flatten()
        vel_disc = cell_disc["vel"]-cvel[np.newaxis,:]
        post_disc = cell_disc.points-center[np.newaxis,:]
        dist_disc = np.linalg.norm(post_disc,axis=1) ##3D distance
        vec_r = post_disc / dist_disc[:,np.newaxis]
        vec_phi = np.cross(vec_z,vec_r)
        vec_phi = vec_phi / np.linalg.norm(vec_phi,axis=1)[:,np.newaxis]
        vec_rc = np.cross(vec_phi,vec_z)
        vec_rc = vec_rc / np.linalg.norm(vec_rc,axis=1)[:,np.newaxis]
        distc_disc = np.sum(post_disc*vec_rc,axis=1) ##horizontal distance
        z_disc = np.sum(post_disc*vec_z,axis=1)
        distc_disc = np.sqrt( dist_disc**2 - z_disc**2)
        vel_r = np.sum(vel_disc * vec_rc, axis=1)
        vel_phi = np.sum(vel_disc * vec_phi, axis=1)
        vel_z = np.sum(vel_disc * vec_z, axis=1)
        B_disc = (cell_disc["Br"]+cell_disc["Bl"])*0.5
        B_r = np.sum(B_disc * vec_rc, axis=1)
        B_phi = np.sum(B_disc * vec_phi, axis=1)
        B_z = np.sum(B_disc * vec_z, axis=1)
        
        dV_disc = cell_disc.get_sizes()**3
        dm_disc = cell_disc["rho"]*dV_disc

#        xp = vec_rc[0,:]
#        yp = vec_phi[0,:]
        yp = np.cross(vec_z, (center_disc_list[iout]-center))
        yp = yp/np.linalg.norm(yp)
        xp = np.cross(yp,vec_z)
        xp = xp/np.linalg.norm(xp)
        x_disc = np.sum(post_disc*xp[np.newaxis,:],axis=1)
        y_disc = np.sum(post_disc*yp[np.newaxis,:],axis=1)
        phi = np.arctan2(y_disc,x_disc)
        pp = np.linspace(-np.pi,np.pi, 65)
        grid_x = lambda x: np.histogram2d(distc_disc, phi, bins=[rr,pp], weights=x)[0]
        grid_m = grid_x(dm_disc)
        grid_V = grid_x(dV_disc)
        grid_PV = grid_x(dV_disc*cell_disc["P"])
        mean_c = np.sqrt(np.sum(grid_PV,axis=1)/np.sum(grid_m,axis=1))
        grid_m2 = grid_x(cell_disc["rho"]**2*dV_disc)
        grid_rho_v = grid_m/grid_V
        grid_rho_m = grid_m2/grid_m
        grid_rho = grid_rho_m
        grid_rho = np.ma.masked_where(np.isnan(grid_rho), grid_rho)
        sig_rho = np.sqrt(grid_m2/grid_V - grid_rho**2)
        grid_mvr = grid_x(dm_disc*vel_r)
        grid_mvphi = grid_x(dm_disc*vel_phi)
        #grid_mvz = grid_x(dm_disc*vel_z)
        grid_mvr2 = grid_x(dm_disc*vel_r**2)
        grid_mvphi2 = grid_x(dm_disc*vel_phi**2)
        #grid_mvz2 = grid_x(dm_disc*vel_z**2)
        sig_vr = np.sqrt(grid_mvr2/grid_m - (grid_mvr/grid_m)**2)
        sig_vphi = np.sqrt(grid_mvphi2/grid_m - (grid_mvphi/grid_m)**2)
        mean_vr = np.sum(grid_mvr,axis=1)/np.sum(grid_m,axis=1)
        mean_vphi = np.sum(grid_mvphi,axis=1)/np.sum(grid_m,axis=1)
        #sig_vz = np.sqrt(grid_mvz2/grid_m - (grid_mvz/grid_m)**2)
        grid_VBr = grid_x(dV_disc*B_r)
        grid_VBr2 = grid_x(dV_disc*B_r**2)
        grid_VBphi = grid_x(dV_disc*B_phi)
        grid_VBphi2 = grid_x(dV_disc*B_phi**2)
        grid_vr = grid_mvr/grid_m
        grid_vr = np.ma.masked_where(np.isnan(grid_vr), grid_vr)
        grid_vphi = grid_mvphi/grid_m
        grid_vphi = np.ma.masked_where(np.isnan(grid_vphi), grid_vphi)
        sig_Br = np.sqrt(grid_VBr2/grid_V - (grid_VBr/grid_V)**2)
        sig_Bphi = np.sqrt(grid_VBphi2/grid_V - (grid_VBphi/grid_V)**2)
        fig, axes = plt.subplots(1,3, subplot_kw=dict(polar=True),figsize=(8,3))
        P_grid, R_grid = np.meshgrid(pp,rr)
        pp_c = 0.5*(pp[:-1]+pp[1:]); rr_c = 0.5*(rr[:-1]+rr[1:])
        #Pc_grid, Rc_grid = np.meshgrid(pp_c,rr_c)
        dphi = (pp[1]-pp[0])*0.5
        Pc_grid, Rc_grid = np.meshgrid(pp+dphi,rr_c)
        vlim = np.maximum(abs(grid_vr).max(), abs(grid_vphi).max())/2.
        axes[0].pcolormesh(P_grid,R_grid,np.log10(grid_rho))
        axes[0].contour(Pc_grid+dphi,Rc_grid,np.log10(np.hstack((grid_rho,grid_rho[:,0,np.newaxis]))),colors='k')
        axes[1].pcolormesh(P_grid,R_grid,grid_vr,cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[1].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vr,grid_vr[:,0,np.newaxis])),colors='k')
        axes[2].pcolormesh(P_grid,R_grid,grid_vphi,cmap='seismic_r',vmin=-vlim,vmax=vlim)
        axes[2].contour(Pc_grid+dphi,Rc_grid,np.hstack((grid_vphi,grid_vphi[:,0,np.newaxis])),colors='k')
        save_fig(path_out+'/disc_visu_faceon_'+str(ioutput).zfill(5),tight=False,ps=True)
        dphi = np.pi/100.
        f, ax = plt.subplots(5,1,figsize=(6,6)) ## plot azimuthal variation for each output
        for ir in range(2,nr/2):
        #    print grid_mvr[ir,:].shape, grid_m[ir,:].shape, mean_vr.shape, mean_c.shape,sig_vr[ir,:].shape
            #ax[0].errorbar(pp[1:]+dphi*ir, grid_mvr[ir,:]/grid_m[ir,:]*vel_ms,yerr=sig_vr[ir,:]*vel_ms, capsize=3)
            #ax[1].errorbar(pp[1:]+dphi*ir, grid_mvphi[ir,:]/grid_m[ir,:]*vel_ms,yerr=sig_vphi[ir,:]*vel_ms, capsize=3)
            ax[0].errorbar(pp[1:]+dphi*ir, (grid_mvr[ir,:]/grid_m[ir,:]-mean_vr[ir])/mean_c[ir],yerr=sig_vr[ir,:]/mean_c[ir], capsize=3)
            ax[1].errorbar(pp[1:]+dphi*ir, (grid_mvphi[ir,:]/grid_m[ir,:]-mean_vphi[ir])/mean_c[ir],yerr=sig_vphi[ir,:]/mean_c[ir], capsize=3)
            ax[3].plot(pp[1:]+dphi*ir, (grid_mvr[ir,:]/grid_m[ir,:]-mean_vr[ir])/mean_c[ir])
            ax[4].plot(pp[1:]+dphi*ir, (grid_mvphi[ir,:]/grid_m[ir,:]-mean_vphi[ir])/mean_c[ir])
            #ax[2].errorbar(pp[1:]+dphi*ir, grid_mvz[ir,:]/grid_m[ir,:]*vel_ms,yerr=sig_vz[ir,:]*vel_ms, capsize=3)
            ax[2].errorbar(pp[1:]+dphi*ir, grid_rho[ir,:]*dens_gcc,yerr=sig_rho[ir,:]*dens_gcc, capsize=3)
            #ax[3].errorbar(pp[1:]+dphi*ir, grid_VBr[ir,:]/grid_V[ir,:]*mag_gauss,yerr=sig_Br[ir,:]*mag_gauss, capsize=3)
            #ax[4].errorbar(pp[1:]+dphi*ir, grid_VBphi[ir,:]/grid_V[ir,:]*mag_gauss,yerr=sig_Bphi[ir,:]*mag_gauss, capsize=3)
        #ax[2].set_ylim([-1000,1000])
        #ax[2].set_ylim([0.,1.e-12])
        ax[0].set_ylabel(r'$v_r$ (m/s)')
        ax[1].set_ylabel(r'$v_\phi$ (m/s)')
        ax[2].set_ylabel(r'$\rho$ (g/cm^3)')
        ax[3].set_ylabel(r'$B_r$ (G)')
        ax[4].set_ylabel(r'$B_\phi$ (G)')
        save_fig(path_out+'/vr_vphi_rho_Br_Bphi_'+str(ioutput).zfill(5),tight=False,ps=True)
        
        ## calculate azimuthally averaged quatities
        dd = [0., thres_rho_list[iout],1.e50]
        rr2 = np.arange(0,rdisc[iout]/len_AU,1./2**lmax)
        sum_x = lambda x: np.cumsum(np.histogram2d(distc_disc,cell_disc["rho"], bins=[rr2,dd], weights=x)[0][:,::-1],axis=1)[:,::-1]
        sum_m = sum_x(dm_disc)
        sum_mvr = sum_x(dm_disc*vel_r)
        sum_mvphi = sum_x(dm_disc*vel_phi)
        #sum_mvz = sum_x(dm_disc*vel_z)
        sum_mvrvphi = sum_x(dm_disc*vel_r*vel_phi)
        sum_mvr2 = sum_x(dm_disc*vel_r**2)
        sum_mvphi2 = sum_x(dm_disc*vel_phi**2)
        #sum_mvz2 = sum_x(dm_disc*vel_z**2)
        sum_Vrhodvrdvphi = sum_mvrvphi - sum_mvr*sum_mvphi/sum_m
        sum_VBrBphi = sum_x(dV_disc*B_r*B_phi)
        PV_disc = dV_disc * cell_disc["P"]
        sum_VP = sum_x(PV_disc)
        sum_V = sum_x(dV_disc)
        alpha_Rey = sum_Vrhodvrdvphi/sum_VP
        alpha_Max = -sum_VBrBphi/sum_VP / (4.*np.pi)
        alpha_stress = alpha_Rey + alpha_Max
        sum_rhodvr2 = sum_mvr2 - sum_mvr**2/sum_m
        sum_rhodvphi2 = sum_mvphi2 - sum_mvphi**2/sum_m
        #sum_rhodvz2 = sum_mvz2 - sum_mvz**2/sum_m
        sum_VPmag = sum_x(dV_disc * np.sum(B_disc**2,axis=1)) / (8.*np.pi)
        surf_disc = np.diff(rr**2) * np.pi
        cs2_disc = sum_VP/sum_m
#        rr_cen = 0.5*(rr[1:]+rr[:-1])
#        Omega_disc = sum_mvphi/sum_m/rr_cen[:,np.newaxis]
#        Qtoomre = Omega_disc*np.sqrt(cs2_disc)*surf_disc[:,np.newaxis]/sum_m * vel_ms**2 * len_cm*1.e-2 / mass_kg / np.pi / 6.67e-11
#        beta = sum_VP / sum_VPmag
#        Mdotr = sum_mvr / sum_V * 4. * np.pi * rr_cen[:,np.newaxis]*zdisc[iout]/len_AU * mass_sol / time_Myr * 1.e-6
#        pdb.set_trace()
        #print 4*np.pi*rr[1:]*zdisc[iout]/len_AU
        #print rad_surf_list[iout]
        #print 4*np.pi*rr[1:]*zdisc[iout]/len_AU/ rad_surf_list[iout]
#        print alpha_stress, np.sum(sum_Vrhodvrdvphi,axis=0)/np.sum(sum_VP,axis=0), -np.sum(sum_VBrBphi,axis=0)/np.sum(sum_VP,axis=0)
        #f, ax = plt.subplots(5,1,figsize=(6,6))
#        print Mdotr[:,0], 'Mdot'
#        print Qtoomre[:,0], 'Q'
#        print beta[:,0], 'beta'
        plt.figure(f_time.number)
#        ax_t[1].plot(rr[1:]*len_AU, alpha_stress[:,0], label=time[iout])
        ax_t[2].plot(rr2[1:]*len_AU, alpha_Rey[:,0], label=time[iout])
#        ax_t[3].plot(rr[1:]*len_AU, alpha_Max[:,0], label=time[iout])
        #ax_t[0].plot(rr[1:nr+1]*len_AU, Rflux_list[iout]*4*np.pi*rr[1:nr+1]*zdisc[iout]*(len_cm/1.e2)**2/len_AU*86400*365/2.e30, label=time[iout])
#        ax_t[0].plot(rr[1:]*len_AU, Mdotr[:,0], label=time[iout])
#        ax_t[4].semilogy(rr[1:]*len_AU, Qtoomre[:,0])
#        ax_t[5].semilogy(rr[1:]*len_AU, beta[:,0])
        #ax_t[4].semilogy(rr[1:]*len_AU, sum_rhodvr2/sum_VP)
        #ax_t[5].semilogy(rr[1:]*len_AU, sum_rhodvphi2/sum_VP)
        #ax[6].semilogy(rr[1:]*len_AU, sum_rhodvz2/sum_VP)
#        ax[0].plot(rr[1:]*len_AU, sum_V*len_AU**3, label='volume')
#        ax[1].plot(rr[1:]*len_AU, sum_m*mass_sol, label='mass')
#        ax[2].plot(rr[1:]*len_AU, sum_mvr/sum_m*vel_ms, label='vr')
#        ax[2].plot(rr[1:]*len_AU, -sum_mvr/sum_m*vel_ms, label='-vr')
#        ax[3].plot(rr[1:]*len_AU, sum_mvphi/sum_m*vel_ms, label='vphi')
#        ax[3].plot(rr[1:]*len_AU, np.sqrt(6.67e-11*(np.cumsum(sum_m)+Msink[iout]/mass_sol)*mass_kg/rr[1:]/len_cm*100),label='vk')
    plt.yscale('log')
    ax_t[0].set_yscale('symlog',linthreshy=1.e-8)
    ax_t[1].set_yscale('symlog',linthreshy=1.e-3)
    ax_t[2].set_yscale('symlog',linthreshy=1.e-3)
    ax_t[3].set_yscale('symlog',linthreshy=1.e-3)
    ax_t[0].legend()
    ax_t[0].set_ylabel(r'$\dot{M} (M_\odot/yr)$')
    ax_t[1].set_ylabel(r'$\alpha$')
    ax_t[2].set_ylabel(r'$\alpha_\mathrm{Reynolds}$')
    ax_t[3].set_ylabel(r'$\alpha_\mathrm{Maxwell}$')
    ax_t[4].set_ylabel(r'$Q$')
    ax_t[5].set_ylabel(r'$\beta$')
    save_fig(path_out+'/alpha_transport',tight=False,ps=True)
    plt.show()
        
#############################################################################################
def disc_characteristics(path,path_out=None,overwrite=True,order='<',scl=20):
    if path_out == None: path_out = path
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    unit_accretion = mass_sol/time_Myr
    unit_flux = dens_gcc * vel_ms * 1.e3 # kg/m^2/s
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    [time,Msink,Mdisc, rdisc, zdisc, Sigma_list, Rflux_list,Vfluxout_list,Vfluxin_list,jRflux_list,jVfluxout_list,jVfluxin_list,center_list,cvel_list,axis_list,outputs] = read_pickle(path_out+'/disc_props.pkl',16) 
    nout = len(outputs)
    f, ax = plt.subplots(4,1,figsize=(8,8))
    for iout in  range(nout)[2::7]:
        ioutput=outputs[iout]
        #ro = pymses.RamsesOutput(path,ioutput,order=order)
        #amr = ro.amr_source(["rho","vel","phi","g","Br","Bl"])
        nr = len(Rflux_list[iout])
        rr = np.linspace(0, rdisc[iout],nr+1)[1:]
        col = 'C'+str(iout%10)
        ax[0].plot(rr,Rflux_list[iout],label='%i yr'%(time[iout]),c=col)
        ax[0].plot(rr,Vfluxout_list[iout],':',c=col)
        ax[0].plot(rr,Vfluxin_list[iout],'--',c=col)
        ax[0].set_yscale('symlog')
        #plt.ylim([10.,1.e2])
        ax[0].set_xlabel('Radius (AU)',fontsize=18)
        ax[0].set_ylabel('Mass flux (kg m$^{-2}$ s$^{-1}$)',fontsize=18)
        ax[0].tick_params(axis='both', which='major', labelsize=15)
        ax[0].legend()

        ax[1].plot(rr,jRflux_list[iout],label='%i yr'%(time[iout]),c=col)
        ax[1].plot(rr,jVfluxout_list[iout],':',c=col)
        ax[1].plot(rr,jVfluxin_list[iout],'--',c=col)
        ax[1].set_yscale('symlog')
        #plt.ylim([10.,1.e2])
        ax[1].set_xlabel('Radius (AU)',fontsize=18)
        ax[1].set_ylabel('J flux (kg s$^{-2}$)',fontsize=18)
        ax[1].tick_params(axis='both', which='major', labelsize=15)
        ax[1].legend()
    plt.show()
    save_fig(path+'/Mass_flux',tight=True)
#        f = plt.figure(figsize=(8,4))
    return
#############################################################################################
def empirical_func(t, M0, t0, tau, n):
    return M0 * (1-np.exp(-((t-t0)/tau)**n))
def Mdot_func(t, M0, t0, tau, n):
    return M0 * np.exp(-((t-t0)/tau)**n) * ((t-t0)/tau)**(n-1) * n/tau
def empirical_fit(time,var):
    tau0 = interp1d(var,time-time[0])(0.5*(var[-1]+var[0]))
    p0 = [var[-1],time[0],tau0*0.8,0.7]
    print p0
    bounds = ([p0[0]*0.1,p0[1]*0.1,p0[2]*0.1,0.1],[p0[0]*10,p0[1],p0[2]*40,3.])
    print bounds
    popt, pcov = curve_fit(empirical_func,time,var,p0=p0,bounds=bounds)
    return popt, pcov
def empirical_func2(t, M0, t0, tau, n):
    return M0 * np.exp(-((t-t0)/tau)**n)
def empirical_fit2(time,var):
    tau0 = interp1d(var,time-time[0])(0.5*(var[-1]+var[0]))
    p0 = [var[0],time[0],tau0,0.7]
    bounds = ([p0[0]*0.1,p0[1]*0.1,p0[2]*0.1,0.1],[p0[0]*10,p0[1],p0[2]*10,0.9])
    print p0
    print bounds
    popt, pcov = curve_fit(empirical_func2,time,var,p0=p0,bounds=bounds)
    return popt, pcov
#############################################################################################
def simple_disc_model(path):
    [time,Msink,Mdisc, rdisc, zdisc, Sigma_list, Rflux_list,Vfluxout_list, Vfluxin_list] = read_pickle(path+'/disc_props_csink.pkl',9)
    Mtot = np.asarray(Msink) + np.asarray(Mdisc)
    #Mfunc = lambda t, M0, t0, tau, n: M0*(1-np.exp(-((t-t0)/tau)**n))
    ## fit Mtot
#    fit_tot, pcov = empirical_fit(time,Mtot)
#    fit_sink, pcov = empirical_fit(time,Msink)
#    fit_r, pcov = empirical_fit(time[:-2],rdisc[:-2])
#    Mtot_emp = empirical_func(time,fit_tot[0],fit_tot[1],fit_tot[2],fit_tot[3])
#    Msink_emp = empirical_func(time,fit_sink[0],fit_sink[1],fit_sink[2],fit_sink[3])
#    Mdisc_emp = Mtot_emp - Msink_emp
#    Rdisc_emp = empirical_func(time,fit_r[0],fit_r[1],fit_r[2],fit_r[3])
#    fit_disc, pcov = empirical_fit2(time,Mdisc)
#    Mdisc_emp2 = empirical_func2(time,fit_disc[0],fit_disc[1],fit_disc[2],fit_disc[3])
#    Msink_emp2 = Mtot_emp - Mdisc_emp2
#    print 'Mtot', fit_tot
#    print 'Msink', fit_sink
#    print 'R', fit_r
#    print 'Mdisc', fit_disc

    f = plt.figure(figsize=(8,4))
    plt.semilogy(time,rdisc,c='C3',lw=2)
#    plt.semilogy(time,Rdisc_emp,':',c='C3') 
    plt.ylim([10.,1.e2])
    plt.xlabel('Time (year)',fontsize=18)
    plt.ylabel('Disk radius (AU)',fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=15)
    save_fig(path+'/Rdisc_time',tight=True) 

    f = plt.figure(figsize=(8,4))
    plt.plot(time,Mtot,c='C0',lw=2,label=r'$M_\mathrm{tot}$')
#    plt.plot(time,Mtot_emp,':',c='C0')
    plt.plot(time,Msink,'--',c='C1',lw=2,label=r'$M_\mathrm{sink}$')
#    plt.plot(time,Msink_emp2,':',c='C1')
    plt.plot(time,Mdisc,'-.',c='C2',lw=2,label=r'$M_\mathrm{disk}$')
#    plt.semilogy(time,Mdisc_emp2,':',c='C2')
    plt.legend(fontsize=18)
    plt.ylim([1.e-3,1.])
    plt.xlabel('Time (year)',fontsize=18)
    plt.ylabel('Mass ($M_\odot$)',fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=15)
    save_fig(path+'/Mdisc_time',tight=True)

#    f = plt.figure(figsize=(8,4))
#    plt.semilogy(time,Mdot_func(time,fit_tot[0],fit_tot[1],fit_tot[2],fit_tot[3]),c='C4',label='Mdot')
#    plt.ylim([1.e-6,1.e-4])
#    plt.xlabel('Time (year)',fontsize=18)
#    plt.ylabel('Mass accretion rate($M_\odot$/year)',fontsize=18)
#    plt.tick_params(axis='both', which='major', labelsize=15)
#    save_fig(path+'/Mdot_time',tight=True)

    f = plt.figure(figsize=(8,10))
    plt.plot(time,Mtot,lw=2,label='Mtot')
#    plt.plot(time,Mtot_emp)
    plt.plot(time,Msink,'--',lw=2,label='Msink')
#    plt.plot(time,Msink_emp2,'--')
    plt.plot(time,Mdisc,lw=2,label='Mdisc')
   # plt.plot(time,Mdisc_emp)
#    plt.plot(time,Mdisc_emp2)
    plt.plot(time,rdisc,':',lw=2,label='Rdisc')
#    plt.plot(time,Rdisc_emp,':')
#    plt.semilogy(time,Mdot_func(time,fit_tot[0],fit_tot[1],fit_tot[2],fit_tot[3]),label='Mdot')
    plt.legend()
    plt.show()
#############################################################################################
def time_series(path,path_out=None,overwrite=True,order='<',disc_thres=1.e-13,core_thres=1.e-11):
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    nr = 40
    outputs = search_ro(path)
    nout = len(outputs)
    f, ax = plt.subplots(1,2)
    f.suptitle(titlepath)

    for iout in range(nout)[::100]:
        ioutput = outputs[iout]
        filename = path_out+'/rho_func_'+str(ioutput).zfill(5)+'.pkl'
        if not os.path.exists(filename):  continue
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        time=ro.info['time'] * time_Myr
        func_M_rho, func_dm_drho, func_d2m_drho2, rho_lims = rho_functions(filename)
        rho_list = np.logspace(np.log10(rho_lims[0])-0.3, np.log10(rho_lims[1])+0.1, nr)
        #r_list = np.logspace(np.log10(r_lims[0])+0.9, np.log10(r_lims[1])-0.5, nr)
        ax[0].loglog(rho_list*dens_gcc, func_M_rho(rho_list)*mass_sol,label=str(time)[:5]+'(Myr)')
        #rho_r_list = func_rho_r(r_list)
        #ax[1].loglog( r_list*len_AU, rho_r_list*dens_gcc, label=str(time)[:5]+'(Myr)')
        #ax[1].loglog( r_list*len_AU, func_M_rho(rho_r_list)*mass_sol,label=r'$M(>\rho)$')
        #ax[1].loglog( r_list*len_AU, func_j_rho(rho_r_list)*dens_gcc*vel_ms*1.e3,label=r'$j$')

    ax[0].set_xlabel(r'$\rho$(g/cm$^{3}$)')
    ax[0].set_ylabel(r'$M(M_\odot)$')
    ax[1].set_xlabel(r'$R$(AU)')
    ax[1].set_ylabel(r'$\rho$(g/cm$^{3}$)')
    ax[0].legend()
    plt.show()
#############################################################################################
def mass_core_disc(path,num,path_out=None,overwrite=True,order='<',scl=20,center_def='None',rotate=None,c_axis=None,c_vel=None,center_img=[0.5,0.5,0.5],disc_thres=1.e-13,core_thres=1.e-11,make_phi=False):
    ro = pymses.RamsesOutput(path,num,order=order)
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    filename = path_out+'/rho_func_'+str(num).zfill(5)+'.pkl'
    if not os.path.exists(filename):  calc_rho_relations(path,num,path_out=None,overwrite=True,order='<',scl=20,center_def=center_def,rotate=rotate,c_axis=c_axis,c_vel=c_vel,center_img=center_img,disc_thres=1.e-13,core_thres=1.e-11,make_phi=False) 
    ro = pymses.RamsesOutput(path,num,order=order)
    time=ro.info['time'] * time_Myr
    func_M_rho, func_J_rho, func_j_rho, func_dm_drho, func_M_r, func_rho_r, rho_lims, r_lims = rho_functions(filename)
    nr = 40
    rho_list = np.logspace(np.log10(rho_lims[0])-0.3, np.log10(rho_lims[1])+0.1, nr)
    r_list = np.logspace(np.log10(r_lims[0])+0.9, np.log10(r_lims[1])-0.5, nr)


#    if center_def=='dmax': center_img = find_center_dmax(ro, c_center=center_img,c_radius=0.005,scl=lmax)
#    elif center_def=='baricenter': center_img = find_baricenter(ro, c_center=center_img,c_radius=200./len_AU,scl=lmax)
#    if c_axis is None and c_vel is None:
#        if rotate == 'AM': c_axis, c_vel = find_axis(ro, c_center=center_img,c_radius=30./len_AU,scl=lmax)
#        else: c_axis = None
#    c_center = center_img.copy()

    disc_thres_code = disc_thres/dens_gcc
    core_thres_code = core_thres/dens_gcc
    disc_select_func = lambda dset: (dset["rho"]>=disc_thres_code)*(dset["rho"]<core_thres_code)
    core_select_func = lambda dset: (dset["rho"]>=core_thres_code)

#    amr = ro.amr_source(["rho","vel"])
#    cell_source = CellsToPoints(amr, smallest_cell_level=scl)
#    cells = cell_source.flatten()
#    r = np.sqrt( np.sum((cells.points-c_center)**2,axis=1) )
#    dm = cells["rho"]*cells.get_sizes()**3
#    dj = np.cross( (cells.points-c_center[np.newaxis,:]), (cells["vel"]-c_vel[np.newaxis,:]) )
#    dj = np.sum( dj*c_axis, axis=1 ) / np.linalg.norm(c_axis)
#    sort = np.argsort(cells["rho"])[-1::-1]
#    sort_r = np.argsort(r)
    
#    M_of_rho = np.cumsum( dm[sort] ) * mass_sol
#    J_of_rho = np.cumsum( dm[sort] * dj[sort] ) * mass_kg * vel_ms * len_cm * 1.e-2
#    M_of_r = np.cumsum( dm[sort_r] )
#    func_M_r = interp1d(r[sort_r], M_of_r)
#    rho_of_r = lambda x: derivative(func_M_r,x,dx=x*1e-1)/(4.*np.pi*x**2)

#    rho = cells["rho"][sort] * dens_gcc
#    rho_list = np.logspace(np.log10(rho[-1])+0.1,np.log10(rho[0])-0.1,40)

#    func_M_rho = interp1d(rho,M_of_rho)
#    func_J_rho = interp1d(rho,J_of_rho)
#    func_j_rho = lambda x: derivative(func_J_rho,x,dx=x*1e-1) / derivative(func_M_rho,x,dx=x*1e-1)
#    m_of_rho = derivative(func_M_rho,rho_list,dx=rho_list*1e-2)

#    r_list = np.logspace(np.log10(r[sort_r[0]])+0.9,np.log10(r[sort_r[-1]])-0.5,40)
#    rho_r_list = rho_of_r(r_list)*dens_gcc

    f, ax = plt.subplots(1,2)
    rho_r_list = func_rho_r(r_list)
    ax[1].loglog( r_list*len_AU, rho_r_list*dens_gcc, label=r'$\rho$')
    ax[1].loglog( r_list*len_AU, func_M_rho(rho_r_list)*mass_sol,label=r'$M(>\rho)$')
    ax[1].loglog( r_list*len_AU, func_j_rho(rho_r_list)*dens_gcc*vel_ms*1.e3,label=r'$j$')
    ax[1].legend(loc='best')
    ax[1].set_xlabel(r'$R$(AU)')

    ax[0].loglog(rho_list, -func_dm_drho(rho_list), label=r'$dM/d\rho$')
    ax[0].loglog(rho_list, func_M_rho(rho_list)*mass_sol,label=r'$M(>\rho)$')
    ax[0].loglog(rho_list, func_j_rho(rho_list)*dens_gcc*vel_ms*1.e3,label=r'$j$')
    ax[0].legend(loc='best')
    ax[0].set_xlabel(r'$\rho$(cm$^{-3}$)')
    #time_series(path,f,ax,path_out=path_out,overwrite=True,order='<',disc_thres=1.e-13,core_thres=1.e-11)
#    plt.show()






