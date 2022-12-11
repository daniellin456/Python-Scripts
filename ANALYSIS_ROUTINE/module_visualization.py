#! /usr/bin/env python
## This module contains functions for imaging the ramses output

##### pymses dependences ##################################################################
import pymses
from pymses.filters import PointFunctionFilter
from pymses import RamsesOutput
from pymses.utils import constants as C
from pymses.utils.regions import Sphere, Cylinder, Box
from pymses.filters import CellsToPoints, RegionFilter
from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
from pymses.sources.ramses.tree_utils import octree_compute_neighbors
from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator

##### python dependences ##################################################################
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle as pck
import pdb
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
from matplotlib import ticker
#import multiprocessing as mp

##### module dependences ##################################################################
from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis, normalisation, read_sink_cvs, cal_M_R
#from module_analysis import cylinder_flux, save_pickle, read_pickle
#############################################################################################


###########################################################################
def save_fig(figname,ps=False,tight=False):
    if tight: plt.tight_layout()
    plt.savefig(figname+'.pdf')
    if ps: plt.savefig(figname+'.eps')

###########################################################################
## produce maps for various properties
def colden_map(ro,amr,camera):
    rho_op = ScalarOperator(lambda dset: dset["rho"], ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    map_colden = rt.process(camera, surf_qty=True)
    map_colden = map_colden.map
    return map_colden.T

def rho_map(ro,amr,camera,dis=0.):
    rho_op = ScalarOperator(lambda dset: dset["rho"], ro.info["unit_density"])
    dmap_rho = slicing.SliceMap(amr, camera, rho_op, z=dis) #,verbose=False)
    map_rho = dmap_rho.map
    return map_rho.T

def vel_map(ro,amr,camera,ax1,ax2,dis=0.):
    vx_op = ScalarOperator(lambda dset: np.sum(dset["vel"]*ax1[np.newaxis,:],axis=1), ro.info["unit_velocity"])
    dmap_vx = slicing.SliceMap(amr, camera, vx_op, z=dis) #,verbose=False)
    map_vx=dmap_vx.map
    vy_op = ScalarOperator(lambda dset: np.sum(dset["vel"]*ax2[np.newaxis,:],axis=1), ro.info["unit_velocity"])
    dmap_vy = slicing.SliceMap(amr, camera, vy_op, z=dis) #,verbose=False)
    map_vy=dmap_vy.map
    return map_vx.T, map_vy.T

def velz_map(ro,amr,camera,ax,dis=0.):
    vz_op = ScalarOperator(lambda dset: np.sum(dset["vel"]*ax[np.newaxis,:],axis=1), ro.info["unit_velocity"])
    dmap_vz = slicing.SliceMap(amr, camera, vz_op, z=dis)
    map_vz = dmap_vz.map
    return map_vz.T

def P_map(ro,amr,camera):
    P_op = ScalarOperator(lambda dset: dset["P"] ,  ro.info["unit_pressure"])
    dmap_P = slicing.SliceMap(amr, camera, P_op, z=0.) #,verbose=False)
    map_P = dmap_P.map
    return map_P.T

def T_map(ro,amr,camera):
    T_op = ScalarOperator(lambda dset: dset["T"] ,  ro.info["unit_temperature"])
    dmap_T = slicing.SliceMap(amr, camera, T_op, z=0.) #,verbose=False)
    map_T = dmap_T.map
    return map_T.T

def Bxy_colden_map(ro,amr,camera,ax1,ax2):
    cells = CellsToPoints(amr).flatten()
    print cells["Br"].shape, cells["rho"].shape
    B_op = ScalarOperator(lambda dset: 0.5*(np.sum((dset["Br"]+dset["Bl"])*ax1[np.newaxis,:],axis=-1))*dset["rho"], ro.info["unit_mag"])
    rt = raytracing.RayTracer(amr,ro.info,B_op)
    map_Bcolden = rt.process(camera, surf_qty=True)
    map_Bxcolden = map_Bcolden.map
    B_op = ScalarOperator(lambda dset: 0.5*(np.sum((dset["Br"]+dset["Bl"])*ax2[np.newaxis,:],axis=-1))*dset["rho"], ro.info["unit_mag"])
    rt = raytracing.RayTracer(amr,ro.info,B_op)
    map_Bcolden = rt.process(camera, surf_qty=True)
    map_Bycolden = map_Bcolden.map
    return map_Bxcolden.T, map_Bycolden.T

def Bxy_slice_map(ro,amr,camera,ax1,ax2,dis=0.):
    B_op = ScalarOperator(lambda dset: 0.5*(np.sum((dset["Br"]+dset["Bl"])*ax1[np.newaxis,:],axis=-1)), ro.info["unit_mag"])
    dmap_B = slicing.SliceMap(amr, camera, B_op, z=dis) #,verbose=False)
    map_Bx = dmap_B.map
    B_op = ScalarOperator(lambda dset: 0.5*(np.sum((dset["Br"]+dset["Bl"])*ax2[np.newaxis,:],axis=-1)), ro.info["unit_mag"])
    dmap_B = slicing.SliceMap(amr, camera, B_op, z=dis) #,verbose=False)
    map_By = dmap_B.map
    return map_Bx.T, map_By.T

def Pscal_map(ro,amr,camera,i):
#    print 'I = ',i
#    for dset in amr.iter_dsets():
#        print dset["vel"],dset["vel"].shape
#        print dset["passive"],dset["passive"].shape, dset["passive"][:,:,i], dset["passive"][:,i].shape
    #cells = CellsToPoints(amr).flatten()
    #print cells["passive"][:,i].shape, cells["passive"][:,i]
    Pscal_op = ScalarOperator(lambda dset: dset["rho"]*dset["passive"][:,:,i], ro.info["unit_density"])
    #Pscal_op = ScalarOperator(lambda dset: dset["rho"]*dset["passive"][:,i], ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,Pscal_op)
    dmap_Pscal = rt.process(camera, surf_qty=True)
#    dmap_Pscal = slicing.SliceMap(amr, camera, Pscal_op, z=0.) #,verbose=False)
    map_Pscal = dmap_Pscal.map
    return map_Pscal.T

def MaxLev_map(ro,amr,camera,read=30):
    level_op = MaxLevelOperator()
    amr.set_read_levelmax(read)
    rt = raytracing.RayTracer(amr,ro.info,level_op)
    datamap = rt.process(camera, surf_qty=True)
    map_level = datamap.map
    return map_level.T

def cpu_map(ro,amr,camera):
    cpu_op = ScalarOperator(lambda dset: dset.icpu*(np.ones_like(dset["P"])), ro.info["unit_pressure"])
    rt = raytracing.RayTracer(amr,ro.info,cpu_op)
    datamap = rt.process(camera, surf_qty=True)
    map_cpu = datamap.map
    return map_cpu.T

#############################################################################################
def cal_rot_mat(c_axis,r_axis=None):
    '''calcultate new axes in rotated frame, define rotation matrix to project vectors back to standard frame'''
    c_axis = np.asarray(c_axis)/np.linalg.norm(c_axis)
    if r_axis is not None: #remove parallel part in r_axis
        r_axis = np.asarray(r_axis) 
        r_axis = r_axis - c_axis*np.inner(c_axis,r_axis); r_axis = r_axis/np.linalg.norm(r_axis) 
    if c_axis[2]>=1: 
        cos_th=1.; sin_th=0.
        if r_axis is not None: cos_ph=r_axis[0]; sin_ph=r_axis[1]
        else: cos_ph=1.; sin_ph=0.; r_axis=np.array([1,0,0])
    else:
        cos_th = c_axis[2];sin_th = np.sqrt(1.-cos_th**2)
        if r_axis is not None: cos_ph = r_axis[0]/cos_th; sin_ph=r_axis[1]/cos_th
        else: cos_ph = -c_axis[0]/sin_th; sin_ph = -c_axis[1]/sin_th; r_axis=np.array([cos_ph*cos_th,sin_ph*cos_th,sin_th])
    #rot_x = np.array([cos_th*cos_ph, cos_th*sin_ph, -sin_ph])
    #rot_y = np.array([-sin_ph, cos_ph, 0])
    #mat_Rot = np.vstack((rot_x, rot_y, c_axis))
    u_axis = np.cross(c_axis,r_axis)
    mat_Rot = np.vstack((r_axis, u_axis, c_axis)).T  # matrix transformation from rotated to standard frame
    return c_axis, r_axis, u_axis, mat_Rot

###########################################################################
def define_cameras(center=[0.5,0.5,0.5],radius=0.5,map_size=512,ax1=None,ax2=None,ax3=None):
    '''define 3 cameras with given center, size, resolution, orientation, axes must be perpendicular if asigned '''
    if ax1 is not None:
        ax_z = np.asarray(ax1)/np.linalg.norm(ax1)
        if ax2 is not None: ax_x = np.asarray(ax2); 
        else: 
            if ax_z[2]>=1: ax_x = np.array([1,0,0])
            else: 
                cos_sin = ax_z[2]/np.sqrt(ax_z[0]**2+ax_z[1]**2)
                ax_x = ax_z.copy()
                ax_x[0:2] *= -cos_sin
                ax_x[2] /= cos_sin
        ax_x = ax_x / np.linalg.norm(ax_x)
        if ax3 is not None: ax_y = ax3
        else: ax_y = np.cross(ax_z,ax_x)
    else: ax_x='x'; ax_y='y'; ax_z='z'
    cam_x = Camera(center=center,line_of_sight_axis=ax_x,region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector=ax_z,map_max_size=map_size)
    cam_y = Camera(center=center,line_of_sight_axis=ax_y,region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector=ax_x,map_max_size=map_size)   
    cam_z = Camera(center=center,line_of_sight_axis=ax_z,region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector=ax_y,map_max_size=map_size)
    return [cam_x, cam_y, cam_z]

###########################################################################
def make_image(amr,ro,center,radius,num,path,pos_sink=None,mass_sink=None,i_im=0,overwrite=False,path_out=None,col_dens_only=False,map_size=20,vel_size=20,tag='',ps=False,all_sink=True,ax1=None,ax2=None,PSCAL=[],add_B=False):
## ax1: new z axis, ax2: new x axis
    if( path_out is not None): directory = path_out
    else: directory = path

    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)   
    len_AU = len_pc * 2.e5
    dx = 0.5**lmax
    print 'len_AU', len_AU
    time=ro.info['time'] * time_Myr
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    titre = 'Time = %.3f kyr'%(time*1.e3)

   # name = directory+'/coldens_z'+'_'+str(i_im)+'_'+format(num,'05')+'.pdf'
   # if (os.path.exists(name) and not overwrite): return

    loss = ['x','y','z']
    if ax1 is not None: z_axis, x_axis, y_axis, mat_Rot = cal_rot_mat(ax1,r_axis=ax2) 
    else: x_axis = np.array([1,0,0]); y_axis=np.array([0,1,0]); z_axis=np.array([0,0,1])
    axes = np.vstack((x_axis,y_axis,z_axis))
    print 'axes', x_axis, y_axis, z_axis

    cameras = define_cameras(center=center,radius=radius,map_size=map_size,ax1=z_axis,ax2=x_axis)
    cameras_vel = define_cameras(center=center,radius=radius,map_size=vel_size,ax1=z_axis,ax2=x_axis)
    ind_selects = [[1,2,0],[2,0,1],[0,1,2]]
    if len(PSCAL)>0:
        if PSCAL[-1]==100: PSCAL_only = True; PSCAL = PSCAL[:-1]
        else: PSCAL_only = False
    for los, cam, cam_v, ind in zip(loss, cameras, cameras_vel, ind_selects):
        print 'center', center, 'ind', ind
        cen1 = center[ind[0]]; cen2 = center[ind[1]]; cen3 = center[ind[2]]
        fig_lim = [(cen1-radius-0.5)*len_AU,(cen1+radius-0.5)*len_AU,(cen2-radius-0.5)*len_AU,(cen2+radius-0.5)*len_AU]
        #fig_lim = [-radius*len_AU,radius*len_AU,-radius*len_AU,radius*len_AU]
#        if los=='y': continue
        if len(PSCAL)>0:
           if PSCAL[0] == -1: # plot all passive scalar on same figure
            figname_sum = directory+'/sum_pscal_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
            if not os.path.exists(figname_sum+'.pdf') or overwrite:
             comp = ['passive_'+str(p) for p in PSCAL[1:]]
             for i in PSCAL[1:]:
                 figname = directory+'/'+comp[i]+'_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
                 map_Pscal = Pscal_map(ro,amr,cam,i)
                 if i==PSCAL[1]: map_Pscal_sum = map_Pscal
                 else: map_Pscal_sum += map_Pscal
                 plt.clf()
                 print 'sum pscal', np.sum(map_Pscal)
                 im = plt.imshow(np.log10(map_Pscal*dens_gcc*len_cm)+1.,extent=fig_lim,origin='lower',interpolation='none',zorder=1,vmin=-4,vmax=5)
                 plt.title(titre)
                 plt.xlabel(loss[ind[0]]+' (AU)')
                 plt.ylabel(loss[ind[1]]+' (AU)')
                 cbar = plt.colorbar(im, extend='both')
                 cbar.set_label(r'$\mathrm{log}_{10}\rho (\mathrm{kg}\/\mathrm{m}^{-2})$')
                 save_fig(figname,ps=ps)
             plt.clf()
             im = plt.imshow(np.log10(map_Pscal_sum*dens_gcc*len_cm)+1.,extent=fig_lim,origin='lower',interpolation='none',zorder=1,vmin=-4,vmax=5)
             plt.scatter((cen1-0.5)*len_AU,(cen2-0.5)*len_AU,color='y', marker='*',zorder=2)
             plt.title(titre)
             plt.xlabel(loss[ind[0]]+' (AU)')
             plt.ylabel(loss[ind[1]]+' (AU)')
             cbar = plt.colorbar(im, extend='both')
             cbar.set_label(r'$\mathrm{log}_{10}\rho (\mathrm{kg}\/\mathrm{m}^{-2})$')
             save_fig(directory+'/sum_pscal_'+los+'_'+tag+str(i_im)+'_'+format(num,'05'),ps=ps)

           else:
            comp = ['passive_'+str(p) for p in PSCAL]
            figname = directory+'/'+comp[i]+'_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
            if not os.path.exists(figname+'.pdf') or overwrite:
             #if len(PSCAL)==3: comp = ['dust','gas','CAI']
             #if len(PSCAL)==2: comp = ['1500K','1650K']
             for i in PSCAL:
               plt.clf()
               map_Pscal = Pscal_map(ro,amr,cam,i)
               print 'sum pscal', np.sum(map_Pscal)
               im = plt.imshow(np.log10(map_Pscal*dens_gcc*len_cm)+1.,extent=fig_lim,origin='lower',interpolation='none',zorder=1,vmin=-4,vmax=5)
               plt.title(titre)
               plt.xlabel(loss[ind[0]]+' (AU)')
               plt.ylabel(loss[ind[1]]+' (AU)')
               cbar = plt.colorbar(im, extend='both')
               cbar.set_label(r'$\mathrm{log}_{10}\rho (\mathrm{kg}\/\mathrm{m}^{-2})$')
               save_fig(figname,ps=ps)
           if PSCAL_only: print "plot only pscal, continue!"
           if PSCAL_only: continue

        figname = directory+'/coldens_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
        if not os.path.exists(figname+'.pdf') or overwrite:
            #print 'center', center, 'ind', ind
            #cen1 = center[ind[0]]; cen2 = center[ind[1]]; cen3 = center[ind[2]]
            #fig_lim = [(cen1-radius-0.5)*len_AU,(cen1+radius-0.5)*len_AU,(cen2-radius-0.5)*len_AU,(cen2+radius-0.5)*len_AU]
            map_colden = colden_map(ro,amr,cam)## column density map
            plt.clf()
            im = plt.imshow(np.log10(map_colden*dens_gcc*len_cm)+1.,extent=fig_lim,origin='lower',interpolation='none',vmin=0,vmax=4,cmap='plasma')
            plt.title(titre)
            plt.xlabel(loss[ind[0]]+' (AU)')
            plt.ylabel(loss[ind[1]]+' (AU)')
            cbar = plt.colorbar(im, extend='both') 
            cbar.set_label(r'$\mathrm{log}_{10}\Sigma (\mathrm{kg}\/\mathrm{m}^{-2})$')
            save_fig(figname,ps=ps)   

            if mass_sink is not None:## add sinks
                mask = np.where( (abs(pos_sink[:,ind[0]]-cen1) <= radius)  & (abs(pos_sink[:,ind[1]]-cen2) <= radius) == True )
                if not all_sink :
                    mask = np.where( (abs(pos_sink[:,ind[0]]-cen1) <= radius)  & (abs(pos_sink[:,ind[1]]-cen2) <= radius) & (abs(pos_sink[:,ind[2]]-cen3) <= radius) == True )
                siz_sym=np.sqrt(mass_sink[mask])*50.
                plt.scatter((pos_sink[mask,ind[0]]-0.5)*len_AU,(pos_sink[mask,ind[1]]-0.5)*len_AU,marker='o',s=siz_sym,color='r',zorder=2)
                circle=plt.Circle(((cen1-0.5)*len_AU,(cen2-0.5)*len_AU),dx*len_AU*4,fill=False,edgecolor='w',linestyle=':',zorder=2)
                ax = plt.gca()
                ax.add_artist(circle)
                print 'SINK POSITION', (cen1-0.5)*len_AU,(cen2-0.5)*len_AU, dx*len_AU*10
                plt.xlim(fig_lim[:2])
                plt.ylim(fig_lim[2:]) 
                #plt.show()
                save_fig(figname+'_sink',ps=ps)

        if col_dens_only: continue

        figname = directory+'/rho_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
        if not os.path.exists(figname+'.pdf') or overwrite:
            cen1 = center[ind[0]]; cen2 = center[ind[1]]; cen3 = center[ind[2]]
            map_rho = rho_map(ro,amr,cam) ## density slice map + velocity
            map_vx, map_vy = vel_map(ro,amr,cam_v,axes[ind[0],:],axes[ind[1],:])
            plt.clf()
            pos_vel = np.linspace(-radius,radius,vel_size+1)
            pos_vel = (pos_vel[:-1] + radius/vel_size) * len_AU
            xx,yy = np.meshgrid(pos_vel+(cen1-0.5)*len_AU,pos_vel+(cen2-0.5)*len_AU)
            q = plt.quiver(xx,yy,map_vx*vel_ms,map_vy*vel_ms,color='w',pivot='middle',zorder=2)
            qk = plt.quiverkey(q, 0.8, -0.08, 1.e4, '10 km/s',color='k', labelpos='W')
            im = plt.imshow(np.log10(map_rho*dens_gcc),extent=fig_lim,origin='lower',interpolation='none',zorder=1,cmap='plasma',vmin=-17,vmax=-10)
            plt.title(titre)
            plt.xlabel(loss[ind[0]]+' (AU)')
            plt.ylabel(loss[ind[1]]+' (AU)')
            cbar = plt.colorbar(im, extend='both')
            cbar.set_label(r'$\mathrm{log}_{10}\rho (\mathrm{g}\/\mathrm{cm}^{-3})$')
            save_fig(figname,ps=ps)
 
  
        if True:
          figname_slice = directory+'/B_slice_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
          figname_int = directory+'/B_int_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
          if not os.path.exists(figname+'.pdf') or overwrite:
            map_Bx_slice, map_By_slice = Bxy_slice_map(ro,amr,cam,axes[ind[0],:],axes[ind[1],:])
            map_Bx_colden, map_By_colden = Bxy_colden_map(ro,amr,cam,axes[ind[0],:],axes[ind[1],:])
            map_Bx_int = map_Bx_colden/map_colden; map_By_int = map_By_colden/map_colden
            nvec = 20
            B_size = map_size / nvec
            plt.clf()
            pos_B = np.linspace(-radius,radius,B_size+1)
            pos_B = (pos_B[:-1] + radius/B_size) * len_AU
            xx,yy = np.meshgrid(pos_B+(cen1-0.5)*len_AU,pos_B+(cen2-0.5)*len_AU)
            print xx.shape ,map_Bx_slice[nvec/2::nvec,nvec/2::nvec].shape, 'shapes!'
            q = plt.quiver(xx,yy,map_Bx_slice[nvec/2::nvec,nvec/2::nvec]*mag_gauss,map_By_slice[nvec/2::nvec,nvec/2::nvec]*mag_gauss,color='w',pivot='middle',zorder=2)
       # qk = plt.quiverkey(q, 0.8, -0.08, 1.e4, '10 km/s',color='k', labelpos='W')
            im = plt.imshow(np.log10(np.sqrt(map_Bx_slice**2+map_By_slice**2)*mag_gauss),extent=fig_lim,origin='lower',interpolation='none',zorder=1,cmap='summer',vmin=-3,vmax=0)
            plt.title(titre)
            plt.xlabel(loss[ind[0]]+' (AU)')
            plt.ylabel(loss[ind[1]]+' (AU)')
            cbar = plt.colorbar(im, extend='both')
            cbar.set_label(r'$\mathrm{log}_{10}B (\mathrm{G})$')
            save_fig(figname_slice,ps=ps)

            plt.clf()
            q = plt.quiver(xx,yy,map_Bx_int[nvec/2::nvec,nvec/2::nvec]*mag_gauss,map_By_int[nvec/2::nvec,nvec/2::nvec]*mag_gauss,color='w',pivot='middle',zorder=2)
       # qk = plt.quiverkey(q, 0.8, -0.08, 1.e4, '10 km/s',color='k', labelpos='W')
            im = plt.imshow(np.log10(np.sqrt(map_Bx_int**2+map_By_int**2)*mag_gauss),extent=fig_lim,origin='lower',interpolation='none',zorder=1,cmap='summer',vmin=-3,vmax=0)
            plt.title(titre)
            plt.xlabel(loss[ind[0]]+' (AU)')
            plt.ylabel(loss[ind[1]]+' (AU)')
            cbar = plt.colorbar(im, extend='both')
            cbar.set_label(r'$\mathrm{log}_{10}B (\mathrm{G})$')
            save_fig(figname_int,ps=ps)

        plt.clf()
        figname = directory+'/P_'+los+'_'+tag+str(i_im)+'_'+format(num,'05')
        if not os.path.exists(figname+'.pdf') or overwrite:
            map_P = P_map(ro,amr,cam)
            im = plt.imshow(np.log10(map_P*pressure_P),extent=fig_lim,origin='lower',cmap='YlGnBu',interpolation='none',zorder=1,vmin=-10,vmax=0)
            plt.title(titre)
            plt.xlabel(loss[ind[0]]+' (AU)')
            plt.ylabel(loss[ind[1]]+' (AU)')
            cbar = plt.colorbar(im, extend='both')
            cbar.set_label(r'$\mathrm{log}_{10}P$ (pascal)')
            #plt.show()
            save_fig(figname,ps=ps)

#        plt.clf()
#        map_T = T_map(ro,amr,cam)
#        im = plt.imshow(np.log10(map_T),extent=fig_lim,origin='lower',cmap='BuPu_r',interpolation='none',zorder=1,vmin=1,vmax=4)
#        plt.title(titre)
#        plt.xlabel(loss[ind[0]]+' (AU)')
#        plt.ylabel(loss[ind[1]]+' (AU)')
#        cbar = plt.colorbar(im, extend='both')
#        cbar.set_label(r'$\mathrm{log}_{10}T$ (K)')
#        save_fig(directory+'/T_'+los+'_'+tag+str(i_im)+'_'+format(num,'05'),ps=ps)

#        map_P = P_map(ro,amr,cam) 
#        map_level = MaxLev_map(ro,amr,cam,read=30) ## max refinement level
#        map_cpu = cpu_map(ro,amr,cam) ## cpu domain map
##############################################################################
def vel_profiles(amr,ro,center,radius,num,path,func_M_R,i_im=0,overwrite=False,path_out=None,nbins=30,tag='',ps=False,ax1=None,ax2=None,c_vel=np.array([0,0,0]),scl=20):
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    time=ro.info['time'] * time_Myr
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    titre = titlepath+' '+ str(time)[0:5] +'(Myr)'
    f, fax = plt.subplots(3,3,sharex=True, sharey=True,figsize=(7.5,7.5))
    if( path_out is not None): directory = path_out
    else: directory = path

    figname = directory+'/v_profile'+'_'+str(i_im)+'_'+format(num,'05')
    if (os.path.exists(figname+'.pdf') and not overwrite): return
    sph = Sphere(center, radius)
    amr = RegionFilter(sph,amr)
    if ax1 is not None: z_axis, x_axis, y_axis, mat_Rot = cal_rot_mat(ax1,r_axis=ax2)
    else: x_axis = np.array([1,0,0]); y_axis=np.array([0,1,0]); z_axis=np.array([0,0,1])
    axes = np.vstack((x_axis,y_axis,z_axis))
    cell_source = CellsToPoints(amr, smallest_cell_level=scl)
    cells = cell_source.flatten()

    dm = cells["rho"]*cells.get_sizes()**3
    dr2 = np.sum((cells.points-center[np.newaxis,:])**2,axis=1); dr = np.sqrt(dr2)
    dz = np.sum((cells.points-center[np.newaxis,:])*z_axis[np.newaxis,:],axis=1)
    drc = np.sqrt(dr2-dz**2)
    r_vec = (cells.points-center[np.newaxis,:])/dr[:,np.newaxis]
    rc_vec = r_vec - z_axis * np.sum(r_vec*z_axis[np.newaxis,:],axis=1)[:,np.newaxis]; rc_vec = rc_vec/np.linalg.norm(rc_vec,axis=1)[:,np.newaxis]
    p_vec = np.cross(z_axis[np.newaxis,:],r_vec)
    vel_rel = cells["vel"]-c_vel[np.newaxis,:]
    vz = np.sum(vel_rel*z_axis[np.newaxis,:],axis=1)
    vr = np.sum(vel_rel*r_vec,axis=1)
    vrc = np.sum(vel_rel*rc_vec,axis=1)
    vp = np.sum(vel_rel*p_vec,axis=1) 

    #nbins = 30 
    dz0 = dz.copy()
    dz = np.absolute(dz)
    hist_m, xbins, ybins = np.histogram2d(drc,dz,bins=nbins,weights=dm)
    hist_vp, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vp*dm)
    hist_vrc, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vrc*dm)
    hist_vz, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vz*dm*np.sign(dz0))
    hist_vp_vrc, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=(vp/vrc)*dm)
    hist_rhovz, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vz*dm*cells["rho"]*np.sign(dz0))
    hist_rhovrc, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vrc*dm*cells["rho"])
    hist_rhovr, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vr*dm*cells["rho"])
    hist_jrhovz, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vp*drc*vz*cells["rho"]*dm*np.sign(dz0))
    hist_jrhovrc, xbins, ybins = np.histogram2d(drc,dz,bins=[xbins,ybins],weights=vp*drc*vrc*cells["rho"]*dm)
    #xx , yy = np.meshgrid(xbins*len_AU, ybins*len_AU)
    hist_vp /= hist_m
    hist_vrc /= hist_m
    hist_vz /= hist_m
    hist_rhovr /= hist_m
    hist_rhovz /= hist_m
    hist_rhovrc /= hist_m
    hist_jrhovz /= hist_m
    hist_jrhovrc /= hist_m
    rbins = np.sqrt(xbins[:,np.newaxis]**2 + ybins[np.newaxis,:]**2)
    xbins = xbins*len_AU; ybins = ybins*len_AU
    vmin=-4;vmax=4
    fmax=10.
    fax[0,0].imshow(1.e-3*vel_ms*hist_vp.T,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]),norm=SymLogNorm(linthresh=fmax/1.e2,linscale=1.,vmin=-fmax, vmax=fmax), origin='lower',cmap='PRGn')
    fax[0,1].imshow(-1.e-3*vel_ms*hist_vrc.T,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]),norm=SymLogNorm(linthresh=fmax/1.e2,linscale=1.,vmin=-fmax, vmax=fmax), origin='lower',cmap='PRGn')
    col1=fax[0,2].imshow(-1.e-3*vel_ms*hist_vz.T,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), norm=SymLogNorm(linthresh=fmax/1.e3,linscale=1.,vmin=-fmax, vmax=fmax),origin='lower',cmap='PRGn')
    col2=fax[1,0].imshow(hist_vp.T/abs(hist_vrc.T),extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), origin='lower',norm=LogNorm(vmin=1.e-3,vmax=1.e3),cmap='seismic')
    fmax = 1.e10
    col4=fax[1,1].imshow(hist_jrhovrc.T*vel_ms**2*len_cm*1.e-2*dens_gcc*1.e3,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), origin='lower',norm=SymLogNorm(linthresh=fmax/1.e6,linscale=1.,vmin=-fmax, vmax=fmax),cmap='PiYG')
    fax[1,2].imshow(hist_jrhovz.T*vel_ms**2*len_cm*1.e-2*dens_gcc*1.e3,extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), origin='lower',norm=SymLogNorm(linthresh=fmax/1.e6,linscale=1.,vmin=-fmax, vmax=fmax),cmap='PiYG')
    fmax = 1.e-7
    col3=fax[2,0].imshow(hist_rhovr.T*vel_ms*dens_gcc*1.e3, extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), origin='lower',norm=SymLogNorm(linthresh=fmax/1.e3,linscale=1.,vmin=-fmax, vmax=fmax),cmap='PuOr')
    fax[2,1].imshow(hist_rhovrc.T*vel_ms*dens_gcc*1.e3, extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), origin='lower',norm=SymLogNorm(linthresh=fmax/1.e3,linscale=1.,vmin=-fmax, vmax=fmax),cmap='PuOr')
    fax[2,2].imshow(hist_rhovz.T*vel_ms*dens_gcc*1.e3, extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), origin='lower',norm=SymLogNorm(linthresh=fmax/1.e3,linscale=1.,vmin=-fmax, vmax=fmax),cmap='PuOr')
    fax[0,0].text(xbins[-1]/1.2, ybins[-1]/1.2, r'$v_\phi$')
    fax[0,1].text(xbins[-1]/1.2, ybins[-1]/1.2, r'$-v_r$')
    fax[0,2].text(xbins[-1]/1.2, ybins[-1]/1.2, r'$-v_z$')
    fax[1,0].text(xbins[-1]/1.3, ybins[-1]/1.2, r'$-v_\phi/v_r$')
    fax[2,0].text(xbins[-1]/1.9, ybins[-1]/1.2, 'flux(kg/m$^2$/s)')
    fax[1,1].text(xbins[-1]/2.8, ybins[-1]/1.2, r'$\rho j v_r$(kg/m$^2$/s m$^2$/s)')
    fax[1,2].text(xbins[-1]/2.8, ybins[-1]/1.2, r'$\rho j v_z$(kg/m$^2$/s m$^2$/s)')
    fax[2,1].text(xbins[-1]/1.2, ybins[-1]/1.2, r'$\rho v_r$')
    fax[2,2].text(xbins[-1]/1.2, ybins[-1]/1.2, r'$\rho v_z$')
    xbins = (xbins[:-1]+xbins[1:])*0.5
    ybins = (ybins[:-1]+ybins[1:])*0.5
    rbins = np.sqrt(xbins[:,np.newaxis]**2 + ybins[np.newaxis,:]**2)
    v_grav = np.sqrt(func_M_R(rbins/len_AU)*len_AU/rbins)*boxlen
    scale_v = 30
    for nh in range(0,nbins,nbins/6):
        fax[0,0].plot(xbins,1.e-3*vel_ms*scale_v*hist_vp[:,nh]+ybins[nh],'b')
        fax[0,1].plot(xbins,-1.e-3*vel_ms*scale_v*hist_vrc[:,nh]+ybins[nh],'b')
        fax[0,2].plot(xbins,-1.e-3*vel_ms*scale_v*hist_vz[:,nh]+ybins[nh],'b')
        fax[0,2].plot(-1.e-3*vel_ms*scale_v*hist_vz[nh,:]+xbins[nh],ybins,'c')
        fax[0,0].plot(xbins,1.e-3*vel_ms*scale_v*v_grav[:,nh]+ybins[nh],'gray')
        fax[0,1].plot(xbins,1.e-3*vel_ms*scale_v*v_grav[:,nh]*xbins/rbins[:,nh]+ybins[nh],'gray')
        fax[0,2].plot(xbins,1.e-3*vel_ms*scale_v*v_grav[:,nh]*ybins[nh]/rbins[:,nh]+ybins[nh],'gray')
        fax[0,2].plot(1.e-3*vel_ms*scale_v*v_grav[nh,:]*ybins/rbins[nh,:]+xbins[nh],ybins,'k')
    fax[0,0].set_xlim([0,xbins[-1]])
    fax[0,0].set_ylim([0,ybins[-1]])
    f.suptitle(titre)
    #f.subplots_adjust(wspace=0, hspace=0)
    f.subplots_adjust(0.07,0.07,0.9,0.95,0,0)
    cax1 = f.add_axes([0.9, 0.73, 0.01, 0.2])
    cax2 = f.add_axes([0.9, 0.51, 0.01, 0.2])
    cax3 = f.add_axes([0.9, 0.07, 0.01, 0.2])
    cax4 = f.add_axes([0.9, 0.29, 0.01, 0.2])
    f.colorbar(col1, cax=cax1)
    #tick_locations = np.logspace(np.log10(fmax)-3,np.log10(fmax),4)
    #tick_locations = np.hstack((-tick_locations[::-1],0,tick_locations))
    f.colorbar(col2, cax=cax2)#,ticks=tick_locations,format=ticker.LogFormatterMathtext())
    f.colorbar(col3, cax=cax3)
    f.colorbar(col4, cax=cax4)
    save_fig(figname,ps=ps,tight=False)

    #plt.show()

##############################################################################
def make_image_slice(amr,ro,center,radius,num,path,func_M_R,i_im=0,overwrite=False,path_out=None,col_dens_only=False,map_size=50,tag='',ps=False,all_sink=True,ax1=None,ax2=None,n_slice=3,step_slice=10,c_vel=np.array([0,0,0])):
## ax1: new z axis, ax2: new x axis
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    dx = 0.5**lmax
    time=ro.info['time'] * time_Myr
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    titre = titlepath+' '+ str(time)[0:5] +'(Myr)'
    print 'title', titre
    #plt.clf()
    f, fax = plt.subplots(7,n_slice*2+1,sharex=True, sharey=True,figsize=(2*n_slice+2,7))
    if( path_out is not None): directory = path_out
    else: directory = path

    name = directory+'/v_slice'+'_'+str(i_im)+'_'+format(num,'05')+'.pdf'
    if (os.path.exists(name) and not overwrite): return


    if ax1 is not None: z_axis, x_axis, y_axis, mat_Rot = cal_rot_mat(ax1,r_axis=ax2)
    else: x_axis = np.array([1,0,0]); y_axis=np.array([0,1,0]); z_axis=np.array([0,0,1])
    axes = np.vstack((x_axis,y_axis,z_axis))
    print 'axes', x_axis, y_axis, z_axis

    cameras = define_cameras(center=center,radius=radius,map_size=map_size,ax1=z_axis,ax2=x_axis)
    #cameras_vel = define_cameras(center=center,radius=radius,map_size=vel_size,ax1=z_axis,ax2=x_axis)
    camera = cameras[2]
    fig_lim = [-radius*len_AU,radius*len_AU,-radius*len_AU,radius*len_AU]
    pos_vel = np.linspace(-radius,radius,map_size+1)[:-1]+radius/map_size
    r_list = pos_vel[map_size/2+1:]
    xx,yy = np.meshgrid(pos_vel*len_AU,pos_vel*len_AU)
    map_r = np.sqrt(xx**2+yy**2)

    c_velz = np.inner(c_vel,z_axis)
    c_velx = np.inner(c_vel,x_axis)
    c_vely = np.inner(c_vel,y_axis)
    c_velr = (c_velx*xx + c_vely*yy) / map_r
    c_velt = (c_vely*xx - c_velx*yy) / map_r

    vmin = -4; vmax = -vmin
    fmin = -1.e-7; fmax = -fmin
    for i_slice in range(-n_slice, n_slice+1):
        step = step_slice / len_AU * i_slice
        map_rho = rho_map(ro,amr,camera,dis=step)
        map_vx, map_vy = vel_map(ro,amr,camera,x_axis,y_axis,dis=step)
        map_vz = velz_map(ro,amr,camera,z_axis,dis=step)-c_velz
        map_vr = ( map_vx*xx + map_vy*yy ) / map_r - c_velr
        map_vt = ( map_vy*xx - map_vx*yy ) / map_r - c_velt

        rz_list = np.sqrt(r_list**2+step**2)
        vf_list = np.sqrt(func_M_R(rz_list)/rz_list)*boxlen
        map_vr[map_size/2,map_size/2+1:] = -vf_list[:]
        map_vt[map_size/2,map_size/2+1:] = vf_list[:]        

        fax[0,i_slice+n_slice].imshow(map_rho*dens_gcc,extent=fig_lim,origin='lower',interpolation='none',norm=LogNorm(vmin=1.e-17,vmax=1.e-11))
        fax[1,i_slice+n_slice].imshow(map_vz*vel_ms/1.e3,extent=fig_lim,origin='lower',interpolation='none',vmin=vmin, vmax=vmax,cmap='PRGn')
        fax[2,i_slice+n_slice].imshow(map_vr*vel_ms/1.e3,extent=fig_lim,origin='lower',interpolation='none',vmin=vmin, vmax=vmax,cmap='PRGn')
        col1 = fax[3,i_slice+n_slice].imshow(map_vt*vel_ms/1.e3,extent=fig_lim,origin='lower',interpolation='none',vmin=vmin, vmax=vmax,cmap='PRGn')
        fax[4,i_slice+n_slice].imshow(map_rho*dens_gcc*map_vz*vel_ms*1.e3,extent=fig_lim,origin='lower',interpolation='none', norm=SymLogNorm(linthresh=fmax/1.e3,linscale=1.,vmin=fmin, vmax=fmax), cmap='seismic_r')
        fax[5,i_slice+n_slice].imshow(map_rho*dens_gcc*map_vr*vel_ms*1.e3,extent=fig_lim,origin='lower',interpolation='none', norm=SymLogNorm(linthresh=fmax/1.e3,linscale=1.,vmin=fmin, vmax=fmax), cmap='seismic_r')
        col2 = fax[6,i_slice+n_slice].imshow(map_rho*dens_gcc*map_vt*vel_ms*1.e3,extent=fig_lim,origin='lower',interpolation='none', norm=SymLogNorm(linthresh=fmax/1.e3,linscale=1.,vmin=fmin, vmax=fmax), cmap='seismic_r')


 
        fax[-1,i_slice+n_slice].set_xlabel(r'$z=$'+str(int(step_slice * i_slice))+'(AU)')
        fax[-1,i_slice+n_slice].set_xticklabels([])
    
    #plt.xlabel(loss[ind[0]]+' (AU)')
        #plt.ylabel(loss[ind[1]]+' (AU)'i)
    for  i in range(4): fax[i,0].set_yticklabels([])
    x_ref = xx[map_size/10,map_size/10]; y_ref = yy[map_size/10,map_size/10]
    fax[-1,0].plot([x_ref,x_ref+50],[y_ref,y_ref],lw=2)
    fax[-1,0].text(x_ref+5, y_ref*0.8, '50 AU')
    
    fax[0,0].set_ylabel(r'$\rho$')
    fax[1,0].set_ylabel(r'$v_z$ (km/s)')
    fax[2,0].set_ylabel(r'$v_r$ (km/s)')
    fax[3,0].set_ylabel(r'$v_\phi$ (km/s)')
    fax[4,0].set_ylabel(r'$\rho v_z$(kg/m$^2$/s)')
    fax[5,0].set_ylabel(r'$\rho v_r$(kg/m$^2$/s)')
    fax[6,0].set_ylabel(r'$\rho v_\phi$(kg/m$^2$/s)')

    f.suptitle(titre)
    #f.subplots_adjust(wspace=0, hspace=0)
    f.subplots_adjust(0.06,0.07,0.91,0.95,0,0)
    cax1 = f.add_axes([0.91, 0.54, 0.01, 0.41])
    cax2 = f.add_axes([0.91, 0.07, 0.01, 0.41])
    f.colorbar(col1, cax=cax1)
    tick_locations = np.logspace(np.log10(fmax)-3,np.log10(fmax),4)
    tick_locations = np.hstack((-tick_locations[::-1],0,tick_locations))
    f.colorbar(col2, cax=cax2,ticks=tick_locations,format=ticker.LogFormatterMathtext())
    save_fig(directory+'/v_slice_'+tag+str(i_im)+'_'+format(num,'05'),ps=ps,tight=False)
    #plt.tight_layout()
    #plt.show()

       # cbar = plt.colorbar(im)


        #cbar.set_label(r'$\mathrm{log}_{10}\rho (\mathrm{g}\/\mathrm{cm}^{-3})$')
       # save_fig(directory+'/rho_'+los+'_'+tag+str(i_im)+'_'+format(num,'05'),ps=ps)

##############################################################################
### do a series of zoom
### path : the path
### num the output number
### zoom_v an arrays which contains the zoom
### sinks particle can be overplotted and used to place the zoom
##############################################################################
def make_image_zoom(path,num,zoom_v,i_ims,path_out=None,plot_sinks=False,overwrite=True,center_img=[0.5,0.5,0.5],make_phi=False,deli=',',col_dens_only=False,tag='',gcomp=True,map_size=512,vel_size=20,xyz_sink=3,center_def='None',order='<',ps=False,ind_sink=0,all_sink=True,rotate=None,make_slice=True,col_dens=True,c_axis=None,PSCAL=[],add_B=False):

    ro = pymses.RamsesOutput(path,num,order=order)
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    dx = 0.5**lmax
    print 'len_AU', len_AU
    vars = ["rho","vel","P"]#,"T"]
    if make_phi: vars = vars + ["phi","g"]
    if len(PSCAL)>0: vars = vars + ["passive"]
    if add_B: vars = vars + ["Br","Bl"]
    if not make_phi: amr = ro.amr_source(vars)
    else: amr = ro.amr_source(vars, grav_compat=gcomp)

    if center_def=='sink': 
        sinks = read_sink_cvs(num,path,deli=deli)
        if sinks.shape[0]>=1:
            Msink = np.sum(sinks[:,1])
            center_img = sinks[np.argmax(sinks[:,1]),xyz_sink:xyz_sink+3]/boxlen
        else:
            Msink = 0. 
            center_def='baricenter'
    if center_def=='dmax': center_img = find_center_dmax(ro, c_center=center_img,c_radius=0.005,scl=lmax)
    elif center_def=='baricenter': 
        if center_img[0]==0.5: center_img = find_baricenter(ro, c_center=center_img,c_radius=150./len_AU,scl=lmax)
        center_img = find_baricenter(ro, c_center=center_img,c_radius=60./len_AU,scl=lmax)
    if rotate == 'AM': c_axis, c_vel = find_axis(ro, c_center=center_img,c_radius=20./len_AU,scl=lmax)
    elif rotate is None and c_axis is not None: c_vel = np.array([0.,0.,0.])
    else: c_axis = None

    if make_slice: func_M_R = cal_M_R(ro, c_center=center_img,c_radius=np.max(zoom_v)*np.sqrt(3),scl=lmax)
    
    mass_sink = None; pos_sink = None
    if plot_sinks: #read sink file
        sinks = read_sink_cvs(num,path,deli=deli)
        if len(sinks)>0:
            mass_sink = sinks[:,1] * mass_sol
            pos_sink = sinks[:,xyz_sink:xyz_sink+3]/boxlen
#            ind_max = np.argmax(sinks[:,1])
            args = np.argsort(mass_sink)[::-1]
            ind_sink = np.min((ind_sink,len(sinks)))
            ind_max = args[ind_sink] #position of the most massive sink
            if center_def=='sink':
                center_img=pos_sink[ind_max,:]
                tag= 'ns'+str(ind_sink)+'_'
    #for i_im in range(len(zoom_v))[1:2]:
    for i_im, rad in zip(i_ims, zoom_v):
#        rad = zoom_v[i_im]
        center = center_img.copy()
        center[center<rad]=rad; center[center>(1-rad)]=1-rad #re-center big box
        print 'center ', center, 'rad ', rad

        if make_phi: make_image_phi(amr,ro,center,rad,num,path,x_v=x_v,y_v=y_v,z_v=z_v,i_im=i_im,force=force,path_out=path_out)
        if col_dens: make_image(amr,ro,center,rad,num,path,pos_sink=pos_sink,mass_sink=mass_sink,i_im=i_im,overwrite=overwrite,path_out=path_out,col_dens_only=col_dens_only,tag=tag,map_size=map_size,vel_size=vel_size,ps=ps,all_sink=all_sink,ax1=c_axis,PSCAL=PSCAL,add_B=add_B)
        if make_slice: make_image_slice(amr,ro,center,rad,num,path,func_M_R,i_im=i_im,overwrite=overwrite,path_out=path_out,col_dens_only=False,map_size=50,tag='',ps=False,all_sink=True,ax1=c_axis,ax2=None,n_slice=3,step_slice=10,c_vel=c_vel,PSCAL=PSCAL)
    return center
##############################################################################
def velocity_ana(path,num,zoom_v,path_out=None,overwrite=True,center_img=[0.5,0.5,0.5],tag='',gcomp=True,center_def='baricenter',order='<',ps=False,rotate='AM',radius=100, nbins=30, i_im=0):

    ro = pymses.RamsesOutput(path,num,order=order)
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    #if not make_phi: 
    amr = ro.amr_source(["rho","vel","P","Bl","Br"])
    #else:  amr = ro.amr_source(["rho","vel","P","phi","g"], grav_compat=gcomp)

    if center_def=='dmax': center_img = find_center_dmax(ro, c_center=center_img,c_radius=0.005,scl=lmax)
    elif center_def=='baricenter': center_img = find_baricenter(ro, c_center=center_img,c_radius=250./len_AU,scl=lmax-1)
    if rotate == 'AM': c_axis, c_vel = find_axis(ro, c_center=center_img,c_radius=50./len_AU,scl=lmax)
    else: c_axis = None

    radius = radius/len_AU
    func_M_R = cal_M_R(ro, c_center=center_img,c_radius=radius*1.5,scl=lmax)
    
    #overwrite = True
    #nbins = 30
    vel_profiles(amr,ro,center_img,radius,num,path,func_M_R,i_im=i_im,overwrite=overwrite,path_out=path_out,nbins=nbins,tag='',ps=False,ax1=c_axis,c_vel=c_vel,scl=20)

###########################################################################

##############################################################################


