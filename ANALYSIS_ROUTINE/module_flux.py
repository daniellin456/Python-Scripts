#! /usr/bin/env python

## The module flux contains analyses of disc flux properties

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
from functools import partial
#import pdb
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.signal import savgol_filter, argrelextrema
#import multiprocessing as mp
##### module dependences ##################################################################
from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis, normalisation, read_sink_cvs, norm_units
from module_visualization import cal_rot_mat, save_fig
from module_analysis import cylinder_flux, save_pickle, read_pickle
#from module_core import normalisation
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
Ms_sink = 1.9889e33

#############################################################################################
#measure flux around in particle in across a cubic surface
def flux_cube(time,ioutput,path_in,path_out):
    ind = path[::-1].find('/')
    if ind>0: titlepath = path_in[len(path_in)-ind:]
    else: titlepath = path_in
    outputs = search_ro(path_in)
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
 
#############################################################################################
#visualize rho, rhovr, rhovz, vz/vr in four quadrants
def plot_grids(time,ioutput,grid_m,grid_V,grid_mvr,grid_mvz,rr,zz,len_AU,path_out):
    #images = [grid_m/grid_V, grid_mvr/grid_V, grid_mvphi/grid_V, grid_mvz/grid_V, grid_mvz/grid_mvr]
    fig, ax = plt.subplots(2,2,figsize=(8,2.7))
    ax[0,0].imshow(np.log10(np.flipud(grid_m/grid_V)).T,origin='lower',extent=np.array([rr[-1],rr[0],zz[0],zz[-1]])*len_AU,aspect='equal')
    ax[0,1].imshow(np.log10(-grid_mvr/grid_V).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU,aspect='equal',cmap='hot')
    ax[1,1].imshow(np.log10(-grid_mvz/grid_V).T,origin='upper',extent=np.array([rr[0],rr[-1],zz[-1],zz[0]])*len_AU,aspect='equal',cmap='hot')
    ax[1,0].imshow(np.log10(np.flipud(grid_mvz/grid_mvr)).T,origin='upper',extent=np.array([rr[-1],rr[0],zz[-1],zz[0]])*len_AU,cmap='seismic',aspect='equal',vmin=-2, vmax = 2)
    ax[0,0].set_xticklabels([])
    ax[0,1].set_xticklabels([])
    ax[0,1].set_yticklabels([])
    ax[1,1].set_yticklabels([])
    fs = 16
    ax[0,0].set_ylabel(r'$\rho$',fontsize=fs)
    ax[1,0].set_ylabel(r'$v_z/v_r$',fontsize=fs)
    ax[0,1].yaxis.set_label_position("right")
    ax[0,1].set_ylabel(r'$\rho v_r$',fontsize=fs)
    ax[1,1].yaxis.set_label_position("right")
    ax[1,1].set_ylabel(r'$\rho v_z$',fontsize=fs)
    fig.subplots_adjust(wspace=0.,hspace=0.)
    fig.suptitle('Time = %i kyr'%(time*1.e3), fontsize=20)
    plt.savefig(path_out+'/fluxrz_'+str(ioutput).zfill(5)+'.pdf')
#############################################################################################
# find disc surface as function of radius
def find_disc_surf(rr,zz,dx,grid_rho,thres_rho,dfact,r_max,iz_max,ir_max,len_AU,overplot=False,contours=False):
    nr = len(rr)
    z_r_disc = np.zeros(nr)
    dsurf = np.zeros(nr)
    if contours: rr_disc = []; zz_disc = []
    ind_z = 0
    for ir in range(nr-1):
        r = rr[ir]
        if r > r_max*1.1 or grid_rho[ir,0]<thres_rho*dfact: break # break if radius too large or density too low
        if z_r_disc[ir] == 0.:
            dd = grid_rho[ir,:iz_max]
            z_loc = np.linspace(dx,zz[iz_max]-dx*1.5,41)
            d_z = interp1d(zz[:iz_max]+dx*0.5,np.log(dd),kind='cubic')
            d2d_dz2 = lambda z: derivative(d_z, z, dx*0.5,n=2)
            d2logd_dlogz2 = d2d_dz2(z_loc)
            ind_z = np.argmax(d2logd_dlogz2)
            #ind_z = argrelextrema(d2logd_dlogz2,np.greater)[0][0]
            z_r_disc[ir] = z_loc[ind_z] #- dx
            ds = np.exp(d_z(z_r_disc[ir]))
            dsurf[ir]=ds
            if ir > 10:
                if z_r_disc[ir]/z_r_disc[ir-1]<0.3: break
#        plt.scatter(z_loc[ind_z],np.exp(d_z(z_loc[ind_z])))
#        plt.semilogy(zz[1:iz_max-1],d2d_dz2(zz[1:iz_max-1]))
#        plt.semilogy(zz[1:iz_max-1],-d2d_dz2(zz[1:iz_max-1]),':')
#        plt.semiilogy(zz[:iz_max],dd)
        if ir > ir_max:  ## find iso-density surface for r>r_max
             if  grid_rho[ir+1,0] < ds: break
             iz_max = min(iz_max+4, len(zz)-1)
             z_d = interp1d(grid_rho[ir+1,:iz_max],range(int(iz_max)))(ds)
             z_d_i = int(np.ceil(z_d))
             dsurf[ir] = ds
             if grid_rho[ir+1,z_d_i-1] < grid_rho[ir,z_d_i-1] and ir > 16:
                 z_r_disc[ir+1] = (z_d+0.5)*dx
             else: iz_max = 10 #min(max(z_d_i+3,4), len(zz)-1)
    if overplot:
        f=overplot[0]; ax=overplot[1]
        plt.figure(f.number)
        ax[0,0].plot(rr*len_AU,z_r_disc*len_AU)
    if contours:
        floc=plt.figure()
        plt.plot(rr,z_r_disc)
        rr_disc.append(rr[z_r_disc>0]+dx*0.5); zz_disc.append(z_r_disc[z_r_disc>0])
        ds = np.mean(dsurf[z_r_disc>0]) # iterate for better disc surface
        iz_max = np.ceil(max(z_r_disc)/dx+0.5)
        facts_d = np.linspace(0.8,1.5,3)
        for fact in facts_d:
          iz_max_loc = np.floor(1.5/len_AU/dx)
          for ir in range(nr):
            r = rr[ir]
            if r > r_max*1.1 or grid_rho[ir,0]<thres_rho*dfact: break # break if radius too large or density too low
            if np.sum(grid_rho[ir,:iz_max] > ds*fact) :
              iz_max_loc = np.argmax(np.append(np.diff(grid_rho[ir,iz_max_loc-np.ceil(1./len_AU/dx):iz_max]),1)>0)+iz_max_loc
              dd = grid_rho[ir,:iz_max_loc]
              z_d = interp1d(np.log(grid_rho[ir,:iz_max_loc]),range(int(iz_max_loc)),kind='linear')(np.log(ds*fact))
              z_r_disc[ir] = (z_d+0.5)*dx
          plt.figure(floc.number)
          plt.plot(rr,z_r_disc,label=str(fact))
          if overplot:
              plt.figure(f.number)
              ax[0,0].plot(rr*len_AU,z_r_disc*len_AU)
              ax[1,1].plot(rr*len_AU,z_r_disc*len_AU)
          rr_disc.append(rr[z_r_disc>0]+dx*0.5)
          zz_disc.append(z_r_disc[z_r_disc>0])
          rr_disc.append(rr[z_r_disc>0]+dx*0.5)
          zz_disc.append(np.ones_like(zz_disc[-1])*max(zz_disc[-1]))
        plt.figure(floc.number); plt.legend()
    else: rr_disc = rr[z_r_disc>0]+dx*0.5; zz_disc = z_r_disc[z_r_disc>0]
    return rr_disc, zz_disc
#############################################################################################
# calculate [r, z] grids of various quantities
def calc_grids(cell_disc,posi_disc,c_vel,vec_z,len_AU,thres_rho_list,dV_disc,dS_disc,rr,zz,plot=False):
    vel_disc = cell_disc["vel"]-c_vel[np.newaxis,:]
    dm_disc = cell_disc["rho"]*dV_disc
    rhoS_disc = cell_disc["rho"]*dS_disc
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
    inds_dens = np.where(cell_disc["rho"] > thres_rho_list*1)
    grid_x = lambda x: np.histogramdd((distc_disc, abs(z_disc),), bins=[rr,zz], weights=x)[0]
    grid_m = grid_x(dm_disc)
    grid_V = grid_x(dV_disc)
    grid_S = grid_x(dS_disc)
    grid_rho = grid_m/grid_V
#        grid_rho = np.ma.masked_where(np.isnan(grid_rho), grid_rho)
    grid_mvr = grid_x(dm_disc*vel_r)
    grid_mvphi = grid_x(dm_disc*vel_phi)
    grid_mvz = grid_x(dm_disc*vel_z*up_low)
    grid_rhoSvr =  grid_x(rhoS_disc*vel_r)
    grid_rhoSvz =  grid_x(rhoS_disc*vel_z*up_low)
    j_ang = np.cross(posi_disc,vel_disc)
#    grid_Jzvr = grid_x(dm_disc*j_ang[:,2]*vel_r)
#    grid_Jzvz = grid_x(dm_disc*j_ang[:,2]*vel_z*up_low)
    grid_Jz = grid_x(dm_disc*j_ang[:,2])
    r_max = max(distc_disc[inds_dens])
    if plot:
        f,ax = plt.subplots(3,2)
        ax[0,0].imshow(np.log10(grid_rho).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        ax[1,0].imshow(np.log10(-grid_mvr/grid_V).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        ax[1,1].imshow(np.log10(-grid_mvz/grid_V).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        #ax[1,0].imshow(np.log10(-grid_rhoSvr/grid_S).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        #ax[1,1].imshow(np.log10(-grid_rhoSvz/grid_S).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        ax[2,0].imshow(np.log10(-grid_mvr/grid_m).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        ax[2,1].imshow(np.log10(-grid_mvz/grid_m).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        ax[0,1].imshow(np.log10(grid_mvz/grid_mvr).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU,cmap='seismic',vmin=-2,vmax=2)
        #ax[3,0].imshow(np.log10(-grid_Jzvr/grid_V).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        #ax[3,1].imshow(np.log10(-grid_Jzvz/grid_V).T,origin='lower',extent=np.array([rr[0],rr[-1],zz[0],zz[-1]])*len_AU)
        return grid_m, grid_V, grid_S, grid_rho, grid_mvr, grid_mvphi, grid_mvz, grid_rhoSvr, grid_rhoSvz, grid_Jz, r_max, f, ax
    else: return grid_m, grid_V, grid_S, grid_rho, grid_mvr, grid_mvphi, grid_mvz, grid_rhoSvr, grid_rhoSvz, grid_Jz, r_max
#############################################################################################
# reduce the resolution of grids from lmax to lnew
def reduce_reso(grid_m, grid_V, grid_mvz,rr, zz, lmax,lnew=14):
    lrd = lmax
    grid_m_rd = grid_m.copy(); grid_V_rd = grid_V.copy()
    grid_mvz_rd = grid_mvz.copy()
    rr_rd = rr.copy(); zz_rd = zz.copy()
    while lrd > lnew:
        dim = grid_m_rd.shape
        nr_rd = (dim[0])/2; nz_rd = (dim[1])/2
        grid_m_rd = grid_m_rd[:nr_rd*2:2,:nz_rd*2:2] + grid_m_rd[1:nr_rd*2:2,:nz_rd*2:2] + grid_m_rd[:nr_rd*2:2,1:nz_rd*2:2] + grid_m_rd[1:nr_rd*2:2,1:nz_rd*2:2]
        grid_V_rd = grid_V_rd[:nr_rd*2:2,:nz_rd*2:2] + grid_V_rd[1:nr_rd*2:2,:nz_rd*2:2] + grid_V_rd[:nr_rd*2:2,1:nz_rd*2:2] + grid_V_rd[1:nr_rd*2:2,1:nz_rd*2:2]
        grid_mvz_rd = grid_mvz_rd[:nr_rd*2:2,:nz_rd*2:2] + grid_mvz_rd[1:nr_rd*2:2,:nz_rd*2:2] + grid_mvz_rd[:nr_rd*2:2,1:nz_rd*2:2] + grid_mvz_rd[1:nr_rd*2:2,1:nz_rd*2:2]
        lrd -= 1;
        rr_rd = rr_rd[:nr_rd*2+1:2]
        zz_rd = zz_rd[:nz_rd*2+1:2]
    grid_rho_rd = grid_m_rd / grid_V_rd
    dx_fac = 2**(lmax-lrd)
    return grid_rho_rd,grid_mvz_rd, dx_fac, rr_rd, zz_rd
#############################################################################################
# fit disc to a parametric function
def best_fit_disc(rr_disc, zz_disc,fix_re=False):
    func_zr = lambda r,a,b,re,n: a * r**b * (re-r)**n
    bounds = ([0,0.6,rr_disc[-1]-rr_disc[0],0.1],[0.5,2,rr_disc[-1]+rr_disc[0],0.3])
    p0 = (0.3,1,rr_disc[-1],0.2)
    if fix_re:
        zz_disc[-1]=0.
        p0 = (p0[0],p0[1],p0[3])
        re = fix_re
        func_zr = lambda r,a,b,n: a * r**b * (re-r)**n
        popt, pcov = curve_fit(func_zr, rr_disc, zz_disc,p0=p0)
        popt = [popt[0], popt[1], re, popt[2]]
    else:
        popt, pcov = curve_fit(func_zr, rr_disc, zz_disc,p0=p0)#bounds=bounds)
    print 'fitting parameters', popt
    return popt
#############################################################################################
# measure mass flux onto a given surface 
def measure_surface_flux(grid_rhovr,grid_rhovz,dx,nr,re,func_zr, func_tan_r,rr_disc,zz,zz_disc,zz_smooth2,len_AU,Mdot_MsMyr,plot_flux=False,plot_cum = False,contours=False,nsmth=2,racc=4.5,figname=None):
    if plot_flux:
        f2, ax2 = plt.subplots(figsize=(7,3))
        ax2b = ax2.twinx()
        ax2.plot(rr_disc*len_AU,zz_disc*len_AU,'gray')
        ax2.plot(rr_disc*len_AU,zz_smooth2*len_AU,'k',lw=2,label='Disk height (AU)')
    Mdot_dict = {'Mdotr':[],'Mdotz':[],'Mdot':[],'mdotr':[],'mdotz':[]}
    grid_rhovr_pad = np.hstack((grid_rhovr[:,0,np.newaxis],grid_rhovr))
    grid_rhovz_pad = np.hstack((grid_rhovz[:,0,np.newaxis],grid_rhovz))
    zz_pad = np.append(0.,zz[:-1]+dx*0.5)
    if contours:
        ncon = len(rrs_disc)
        iter = ncon
        nplot = 1
        scale = ncon
    else:
        nz = 5; iter = nz+1
        nplot = 1
        fact_max = (max(zz)-dx*0.5)/max(zz_smooth2)
        facts_z = np.linspace(0.8,min(2.,fact_max),nz)
        facts_z = np.sort(np.append(facts_z,[1.]))
        scale = facts_z
    for ii in range(iter):
        if contours:
            fact_z = 0.; iz = icon = ii
            if icon == 0: zz_disc = zz_smooth
            else:
                rr_disc = rrs_disc[icon]
                nr =  len(rr_disc)
                func_zr = interp1d(rrs_disc[icon],zzs_disc[icon],kind='linear')
                print rrs_disc[icon],zzs_disc[icon]
                rr_disc = np.linspace(0,rr_disc[-1]+dx*0.5, nr*2+1)
                zz_disc = np.zeros_like(rr_disc)
                zz_disc[1:-1] = func_zr(rr_disc[1:-1])
                zz_disc[0] = zz_disc[1]; zz_disc[-1]=0.
        else:
            iz = ii
            fact_z = facts_z[iz]
            zz_disc = zz_smooth2 * fact_z
#          zz_disc = np.maximum(zz_disc,np.sqrt(np.maximum(20.25*dx**2-rr_disc**2,np.zeros_like(rr_disc))))
        zz_disc0 = zz_disc.copy()
        tan_r0 = func_tan_r(rr_disc)
        tan_r = tan_r0[1::2]
        #zz_disc = ((racc*dx)**nsmth + zz_disc**nsmth)**(1./nsmth) ## limit to > sink r_acc
        #tan_r = func_tan_r(rr_disc[1::2]) * (zz_disc0[1::2]/zz_disc[1::2])**(nsmth-1) * fact_z
        negtan = np.where(tan_r0<0)[0]
        if len(negtan)>0: ide = negtan[0]
        else: ide = len(zz_disc)
        zz_disc[:ide] = np.maximum(racc*dx,zz_disc[:ide])
        func_zr2 = interp1d(rr_disc,zz_disc,kind='linear')
        tan_r[:ide/2][zz_disc[1:ide:2]==racc*dx] = 0.
        mdotr = np.zeros(nr)
        mdotz = np.zeros(nr)
        for ir in range(nr):
            ir2 = ir * 2 + 1
            fac = np.ceil(abs(tan_r[ir]))
            zmin = 0; iz_max = int(np.ceil(zz_disc[ir2]/dx)+2+fac)
            func_mvr = interp1d(zz_pad[zmin:iz_max], grid_rhovr_pad[ir,zmin:iz_max],kind='linear')
            func_mvz = interp1d(zz_pad[zmin:iz_max], grid_rhovz_pad[ir,zmin:iz_max],kind='linear')
            if abs(tan_r[ir]) <= 1.:
                surf = 2. * np.pi * rr_disc[ir2] * dx
                mdotr[ir] = func_mvr(zz_disc[ir2]) * tan_r[ir] * surf
                mdotz[ir] = - func_mvz(zz_disc[ir2]) * surf
            else:
                nbin = np.ceil(abs(tan_r[ir]))
                sub_tan = np.ones(int(nbin)) * nbin
                while np.sum(abs(sub_tan) >= nbin):
                    nbin *= abs(sub_tan[-1])/nbin
                    if ir == nr-1:
                        dxf = rr_disc[ir2+1]-rr_disc[ir2-1]
                        sub_r = np.linspace(rr_disc[ir2-1],rr_disc[ir2+1],nbin+1)[:-1] + 0.5*dxf/nbin
                        sub_surf = 2. * np.pi * sub_r * dxf/ nbin
                    else:
                        sub_r = np.linspace(rr_disc[ir2-1],rr_disc[ir2+1],nbin+1)[:-1] + 0.5*dx/nbin
                        sub_surf = 2. * np.pi * sub_r * dx / nbin
                    sub_z0 = func_zr(sub_r) * fact_z
                    #sub_z = ((racc*dx)**nsmth + sub_z0**nsmth)**(1./nsmth)
                    #sub_tan = func_tan_r(sub_r)*(sub_z0/sub_z)**(nsmth-1) * fact_z
                    sub_z = sub_z0; sub_tan = func_tan_r(sub_r)
                mdotr[ir] = np.sum( func_mvr(sub_z) * sub_tan * sub_surf )
                mdotz[ir] = np.sum( - func_mvz(sub_z) * sub_surf )
## add the last bin that have only radial contribution
        if abs(tan_r[ir]) <= 1.: zlim = zz_disc[ir2] + 0.5*dx * tan_r[ir]
        else: zlim = sub_z[-1] + 0.5*dxf/nbin * sub_tan[-1]
        nzlim = int(np.ceil(zlim/dx))
#        print ir, racc, zlim/dx,np.arange(0.5*dx,zlim,dx)/dx
        mdotr_v  = - grid_rhovr_pad[ir,1:nzlim+1]
        mdotr_v[-1] = mdotr_v[-1] * (zlim/dx-nzlim+1) ## corect the last bin
        mdotr_v = np.sum(mdotr_v)*2*np.pi*re*dx
        mdotr[ir] += mdotr_v

        Mdotr_tot = np.sum(mdotr)#*Mdot_MsMyr
        Mdotz_tot = np.sum(mdotz)#*Mdot_MsMyr
        Mdot_tot = Mdotr_tot + Mdotz_tot

        mdotr_c = np.cumsum(mdotr)
        mdotz_c = np.cumsum(mdotz)
        
        Mdot_dict['Mdotr'].append(Mdotr_tot)
        Mdot_dict['Mdotz'].append(Mdotz_tot)
        Mdot_dict['Mdot'].append(Mdot_tot)
        Mdot_dict['mdotr'].append(mdotr)
        Mdot_dict['mdotz'].append(mdotz)

        if  plot_flux:
          if iz%nplot==0 or fact_z == 1.:
            if plot_cum: mdotrp = mdotr_c; mdotzp = mdotz_c
            else: mdotrp = mdotr; mdotzp = mdotz
            if fact_z != 1.:
                ax2.plot(rr_disc*len_AU,zz_disc*len_AU, 'C'+str((iz/nplot)%10),lw=0.5)
                ax2b.plot(rr_disc[1::2]*len_AU,(mdotrp+mdotzp)*Mdot_MsMyr, 'C'+str((iz/nplot)%10), label='%.2f Ms/Myr(%.2f z)'%(Mdot_tot*Mdot_MsMyr,fact_z))
            if fact_z ==1:
                ax2.plot(rr_disc*len_AU,zz_disc*len_AU, 'k',lw=0.5)
                ax2b.plot(rr_disc[1::2]*len_AU,(mdotrp+mdotzp)*Mdot_MsMyr, 'b', label='%.2f Ms/Myr'%(Mdot_tot*Mdot_MsMyr))
                #ax2b.plot(rr_disc[1::2]*len_AU,(mdotrp)*Mdot_MsMyr, 'k--',label='$\dot{M}_r(<R)$') #label='%.1f Ms/Myr(%.2f z)'%(Mdotr_tot,fact_z))
                #ax2b.plot(rr_disc[1::2]*len_AU,(mdotzp)*Mdot_MsMyr, 'k:',label='$\dot{M}_z(<R)$') #label='%.1f Ms/Myr(%.2f z)'%(Mdotz_tot,fact_z))
          ax2.set_xlabel(r'$R$ (AU)',fontsize=15)
          ax2.set_xlim([0,40])
          ax2.set_ylim([0,10])
          ax2.set_ylabel(r'$H$ (AU)',fontsize=15)
          ax2b.set_yscale('log')
          ax2b.set_ylim([1.e-3,1.e1])
          ax2b.set_ylabel(r'$S$ (Kg s$^{-1}$ m$^{-2})$',fontsize=15)
        #  plt.savefig(figname)
    if plot_flux: return Mdot_dict, scale, f2,ax2,ax2b
    else: return Mdot_dict, ncon, facts_z
#############################################################################################
def write_disc_history(path_out,ioutput,time,Msink,grid_m,grid_rhovr,grid_rhovz,grid_Jz,dx,nr,re,func_zr,func_tan_r,rr_disc,zz,zz_disc,len_AU,mass_sol,Mdot_MsMyr,nsmth=2,racc=4.5,overwrite=True):
    filename = (path_out+'/disc_history.npz') 
    if not os.path.exists(filename): 
        outputs=[]; times=[]; Msinks=[]; Mdiscs=[]; Mdots=[]; rs=[]; zs=[]
    else: 
        [outputs, times, Msinks, Mdiscs, Mdots, rs, zs] = read_pickle(filename,7)
    if ioutput in outputs and not overwrite: return
    else:
        if ioutput in outputs:
            iout = np.where(outputs==ioutput)[0]
            print iout, outputs,ioutput
        else:
            outputs.append(ioutput)
            times.append(time)
            Msinks.append(Msink)
            iout = len(outputs)
    print 'disc_history', iout, ioutput
    print outputs

    grid_rhovr_pad = np.hstack((grid_rhovr[:,0,np.newaxis],grid_rhovr))
    grid_rhovz_pad = np.hstack((grid_rhovz[:,0,np.newaxis],grid_rhovz))
    zz_pad = np.append(0.,zz[:-1]+dx*0.5)

    tan_r0 = func_tan_r(rr_disc)
    tan_r = tan_r0[1::2]
    ide = np.where(tan_r0<0)[0][0]
    zz_disc[:ide] = np.maximum(racc*dx,zz_disc[:ide])
    tan_r[:ide/2][zz_disc[1:ide:2]==racc*dx] = 0.
    mdotr = np.zeros(nr)
    mdotz = np.zeros(nr)
    m_r = np.zeros(nr)
    for ir in range(nr):
        ir2 = ir * 2 + 1
        fac = np.ceil(abs(tan_r[ir]))
        zmin = 0; iz_max = int(np.ceil(zz_disc[ir2]/dx)+2+fac)
        func_mvr = interp1d(zz_pad[zmin:iz_max], grid_rhovr_pad[ir,zmin:iz_max],kind='linear')
        func_mvz = interp1d(zz_pad[zmin:iz_max], grid_rhovz_pad[ir,zmin:iz_max],kind='linear')
        func_m = interp1d(zz[1:iz_max+1], np.cumsum(grid_m[ir,0:iz_max]),kind='linear')
        if abs(tan_r[ir]) <= 1.:
            surf = 2. * np.pi * rr_disc[ir2] * dx
            mdotr[ir] = func_mvr(zz_disc[ir2]) * tan_r[ir] * surf
            mdotz[ir] = - func_mvz(zz_disc[ir2]) * surf
            m_r[ir] = func_m(zz_disc[ir2])
        else:
            nbin = np.ceil(abs(tan_r[ir]))
            sub_tan = np.ones(int(nbin)) * nbin
            while np.sum(abs(sub_tan) >= nbin):
                nbin *= abs(sub_tan[-1])/nbin
                if ir == nr-1:
                    dxf = rr_disc[ir2+1]-rr_disc[ir2-1]
                    sub_r = np.linspace(rr_disc[ir2-1],rr_disc[ir2+1],nbin+1)[:-1] + 0.5*dxf/nbin
                    sub_surf = 2. * np.pi * sub_r * dxf/ nbin
                else:
                    sub_r = np.linspace(rr_disc[ir2-1],rr_disc[ir2+1],nbin+1)[:-1] + 0.5*dx/nbin
                    sub_surf = 2. * np.pi * sub_r * dx / nbin
                sub_z = func_zr(sub_r) 
                sub_tan = func_tan_r(sub_r)
            mdotr[ir] = np.sum( func_mvr(sub_z) * sub_tan * sub_surf )
            mdotz[ir] = np.sum( - func_mvz(sub_z) * sub_surf )
            m_r[ir] = np.mean(func_m(sub_z))
## add the last bin that have only radial contribution
    if abs(tan_r[ir]) <= 1.: zlim = zz_disc[ir2] - 0.5*dx * tan_r[ir]
    else: zlim = sub_z[-1] - 0.5*dxf/nbin * sub_tan[-1]
    nzlim = int(np.ceil(zlim/dx))
    mdotr_v  = - grid_rhovr_pad[ir,1:nzlim+1]
    mdotr_v[-1] = mdotr_v[-1] * (zlim/dx-nzlim+1) ## corect the last bin
    mdotr_v = np.sum(mdotr_v)*2*np.pi*re*dx
    mdotr[ir] += mdotr_v

    if iout == len(outputs):
        Mdiscs.append(m_r*mass_sol)
        Mdots.append((mdotr+mdotz)*Mdot_MsMyr)
        rs.append(rr_disc*len_AU)
        zs.append(zz_disc*len_AU)
    else:
        Mdiscs[iout] = m_r*mass_sol
        Mdots[iout] = (mdotr+mdotz)*Mdot_MsMyr
        rs[iout] = rr_disc*len_AU
        zs[iout] = zz_disc*len_AU

    print outputs
    save_pickle(filename,[outputs, times, Msinks, Mdiscs, Mdots, rs, zs])
#############################################################################################
def plot_disc_history(path_out):
    filename = (path_out+'/disc_history.npz')
    if not os.path.exists(filename): return
    #else: disc_history = np.load(filename)
    else: [outputs, times, Msinks, Mdiscs, Mdots, rs, zs] = read_pickle(filename,7)
    order = np.argsort(outputs)
    times = [times[i] for i in order]
    Msinks = [Msinks[i] for i in order]
    Mdiscs = [Mdiscs[i] for i in order]
    Md = [np.sum(M) for M in Mdiscs]
    Mdots = [Mdots[i] for i in order]
    Macc = [np.sum(M) for M in Mdots]
    rs = [rs[i] for i in order]
    zs = [zs[i] for i in order]
    func_Ms_t = interp1d(times,Msinks)
    dt = np.diff(times)
    print dt
    dtmin = np.minimum(dt[1:],dt[:-1])/5.
    #dt = 0.5*(dt[1:]+dt[:-1])
    Mdot_sink = derivative(func_Ms_t,times[1:-1],dx=dtmin)
    nout = len(times)
#    print disc_history['outputs']
#    print disc_history['times']
    print 'Msinks', Msinks
    print 'Mdiscs', Mdiscs
    print 'rs', rs
    print 'zs', zs
    print order
#    print disc_history['Msinks'].shape
#    print disc_history['Mdiscs'].shape
#    print disc_history['rs'].shape
#    print disc_history['zs'].shape
    fig_M = plt.figure()
    plt.plot(times,Msinks,label='Msink')
    plt.plot(times,Md,label='Mdisc')
    plt.plot(times,np.asarray(Msinks)+np.asarray(Md),label='Mtot')
    plt.plot(times[1:-1],Mdot_sink,label='Mdot_sink')
    plt.plot(times,Macc,label='Mdot_disc')
    plt.yscale('log')
    plt.legend()
    plt.xlabel('time (Myr)',fontsize=20)
    plt.ylabel(r'$M$ (Ms), $M_{dot}$ (Ms/Myr)',fontsize=20)
    plt.savefig(path_out+'/history_mass.pdf')

    fig_surf = plt.figure()
    for iout in range(nout):
        plt.plot(rs[iout],zs[iout],label='%i kyr, %.2f Ms'%(times[iout]*1.e3, Msinks[iout]))
    plt.legend()
    plt.xlabel(r'$R$ (AU)',fontsize=20)
    plt.ylabel(r'$H$ (AU)',fontsize=20)
    plt.savefig(path_out+'/history_height.pdf')

    fig_Mdot = plt.figure()
    for iout in range(nout):
        plt.plot(rs[iout][1::2],np.asarray(Mdots[iout]),label='%i kyr, %.2f Ms'%(times[iout]*1.e3, Msinks[iout]))
    plt.yscale('log')
    plt.legend()
    plt.xlabel(r'$R$ (AU)',fontsize=20)
    plt.ylabel('flux (kg/m$^2$)',fontsize=20)
    plt.savefig(path_out+'/history_flux.pdf')
    plt.show()

#############################################################################################
def detail_disc_flux(path,path_out=None,overwrite=True,order='<',scl=20,contours=False,Msink_sams=None,time_sams=None):
    if path_out == None: path_out = path
    units = norm_units(path)
    boxlen = units["boxlen"]; lmax = units["lmax"]
    len_AU = units["len_AU"]
    mass_sol = units["mass_sol"]
    Mdot_MsMyr = units["mass_sol"]/units["time_Myr"]/boxlen
    unit_flux = units["dens_gcc"] * units["vel_ms"] * 1.e3  # SI
    unit_jflux = unit_flux * units["vel_ms"] * units["len_m"] # SI 
    Ms_sink = 1.9889e33
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
#np.savez(path_out+'/disc_basic_parameters', ioutput_list=ioutput_list, time_list=time_list, Msink_list=Msink_list, Mdisc_list=Mdisc_list, thres_rho_list=thres_rho_list, center_sink_list=center_sink_list, center_disc_list=center_disc_list, cvel_sink_list=cvel_sink_list, cvel_disc_list=cvel_disc_list, ax_csink_list=ax_csink_list, ax_ctot_list=ax_ctot_list, r_csink_list=r_csink_list, r_ctot_list=r_ctot_list)
    discload = np.load(path_out+'/disc_basic_parameters.npz')
    outputs = discload['ioutput_list']
    Msink=np.asarray(discload['Msink_list'])
    thres_rho_list = discload['thres_rho_list']
    center_sink_list=np.asarray(discload['center_sink_list'])
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])
    center_disc_list=np.asarray(discload['center_disc_list'])
    cvel_disc_list=np.asarray(discload['cvel_disc_list'])
#    Mdisc=np.asarray(discload['Mdisc_list'])
    ax_csink = np.asarray(discload['ax_csink_list'])
    ax_ctot = np.asarray(discload['ax_ctot_list'])
    r_csink = np.asarray(discload['r_csink_list'])
    time_list = np.asarray(discload['time_list'])
    nout = len(outputs)
    for i in range(nout):
        print i, outputs[i], time_list[i]*1.e3
    if time_sams: isams = [np.argmin(abs(time_list-time_sam)) for time_sam in time_sams]
    elif Msink_sams: isams = [np.argmin(abs(Msink-Msink_sam)) for Msink_sam in Msink_sams]
    else: isams = range(nout/4,nout,nout/4)
    print isams, [outputs[isam] for isam in isams], [time_list[isam]*1.e3 for isam in isams], [Msink[isam] for isam in isams]
    r_sel = 45 / len_AU
    z_sel = 10 / len_AU
    dx = 1./2**lmax
    rr = np.arange(0,r_sel,dx)
    zz = np.arange(0,z_sel,dx)
    plot_im = True
    plot_flux = True
    for iout in isams:
        ioutput=outputs[iout]
        c_center = center_sink_list[iout]
        c_vel = cvel_sink_list[iout]
        vec_z = ax_csink[iout] 
        r_int = r_csink[iout]
        if np.sum(c_center)==0.:
            c_center = center_disc_list[iout]
            c_vel = cvel_disc_list[iout]
            vec_z = ax_ctot[iout] 
        print c_center, vec_z
        filecell = path_out+'/selected_cells_'+str(ioutput).zfill(5)+'.npz'
        if os.path.exists(filecell):
            cell_disc = np.load(filecell)
            dV_disc = cell_disc["dx"]**3 
            dS_disc = cell_disc["dx"]**2
            posi_disc = cell_disc["x"]-c_center[np.newaxis,:]
        else:
            ro = pymses.RamsesOutput(path,ioutput,order=order)
            amr = ro.amr_source(["rho","vel"])
            disc_region = Cylinder(c_center,vec_z,r_sel,z_sel*2)
            amr_disc = RegionFilter(disc_region,amr)
            cell_disc = CellsToPoints(amr_disc, smallest_cell_level=scl).flatten()
            np.savez(filecell,rho=cell_disc["rho"],vel=cell_disc["vel"],x=cell_disc.points,dx=cell_disc.get_sizes())
            dV_disc = cell_disc.get_sizes()**3 
            dS_disc = cell_disc.get_sizes()**2
            posi_disc = cell_disc.points-c_center[np.newaxis,:]
## calculate grids of variables
        grid_m, grid_V, grid_S, grid_rho, grid_mvr, grid_mvphi, grid_mvz, grid_rhoSvr, grid_rhoSvz, grid_Jz, r_max, f, ax = calc_grids(cell_disc,posi_disc,c_vel,vec_z,len_AU,thres_rho_list[iout],dV_disc,dS_disc,rr,zz,plot=plot_im)
        #plt.show()
        #continue
        if plot_im: overplot=(f,ax)
        else: overplot=False
        #grid_rhovr = grid_rhoSvr/grid_S
        #grid_rhovz = grid_rhoSvz/grid_S
        grid_rhovr = grid_mvr/grid_V
        grid_rhovz = grid_mvz/grid_V
        grid_jz = grid_Jz/grid_m
         
## reduce resolution to 1 AU for disc surface identification
        lnew=lmax
        grid_rho_rd, grid_mvz_rd,dx_fac, rr_rd, zz_rd = reduce_reso(grid_m, grid_V,grid_mvz, rr, zz, lmax,lnew=lnew)
        print 'dx_fac = ', dx_fac, grid_rho_rd.shape
        n_rd = len(rr_rd)
        #nskip = int(21/len_AU/(dx*dx_fac)) # impose minimum r_disc (AU)
        nskip = int(r_int*0.8/(dx*dx_fac))
        posvz = np.where(np.sum(grid_mvz_rd[nskip:,:4],axis=1)>0)[0]
        if len(posvz)>0: indvz = posvz[0] + nskip
        else: indvz = n_rd
        #indvz = n_rd
        print indvz, n_rd
        print grid_mvz_rd[:,0]
        grid_rho_rd = grid_rho_rd[:indvz,:]; rr_rd = rr_rd[:indvz]

## identify disc surface
        dfact = 0.5 ## factor of threshold density to find isodensity contour
        iz_max = int(np.ceil(5/len_AU/dx)) ## maximal height in AU
        ir_max = 200# np.floor(r_max/dx)-8
        if contours:
            rrs_disc, zzs_disc = find_disc_surf(rr_rd,zz_rd,dx*dx_fac,grid_rho_rd,thres_rho_list[iout],dfact,r_max,int(iz_max/dx_fac),int(ir_max/dx_fac),len_AU,overplot=overplot,contours=True)
            rr_disc = rrs_disc[0]; zz_disc = zzs_disc[0]
        else: rr_disc, zz_disc = find_disc_surf(rr_rd,zz_rd,dx*dx_fac,grid_rho_rd,thres_rho_list[iout],dfact,r_max,int(iz_max/dx_fac),int(ir_max/dx_fac),len_AU,overplot=overplot)
        nr =  len(rr_disc)
        #rr_disc = (rr_rd[:-1] + rr_rd[1:])*0.5 
        zz_cen = (zz_rd[:-1] + zz_rd[1:])*0.5
        print grid_rho_rd.shape, rr_rd.shape, zz_rd.shape
        zh_disc = np.sqrt( np.sum(grid_rho_rd[:nr,:] * zz_cen[np.newaxis,:]**2,axis=1)/np.sum(grid_rho_rd[:nr,:],axis=1) )
        if plot_im:
            for aa in ax.flatten(): aa.plot(rr_disc*len_AU, zz_disc*len_AU,'r')
            for aa in ax.flatten(): aa.plot(rr_disc*len_AU, zh_disc*len_AU,'w')
        #plt.show() 
## fit disc surface to empirical function
        if indvz < n_rd: [a,b,re,n] = best_fit_disc(rr_disc, zz_disc,fix_re=rr_disc[-1]+dx*dx_fac*0.5)
        else: [a,b,re,n] = best_fit_disc(rr_disc, zz_disc) 
## interpolate disc surface and restore to originate resolution + 1 level to have bin-centered points
        func_zr = interp1d(rr_disc,zz_disc,kind='linear')
        rr_disc = np.linspace(0,rr_disc[-1]+dx*dx_fac*0.5, nr*2*dx_fac+1)
        zz_disc = np.zeros(nr*2*dx_fac+1)
        zz_disc[dx_fac:-dx_fac] = func_zr(rr_disc[dx_fac:-dx_fac])
        zz_disc[:dx_fac] = zz_disc[1]; zz_disc[-dx_fac:]=0.
#        zz_smooth = savgol_filter(zz_disc,17,3)
#        zz_smooth2 = savgol_filter(zz_smooth,13,3)
#        zz_smooth2[-1] = 0.

## reconstruct empirical disc surface
        func_zr = lambda r: a * r**b * (re-r)**n
        func_tan_r = lambda r: a * r**b * (re-r)**n * (b/r - n/(re-r))
        if indvz==n_rd: rr_disc[-1] = re
        zz_smooth2 = func_zr(rr_disc) ## empirical disc contour
        print zz_smooth2
        print rr_disc/dx
        if plot_im: ax[0,0].plot(rr_disc*len_AU,zz_smooth2*len_AU,'y')
## measure flux on a series of surfaces
        plot_cum = False
        #figname =  path_out+'/surface_flux_'+str(ioutput).zfill(5)+'.pdf'
        Mdot_dict, scale, f2, ax2, ax2b = measure_surface_flux(grid_rhovr,grid_rhovz,dx,nr*dx_fac,re,func_zr, func_tan_r,rr_disc,zz,zz_disc,zz_smooth2,len_AU,Mdot_MsMyr,plot_flux=plot_flux,plot_cum = plot_cum,contours=contours,nsmth=2,racc=4.5)

        Mdot_tot = Mdot_dict['Mdot']
#        ax2b.plot(rr_disc[1::2]*len_AU,-grid_rhovr[:nr,int(6/len_AU/dx)]*2*np.pi*rr_disc[1::2]*dx*Mdot_MsMyr,'grey',lw=4)
#        ax2b.plot(rr_disc[1::2]*len_AU,-grid_rhovr[:nr,int(4/len_AU/dx)]*2*np.pi*rr_disc[1::2]*dx*Mdot_MsMyr,'grey',lw=3)

        write_disc_history(path_out,ioutput,time_list[iout],Msink[iout],grid_m,grid_rhovr,grid_rhovz,grid_Jz,dx,nr,re,func_zr,func_tan_r,rr_disc,zz,zz_smooth2,len_AU,mass_sol,Mdot_MsMyr,nsmth=2,racc=4.5,overwrite=False)
        ax2b.set_yscale('log')
        #ax2b.set_yscale('linear')
        plt.legend(fontsize=10,labelspacing=0.03)
        #plt.axes().set_aspect('equal')
        plt.xlabel(r'$R$ (AU)',fontsize=20)
        plt.title('Time = %i kyr'%(time_list[iout]*1.e3), fontsize=24)
        plt.gcf().subplots_adjust(bottom=0.3)
        if contours: plt.savefig(path_out+'/Mdot_r_rhocon_'+str(ioutput).zfill(5)+'.pdf')
        plt.savefig(path_out+'/Mdot_r_'+str(ioutput).zfill(5)+'.pdf',bbox_inches='tight')
        f3 = plt.figure()
        if contours: plt.semilogy(range(scale),np.asarray(Mdot_tot)*Mdot_MsMyr)
        else: plt.plot(scale,np.asarray(Mdot_tot)*Mdot_MsMyr)
        plt.xlabel('z factor',fontsize=20); plt.ylabel(r'$\dot{M}_\mathrm{tot} (Ms/Myr)$',fontsize=20)
        plt.title('Time = %i kyr'%(time_list[iout]*1.e3), fontsize=24)
        if contours: plt.savefig(path_out+'/Mdot_zfactor_rhocon_'+str(ioutput).zfill(5)+'.pdf')
        plt.savefig(path_out+'/Mdot_zfactor_'+str(ioutput).zfill(5)+'.pdf')
        plt.show()
#    plt.show()
