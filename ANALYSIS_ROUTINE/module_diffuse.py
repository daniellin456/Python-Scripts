#! /usr/bin/env python

## The module flux contains analyses of disc flux properties

##### pymses dependences ##################################################################
import pymses
from pymses.filters import PointFunctionFilter
from pymses import RamsesOutput
from pymses.utils import constants as C
#from pymses.utils.regions import Sphere, Cylinder
from pymses.filters import CellsToPoints #, RegionFilter
#from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
##### python dependences ##################################################################
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle as pck
from scipy.interpolate import interp1d, interp2d
#from scipy.misc import derivative
from scipy.optimize import curve_fit
from functools import partial
import scipy as sp
#import pdb
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import LineCollection
#from matplotlib.colors import ListedColormap, BoundaryNorm
#from scipy.signal import savgol_filter, argrelextrema
#import multiprocessing as mp
##### module dependences ##################################################################
from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis, normalisation, read_sink_cvs
from module_visualization import cal_rot_mat, save_fig
from module_analysis import cylinder_flux, save_pickle, read_pickle
#from module_core import normalisation
## module usage:
## boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
Ms_sink = 1.9889e33
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
def gaussian(n0, x,c, sig):
    val = n0 * np.exp(-(x-c)**2 / sig**2/ 2.) / np.sqrt(2*np.pi) / sig
    return val
def inte_gau(n0,x,c,sig):
    val = n0 / 2. * sp.special.erfc((c-x)/sig/np.sqrt(2))
    return val
#############################################################################################
def gaussian_fit(mass,posi,center,ax,m0,dx,len_AU,writedir,tag,plot=True,plot_show=False):
    ax = ax / np.linalg.norm(ax)
    print 'ax_sink', ax
    posi = posi-center[np.newaxis,:]
    z = np.sum(posi * ax[np.newaxis,:], axis=1)
    r = np.sqrt( np.sum(posi**2,axis=1) - z**2 )
    fit = []
#    posi_mean = np.sum(mass[:,np.newaxis] * posi,axis=0) / np.sum(mass)
#    posi_var = np.sum(mass[:,np.newaxis,np.newaxis] * posi[:,np.newaxis,:] * posi[:,:,np.newaxis], axis=0) / np.sum(mass)
#    w, v = np.linalg.eigh(posi_var)
#    ind_eig = np.argmax(w)
#    print 'eigenvalues', w
#    ax_d = np.cross(posi_mean, v[:,ind_eig])
#    ax_d = ax_d/np.linalg.norm(ax_d)
#    print 'ax_pas', ax_d, np.sum(ax*ax_d)
#    z = np.sum(posi * ax_d[np.newaxis,:], axis=1)
#    r = np.sqrt( np.sum(posi**2,axis=1) - z**2 )


    for dist, dir in zip([r,z],['r_','z_']): 
        order = np.argsort(dist)
        mass_o = mass[order]/m0; dist = dist[order]
        mass_cum = np.cumsum(mass_o)
        mass_acc = 1. - mass_cum[-1]
        mass_cum += mass_acc
        mass_o = np.append(mass_acc, mass_o)
        mass_cum = np.append(mass_acc, mass_cum); dist = np.append(dist[0]*2-dist[1],dist)
        dist_mass = interp1d(mass_cum,dist)
    #    print mass_cum, m0/2.
        if m0/2.> mass_cum[0]:
            dist_mass = interp1d(mass_cum,dist)
            d0 = dist_mass(1./2.)
        else: d0 = mass_cum[0]
        sig = np.sqrt( np.sum(mass_o*(dist-d0)**2)) 
        sig_low = 0.
        if sig/dx > sig_low:  popt, pcov = curve_fit(partial(inte_gau,1.),dist,mass_cum,p0=[d0,sig])
        else:
            zmax = np.ceil((np.amax(dist)/dx)/0.5)
            zmin = np.floor((np.amin(dist)/dx)/0.5)    
            dist_bins = np.arange(zmin,zmax+1,1) * dx
            #dist_bins[0] = np.amin(dist) - dx; dist_bins[-1] = np.amax(dist) + dx
            hist_mass = np.histogram(dist, bins = dist_bins, weights=mass_o)[0]
            mass_cum_bin = np.cumsum(hist_mass)
            popt2, pcov = curve_fit(partial(inte_gau,1.),dist_bins[1:],mass_cum_bin,p0=[d0,sig])
            print popt2[0]*len_AU, d0*len_AU
            print popt2[1]*len_AU, sig*len_AU
            popt, pcov = curve_fit(partial(inte_gau,1.),dist,mass_cum,p0=[d0,sig])
        print popt[0]*len_AU, d0*len_AU
        print popt[1]*len_AU, sig*len_AU
        if plot:
            plt.clf()
            plt.plot(dist/dx, mass_cum)
            plt.plot(dist/dx, inte_gau(1.,dist,popt[0],popt[1]))
            if sig/dx <= sig_low: plt.plot(dist_bins[1:]/dx, mass_cum_bin)
            if sig/dx <= sig_low: plt.plot(dist/dx, inte_gau(1.,dist,popt2[0],popt2[1]))
            plt.xlabel('dx')
            plt.ylabel('CDF(m)')
            plt.savefig(writedir+'/CDF_diffusion_'+dir+tag+'.pdf')
        #plt.plot(dist/dx, mass_o/np.append(dist[1]-dist[0],np.diff(dist)))
        #plt.plot(dist/dx, gaussian(1.,dist,popt[0],popt[1]))
#        if popt[1] < dx / 2.:
#            dist = np.hstack((dist,dist[1:]-0.25*dx,dist[1:]+0.25*dx))
#            mass_o = np.hstack(( mass_o[0], mass_o[1:]*0.5, mass_o[1:]*0.25, mass_o[1:]*0.25))
#            order = np.argsort(dist)
#            mass_o = mass_o[order]; dist = dist[order]
#            mass_cum = np.cumsum(mass_o)
#            popt, pcov = curve_fit(partial(inte_gau,1.),dist,mass_cum,p0=[d0,sig])
#            print popt[0]*len_AU, d0*len_AU
#            print popt[1]*len_AU, sig*len_AU
#            plt.plot(dist/dx, mass_cum)
#            plt.plot(dist/dx, inte_gau(1.,dist,popt[0],popt[1]))
#            #plt.plot(dist/dx, mass_o/np.append(dist[1]-dist[0],np.diff(dist)))
#            #plt.plot(dist/dx, gaussian(1.,dist,popt[0],popt[1]))
        if plot_show: plt.show()
        fit.append(popt)
    return fit

#############################################################################################
def calc_alpha(amr,mass_sink,pos_sink,vel_sink,thres_rho,dx,len_AU,boxlen):
    rhofunc = lambda dset: dset["rho"] > thres_rho
    cells = CellsToPoints(amr)
    cells = PointFunctionFilter(rhofunc,cells)
    cells = cells.flatten()
    posi = cells.points - pos_sink[np.newaxis,:]
    vel = cells["vel"] - vel_sink[np.newaxis,:]
    dx_cell = cells.get_sizes()
    dV = dx_cell**3
    dm = cells["rho"] * dV
    ax_d = np.sum(dm[:,np.newaxis] * np.cross(posi,vel),axis=0)
    ax_d = ax_d / np.linalg.norm(ax_d)
    zd = np.sum(posi * ax_d[np.newaxis,:],axis=1)
    rd = np.sqrt(np.sum(posi**2,axis=1)-zd**2)
    rmax = np.amax(rd)
    print 'rmax', rmax, dx
    rr = np.arange(0,rmax+dx,dx)
    rc = 0.5 * ( rr[:-1]+rr[1:] )
    hist_x = lambda x: np.histogram(rd, bins=rr, weights=x)[0] 
    mzz = hist_x(dm*zd**2)
    mass = hist_x(dm)
    hd = np.sqrt(mzz/mass)  ## disc scale height
    mcc = hist_x(dV * cells["P"])
    cs = np.sqrt(mcc/mass)
    print 'boxlen = ', boxlen
    hp = rc**1.5 * cs / boxlen / np.sqrt(mass_sink)
    print hp / hd
    hmask = hp.copy() * 2
    for ir in range(len(rr)):
        if hp[ir]<dx: hmask[ir] =hd[ir] * 2
        else: break
    print hp * 2 * len_AU
    print hmask * len_AU
    func_h = interp1d(rc,hmask,fill_value='extrapolate')
    #plt.plot(hp*len_AU,hd*len_AU)
    #plt.show()
    vecr = posi / np.linalg.norm(posi,axis=1)[:,np.newaxis]
    vecphi = np.cross(ax_d, vecr)
    vr = np.sum( vel * vecr, axis=1)
    vphi = np.sum( vel * vecphi, axis=1)
    #mask = np.where(abs(zd)<=func_h(rd))
    #print mask
    mask = np.where(abs(zd)>func_h(rd))
    dm2 = dm.copy(); dm2[mask]=0.
    dV2 = dV.copy(); dV2[mask]=0.
    mvrvphi = hist_x( (dm2 * vr * vphi) )
    mvr = hist_x( (dm2 * vr) ) 
    mvphi = hist_x( (dm2 * vphi) )
    mass2 = hist_x(dm2)
    mdvrdvphi = mvrvphi - mvr * mvphi / mass2
    mcc2 = hist_x(dV2 * cells["P"])
    B_disc = (cells["Br"]+cells["Bl"])*0.5
    B_r = np.sum(B_disc * vecr, axis=1)
    B_phi = np.sum(B_disc * vecphi, axis=1)
    B_z = np.sum(B_disc * ax_d, axis=1)
    g_r = np.sum(cells["g"] * vecr, axis=1)
    g_phi = np.sum(cells["g"] * vecphi, axis=1)
    VBrBphi = hist_x(dV2 * B_r * B_phi)
    VBzBphi = hist_x(dV2 * B_z * B_phi)
    Vgrgphi = hist_x(dV2 * g_r * g_phi)

#    mvr2 = hist_x( dm * vr**2)
#    sig2r = (mvr2 - mvr**2/mass) / mcc
#    mvphi2 = hist_x( dm * vphi**2)
#    sig2phi = (mvphi2 - mvphi**2/mass) / mcc
    alpha = mdvrdvphi/mcc2
    alpha_Max = -VBrBphi / mcc2 #/ (4.*np.pi)
    alpha_Maxz = -VBzBphi / mcc2 #/ (4.*np.pi)
    alpha_Poi = Vgrgphi / mcc2 / (4.*np.pi)
#    pp = np.linspace(-np.pi,np.pi,16)
#            yp = np.cross(vec_z, (center_disc_list[iout]-center))
#        yp = yp/np.linalg.norm(yp)
#        xp = np.cross(yp,vec_z)
#        xp = xp/np.linalg.norm(xp)
#        x_disc = np.sum(post_disc*xp[np.newaxis,:],axis=1)
#        y_disc = np.sum(post_disc*yp[np.newaxis,:],axis=1)
#        phi = np.arctan2(y_disc,x_disc)
    plt.figure()
    plt.plot(rr[:-1]*len_AU,alpha,label=r'$\alpha_{Rey}$')
    plt.plot(rr[:-1]*len_AU,alpha_Max,label=r'$\alpha_{Max}$')
#    plt.plot(rr[:-1]*len_AU,alpha_Maxz,label=r'$\alpha_{Max,z}$')
#    plt.plot(rr[:-1]*len_AU,alpha_Poi,label=r'$\alpha_{Grav}$')
#    plt.plot(rr[:-1]*len_AU,alpha+alpha_Max+alpha_Poi)
#    plt.plot(rr[:-1]*len_AU,sig2r,label=r'$\sigma^2_r$')
#    plt.semilogy(rr[:-1]*len_AU,sig2phi,label=r'$\sigma^2_\phi$')
    plt.yscale('symlog',linthreshy=1.e-3)
    plt.legend()
    plt.show()
    return  ax_d, rr, alpha, cs, hd
#############################################################################################
def evaluate_diffusion(path,pscal,path_out=None,overwrite=True,order='<',scl=20,plot_CDF=False):
    if path_out == None: path_out = path
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    len_AU = len_pc * 2.e5
    Ms_sink = 1.9889e33
    ind = path[::-1].find('/')
    if ind>0: titlepath = path[len(path)-ind:]
    else: titlepath = path
    npscal = len(pscal)
    m0 = np.zeros(npscal)

    outputs = search_ro(path)
    nout = len(outputs)
    ro = pymses.RamsesOutput(path,outputs[0],order=order)
    time0 = ro.info["time"]*time_Myr
    dx = 1./2**lmax
    
    ind_seed = path_out.find('_seed')
    ind_ori = path_out[ind_seed-1::-1].find('_')
    path_origin = path_out[:ind_seed-ind_ori-1]
    print path_origin
    discload = np.load(path_origin+'/disc_basic_parameters.npz')
    outputs_ori = discload['ioutput_list']
    thres_rho_list = discload['thres_rho_list']
    thres_rho = thres_rho_list[np.where(outputs_ori==outputs[0])]
    nend = min(nout,10)
    time = np.zeros(nend)
    time[0] = time0
    dist = np.zeros((nend,npscal,2))
    dist[0,:,0] = np.array([4,6,8,10,14,18]) / len_AU
    sig = np.zeros((nend,npscal,2))
    sig[0,:,:] = dx
    ax_d = np.zeros((nend,3))
    rr = [[]]*nend
    alpha = [[]]*nend
    cs = [[]]*nend
    hd = [[]]*nend
    istart = 0
    nend0a = 0
    for rad in ['0','90','180','270']:
        file_alpha = path_out[:ind_seed+6] + rad + '/new_alpha_r_t.npz' 
#        file_alpha = path_out[:ind_seed+6] + rad + '/alpha_r_t.npz'
        print file_alpha
        if os.path.exists(file_alpha):
            alpha_record = np.load(file_alpha)
            nend0a = len(alpha_record["rr"])
            nend_cp = min(nend, nend0a)
            rr[:nend_cp] = alpha_record["rr"][:nend_cp]
            alpha[:nend_cp] = alpha_record["alpha"][:nend_cp]
            cs[:nend_cp] = alpha_record["cs"][:nend_cp]
            hd[:nend_cp] = alpha_record["hd"][:nend_cp]
            ax_d[:nend_cp,:] = alpha_record["ax_d"][:nend_cp,:]
            istart = nend_cp
            break
    for iout in range(istart,nend):
        ioutput = outputs[iout]
        sink = read_sink_cvs(ioutput,path,deli=',')
        mass_sink = sink[0,1] * Ms_sink /  mass_kg * 1.e-3
        print 'Msink =', mass_sink, sink[0,1]
        print (len_cm*1.e-2)**0.5 * vel_ms / np.sqrt(Ms_sink * 1.e-3*6.67e-11) 
        print  (len_cm*1.e-2)**0.5 * vel_ms*boxlen / np.sqrt(mass_kg * 6.67e-11)
        print boxlen
        print Ms_sink, mass_kg/mass_sol
        pos_sink = sink[0,3:6]/boxlen
        vel_sink = sink[0,6:9]
        ax_sink = sink[0,10:13]
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        time[iout] = ro.info["time"] * time_Myr
        amr = ro.amr_source(["rho","vel","P","Bl",'Br',"g"])
        ax_d[iout,:],rr[iout], alpha[iout], cs[iout], hd[iout] = calc_alpha(amr,mass_sink,pos_sink,vel_sink,thres_rho,dx,len_AU,boxlen)
    if nend > nend0a: np.savez(file_alpha,time=time,ax_d=ax_d,rr=rr,alpha=alpha,cs=cs,hd=hd)


    file_diff = path_out+'/diffusion_t.npz'
    istart = 1
    nend0 = 0
    if os.path.exists(file_diff): 
        diff_record = np.load(file_diff)
        nend0 = len(diff_record["time"])
        nend_cp = min(nend,nend0)
        time[1:nend_cp] = diff_record["time"][1:nend_cp]
        dist[1:nend_cp,:,:] = diff_record["dist"][1:nend_cp,:,:]
        sig[1:nend_cp,:,:] = diff_record["sig"][1:nend_cp,:,:]
        istart = nend_cp
       
       
    for iout in range(istart,nend):
        ioutput = outputs[iout]
        sink = read_sink_cvs(ioutput,path,deli=',')
        mass_sink = sink[0,1] 
        pos_sink = sink[0,3:6]/boxlen
        vel_sink = sink[0,6:9]
        ax_sink = sink[0,10:13]
        ro = pymses.RamsesOutput(path,ioutput,order=order)
        time[iout] = ro.info["time"] * time_Myr
        amr = ro.amr_source(["rho","passive"])
        psfunc = lambda dset: np.sum(dset["passive"],axis=1) >1.e-8
        #amr = PointFunctionFilter(psfunc,amr)
        cells = CellsToPoints(amr)#.flatten()
        cells = PointFunctionFilter(psfunc,cells)
        cells = cells.flatten()
#        posi = cells.points - pos_sink[np.newaxis,:]
#        dx_cell = cells.get_sizes()
#        dm = cells["rho"] * dx_cell**3
#        mass = dm * np.sum(cells["passive"],axis=1)
#        posi_var = np.sum(mass[:,np.newaxis,np.newaxis] * posi[:,np.newaxis,:] * posi[:,:,np.newaxis], axis=0) / np.sum(mass)
#        w, v = np.linalg.eigh(posi_var)
#        ind_eig = np.argmin(w)
#        print 'eigenvalues', w
#        ax_d = v[:,ind_eig]
#        ax_d = ax_d/np.linalg.norm(ax_d)
#        print 'ax_pas', ax_d, ax_sink, np.sum(ax_d*ax_sink)
#        ax_d = ax_d * np.sign(np.sum(ax_d*ax_sink))
    #    print cells["passive"]
    #    print cells["passive"].shape
        for ip in pscal:
            mask = (cells["passive"][:,ip]>1.e-8)
        #    print mask
            psfunc = lambda dset: dset["passive"][:,ip] >1.e-8
            cells_i = PointFunctionFilter(psfunc,cells).flatten()
            #cells_i = cells*mask
            #print cells_i, cells_i.shape
    #        print cells_i["passive"]
    #        print cells_i["passive"].shape
            dx_cell = cells_i.get_sizes()
            dm = cells_i["rho"] * dx_cell**3
            mass = dm * cells_i["passive"][:,ip] 
            if m0[ip]==0: m0[ip] = np.sum(mass)
            posi = cells_i.points
            tag = str(ioutput).zfill(5)+'_'+str(ip)
            fit =  gaussian_fit(mass, posi, pos_sink, ax_d[iout,:], m0[ip],dx,len_AU,path_out,tag,plot=plot_CDF,plot_show=False)
            dist[iout,ip,0] = fit[0][0]; dist[iout,ip,1] = fit[1][0]
            sig[iout,ip,0] = fit[0][1]; sig[iout,ip,1] = fit[1][1]
    if nend > nend0: np.savez(file_diff,time=time,dist=dist,sig=sig)
#    return
   
    Ddiff = np.zeros((nend,npscal,2)) 
    istart = 1
    for i in range(istart,nend):
        for ip in range(npscal):
            fitD = np.polyfit(time[istart-1:i+1],sig[istart-1:i+1,ip,0]**2,1)
        #    print -fitD[1]/fitD[0]-time0, time0
            Ddiff[i, ip,0] = fitD[0]
            fitD = np.polyfit(time[istart-1:i+1],sig[istart-1:i+1,ip,1]**2,1)
        #    print -fitD[1]/fitD[0]-time0, time0
            Ddiff[i, ip,1] = fitD[0]
    #return
    linear = False
    plt.figure()
    plt.loglog((time-time0)*1.e3,sig[:,:,0]*len_AU,label='radial')
    plt.loglog((time-time0)*1.e3,np.sqrt(time-time0)*1.e2,'k')
    if linear: plt.xscale('linear');plt.yscale('linear')
    plt.xlabel('time(kyr)')
    plt.ylabel(r'$\sigma_r$(AU)')
    plt.figure()
    plt.semilogx((time-time0)*1.e3,dist[:,:,0]*len_AU)
    if linear: plt.xscale('linear');plt.yscale('linear')
    plt.xlabel('time(kyr)')
    plt.ylabel('r(AU)')
    plt.figure()
    plt.loglog((time-time0)*1.e3,sig[:,:,1]*len_AU,label='vertical')
    plt.loglog((time-time0)*1.e3,np.sqrt(time-time0)*1.e2,'k')
    if linear: plt.xscale('linear');plt.yscale('linear')
    plt.xlabel('time(kyr)')
    plt.ylabel(r'$\sigma_z$(AU)')
    plt.figure()
    plt.semilogx((time-time0)*1.e3,dist[:,:,1]*len_AU)
    if linear: plt.xscale('linear');plt.yscale('linear')
    plt.xlabel('time(kyr)')
    plt.ylabel('z(AU)')
    rr_c = []
    for r in rr:
        rr_c.append(0.5 * (r[:-1]+r[1:]) * len_AU)
    plt.figure()
    for i in range(nend):
        plt.plot(rr_c[i],alpha[i])
    plt.yscale('symlog',linthreshy=1.e-2)
    plt.xlabel('r(AU)'); plt.ylabel(r'$\alpha$')
    plt.figure()
    for i in range(nend):
        plt.plot(rr_c[i],cs[i]*vel_ms)
    plt.xlabel('r(AU)'); plt.ylabel(r'$c_s$')
    plt.figure()
    for i in range(nend):
        plt.plot(rr_c[i],hd[i]*len_AU)
    plt.xlabel('r(AU)'); plt.ylabel(r'$H_d$')

    plt.figure()
    for i in range(nend):
        plt.plot(rr_c[i],-hd[i]*cs[i]*alpha[i]*vel_ms*len_cm*1.e-2)
    plt.gca().set_color_cycle(None)
    for i in range(1,nend):
#        plt.plot(dist[i,:,0]*len_AU,sig[i,:,0]**2/(time[i,np.newaxis]-time0)*time_Myr*vel_ms*len_cm*1.e-2,'--')
        plt.plot(dist[i,:,0]*len_AU,Ddiff[i,:,0]*time_Myr*vel_ms*len_cm*1.e-2*boxlen,'--')
    plt.gca().set_color_cycle(None)
    for i in range(1,nend):
#        plt.plot(dist[i,:,0]*len_AU,sig[i,:,1]**2/(time[i,np.newaxis]-time0)*time_Myr*vel_ms*len_cm*1.e-2,':')
        plt.plot(dist[i,:,0]*len_AU,Ddiff[i,:,1]*time_Myr*vel_ms*len_cm*1.e-2*boxlen,':')
    plt.xlabel('r(AU)'); plt.ylabel(r'$\alpha c_s h_d$')
    plt.yscale('log')
#    print Ddiff*time_Myr*vel_ms*len_cm*1.e-2*boxlen
#    print Ddiff*time_Myr*vel_ms*len_cm*1.e-2/boxlen
#    print -hd[i]*cs[i]*alpha[i]*vel_ms*len_cm*1.e-2
    plt.show()

