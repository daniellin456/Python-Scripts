#! /usr/bin/env python

## This module calculates statistical properties of a ramses output

#def var_op(cells, var='density'):
#def weight_op(cells, weight='volume'):
#def sum_op(cells, var='mass'):
#def do_pdf(cells, bins, var='density', weight='volume'):  # 1D pdf with pre-defined bins, weights could be multi-dimensional
#def pdf_var(ro, lmax=20, scale='log', var='density', weight='volume', vmin=None, vmax=None, nbins=40, ncpu=1, filt_frac=1, region=None): # cpu-wise pdf for large output
#def plot_hists(readdir, outputs, writedir=None, scl=20, scale='log',var='density',weight='volume', nbins=40, overwrite=False, ncpu=1,filt_frac=0.01,cumulated=False):
#def pdf_Ekin(ro, scl=20, scale='log', var='density', nbins=40, ncpu=1, filt_frac=1, region=None):
#def sum_var(ro, lmax=20, var='mass', ncpu=1, region=None): # cpu-wise sum for large output
#def calc_sigma_rho(readdir, outputs, writedir=None, scl=20, scale='log',var='density',weight='volume', nbins=40, overwrite=False, ncpu=1,filt_frac=0.01,cumulated=False):

##### pymses dependences ##################################################################
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
from pymses import RamsesOutput
from pymses.utils import constants as C
from pymses.utils.regions import Sphere, Cylinder
from pymses.filters import CellsToPoints, RegionFilter, PointFunctionFilter, PointRandomDecimatedFilter
#from mypymses_region import  RegionDifferentiator
#from mypymses import Ellipsoid
#from pymses.analysis.point_sampling import PointSamplingProcessor

##### python dependences ##################################################################
from functools import partial
import argparse
import os
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt

##### module dependences ##################################################################
from module_core import renormalize_units
from module_visualization import save_fig
#from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis
#from module_visualization import cal_rot_mat
#from module_core import normalisation

## module usage:
## boxlen, len_pc, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 

#############################################################################################
def var_op(cells, var='density',i=0,cen=np.array([[0.5,0.5,0.5]])):
    if var == 'density': return cells["rho"]
    elif var == 'temperature': return cells["T"]
    elif var == 'temperature2': return cells["P"]/cells["rho"]
    elif var == 'ekin': return cells["rho"]*cells["vel"]**2
    elif var == 'passive': return cells["rho"]*cells["passive"][:,i]
    elif var == 'radius': return np.linalg.norm(cells.points-cen,axis=1)
    elif var == 'B': return 0.5*np.linalg.norm(cells["Br"]+cells["Bl"],axis=1)
#############################################################################################
def weight_op(cells, weight='volume'):
    dx = cells.get_sizes()
    vol = dx**3
    if weight == "mass": return cells["rho"] * vol
    elif weight == "volume": return vol
    elif weight == "Ekin": return cells["rho"][:,np.newaxis] * vol[:,np.newaxis] * cells["vel"]**2
    elif weight == "Ekin0": return cells["rho"] * vol * cells["vel"][:,0]**2
    elif weight == "Ekin1": return cells["rho"] * vol * cells["vel"][:,1]**2
    elif weight == "Ekin2": return cells["rho"] * vol * cells["vel"][:,2]**2
    elif weight == "vel": return cells["rho"][:,np.newaxis] * vol[:,np.newaxis] * cells["vel"]
    else: return np.ones_like(vol)
#############################################################################################
def sum_op(cells, var='mass'):
    dx = cells.get_sizes()
    vol = dx**3.
    v2 = cells["vel"]**2
    if var == 'mass': return np.sum(cells["rho"] * vol)
    elif var == 'Ekin': return np.sum(cells["rho"] * vol * np.sum(v2,axis=1))
    elif var == 'Ekin3': return np.sum(cells["rho"][:,np.newaxis] * vol[:,np.newaxis] * v2, axis=0)
#############################################################################################
def make_bins(vmin,vmax,nbins,scale='log'):
    if scale=='log': return np.logspace(np.log10(vmin),np.log10(vmax),nbins+1)
    else: return np.linspace(vmin,vmax,nbins+1)
#############################################################################################
def do_pdf(cells, bins, var='density', weight='volume',cen=np.array([[0.5,0.5,0.5]])):  # 1D pdf with pre-defined bins, weights could be multi-dimensional
    data = var_op(cells, var=var,cen=cen)
    weights = weight_op(cells, weight=weight)
    if len(weights.shape)==1: pdf, bins = np.histogram(data,bins=bins,weights=weights)
    else:
        ndim = weights.shape[1] 
        pdf = []
        for idim in range(ndim):
            pdfi, bins = np.histogram(data,bins=bins,weights=weights[:,idim])     
            pdf.append(pdfi)
        pdf = np.asarray(pdf).T
    return pdf    
#############################################################################################
def do_pdf_2d(cells, bins, var='density', var2='radius', weight='volume',cen=np.array([[0.5,0.5,0.5]])):  # 2D pdf with pre-defined bins, weights could be multi-dimensional
    data = var_op(cells, var=var,cen=cen)
    data2 = var_op(cells, var=var2, cen=cen)
    weights = weight_op(cells, weight=weight)
    if len(weights.shape)==1: pdf, b1, b2 = np.histogram2d(data,data2,bins=bins,weights=weights)
    else:
        ndim = weights.shape[1]
        pdf = []
        for idim in range(ndim):
            pdfi, b1, b2 = np.histogram2d(data,data2,bins=bins,weights=weights[:,idim])
            pdf.append(pdfi)
        pdf = np.asarray(pdf)
        pdf = np.rollaxis(pdf,0,3)
    return pdf
#############################################################################################
def pdf_var(ro, lmax=20, scale='log', var='density', weight='volume', bin_edges=None, vmin=None, vmax=None, nbins=40, ncpu=1, filt_frac=1, region=None,cen=np.array([[0.5,0.5,0.5]]), var2=None, vmin2=None, vmax2=None, nbins2=5, scale2='log',bin_edges2=None): # cpu-wise pdf for large output
    # var: 'density', 'temperature', 'ekin'
    # weight: 'volume', 'mass', 'Ekin'
    svar = ["rho"] #,"passive"]
    if var in ['temperature2', 'pressure']: svar += ["P"]
    if var in ['temperature']: svar += ["T"]
    if var in ['ekin']: svar += ["vel"]
    if var in ['B']: svar += ["Bl","Br"]
    if var2 in ['temperature2', 'pressure']: svar += ["P"]
    if var2 in ['temperature']: svar += ["T"]
    if var2 in ['ekin']: svar += ["vel"]
    if var2 in ['B']: svar += ["Bl","Br"]
    if weight in ['Ekin', 'vel']: svar += ["vel"]
    if weight[:-1] in ['Ekin', 'vel']: svar += ["vel"]
    if weight in ['phi', 'g']: svar += ["g"]
    svar = list(set(svar))
    print svar
    amr = ro.amr_source(svar)
    if region is not None: amr = RegionFilter(region, amr) 
    cell_source = CellsToPoints(amr, smallest_cell_level=lmax)
    if bin_edges is None:
        if vmin is None or vmax is None: #undefined histogram range
            sub_cells = PointRandomDecimatedFilter(filt_frac, cell_source) 
            data = var_op(sub_cells.flatten(),var=var,cen=cen)
        #print 'maximum temperature: %.1f'%data.max()
        #data = var_op(sub_cells.flatten(),var='passive',i=0)
        #print 'maximum PSCAL: %.1f'%data.max()
        #return np.arange(5),np.arange(6)
            if vmin is None: vmin = data.min() / 1.01
            if vmax is None: vmax = data.max() * 1.01
        bin_edges = make_bins(vmin,vmax,nbins,scale=scale)
    else: nbins = len(bin_edges)-1
    if var2 is not None and bin_edges2 is None:
        if vmin2 is None or vmax2 is None: #undefined histogram range
            sub_cells = PointRandomDecimatedFilter(filt_frac, cell_source)
            data = var_op(sub_cells.flatten(),var=var2,cen=cen)
            if vmin2 is None: vmin = data.min() / 1.01
            if vmax2 is None: vmax = data.max() * 1.01
        bin_edges2 = make_bins(vmin2,vmax2,nbins2,scale=scale2)
    elif var2 is not None: nbins2 = len(bin_edges2)-1
    if ncpu<=1: 
        pool = None
        if var2 is None: func = lambda x: do_pdf(x, bin_edges, var=var, weight=weight)
        else: func = lambda x: do_pdf_2d(x, [bin_edges,bin_edges2], var=var, var2=var2, weight=weight)
        get_f = lambda x: x
    else: 
        pool = mp.Pool(ncpu)
        if var2 is None: func = lambda x: pool.apply_async(partial(do_pdf, bins=bin_edges, var=var, weight=weight),x)
        else: func = lambda x: pool.apply_async(partial(do_pdf_2d, bins=[bin_edges,bin_edges2], var=var, var2=var2, weight=weight),x)
        get_f = lambda x: x.get()
    if var2 is None: 
        pdf = np.zeros(nbins)
        if weight in ['Ekin','vel']: pdf = np.zeros((nbins,3))
    else: 
        pdf = np.zeros((nbins,nbins2))
        if weight in ['Ekin','vel']: pdf = np.zeros((nbins,nbins2,3))
    sub_pdf = map(func, cell_source.iter_dsets())
    for sub in sub_pdf:
        pdf += get_f(sub)
    if var2 is None: return pdf, bin_edges
    else: return pdf, bin_edges, bin_edges2
#############################################################################################
def plot_hists(readdir, outputs, writedir=None, scl=20, scale='log',var='density',weight='volume', nbins=40, overwrite=False, ncpu=1,filt_frac=0.01,cumulated=False,show=False):
    if writedir is None: writedir=readdir
    figname = writedir+'/pdf_'+weight+'_'+var
    if os.path.exists(figname+'.pdf') and not overwrite: return
    ind = readdir[::-1].find('/')
    if ind>0: titlepath = readdir[len(readdir)-ind:]
    else: titlepath = readdir
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
    unit_bin = {'density': 1, 'temperature': 1, 'temperature2': temperature_K, 'ekin': mass_kg * vel_ms**2 / (len_cm/1.e2)**3 * 0.5}
    unit_bin = unit_bin[var]
    unit_pdf = {'volume': (len_cm/1.e2)**3, 'mass': mass_kg, 'Ekin': mass_kg * vel_ms**2 * 0.5, 'vel': vel_ms, 'count': 1.}
    unit_pdf = unit_pdf[weight]
    label_x = {'density': 'cm$^{-3}$', 'temperature': 'temperature (K)', 'temperature2': 'temperature (K)', 'ekin': 'energy (J/m$^3$)'}
    label_y = {'volume': 'm$^3$', 'mass': 'mass (kg)', 'Ekin': 'kinetic energy (J)', 'vel': 'velocity (m/s)', 'count': 'N'}
    f_pdf = plt.figure()
    if cumulated: f_cdf = plt.figure()
    for ioutput in outputs:
        ro = RamsesOutput(readdir, ioutput, order='<')
        time = ro.info["time"]*time_Myr
        pdf, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var=var, weight=weight, nbins=nbins, ncpu=ncpu, filt_frac=filt_frac) 
        if scale=='log': bin_centers = np.sqrt(bin_edges[:-1]*bin_edges[1:])
        else: bin_centers = 0.5*(bin_edges[:-1]+bin_edges[1:])
        plt.figure(f_pdf.number)
        plt.loglog(bin_centers*unit_bin,pdf*unit_pdf,label='%i year'%(time*1.e6))
        if cumulated: 
            pdf = np.cumsum(pdf[::-1])[::-1]
            plt.figure(f_cdf.number)
            plt.semilogx(bin_centers*unit_bin,pdf*unit_pdf,label='%i year'%(time*1.e6))
    plt.figure(f_pdf.number)
    plt.legend(loc='best')
    plt.xlabel(label_x[var])
    plt.ylabel(label_y[weight])
    plt.title(titlepath+' '+weight+'-'+var)
    if not os.path.exists(figname+'.pdf') or overwrite: save_fig(figname, tight=False)
    if cumulated: 
        figname = figname+'_cum'
        plt.figure(f_cdf.number)
        plt.legend(loc='best')
        plt.title(titlepath+' '+weight+'-'+var)
        if not os.path.exists(figname+'.pdf') or overwrite: save_fig(figname, tight=False)
    if show: plt.show() 
#############################################################################################
def pdf_Ekin(ro, scl=20, scale='log', var='density', nbins=40, ncpu=1, filt_frac=1, region=None):
    pdf_Ekin, bin_edges = pdf_var(ro, lmax=scl, scale=scale, var=var, weight='Ekin', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac, region=region)
    pdf_Mv, bin_edges = pdf_var(ro, lmax=scl, scale=scale, var=var, weight='vel', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac,vmin=bin_edges[0],vmax=bin_edges[-1], region=region)
    pdf_M, bin_edges = pdf_var(ro, lmax=scl, scale=scale, var=var, weight='mass', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac,vmin=bin_edges[0],vmax=bin_edges[-1], region=region)
    return pdf_Ekin, pdf_Mv, pdf_M, bin_edges
#############################################################################################
def sum_var(ro, lmax=20, var='mass', ncpu=1, region=None): # cpu-wise sum for large output
    svar = ["rho"]
    if var in ['temperature', 'pressure']: svar += ["P"]
    if var in ['Ekin']: svar += ["vel"]
    if var in ['mag']: svar += ["Br", "Bl"]
    if var in ['phi', 'g']: svar += ["g"]
    svar = list(set(svar))
    amr = ro.amr_source(svar)
    if region is not None: amr = RegionFilter(region, amr)
    cell_source = CellsToPoints(amr, smallest_cell_level=lmax)
    if ncpu<=1:
        pool = None
        func = lambda x: sum_op(x, var=var)
        get_f = lambda x: x    
    else:
        pool = mp.Pool(ncpu)
        func = lambda x: pool.apply_async(partial(sum_op, var=var),x)
        get_f = lambda x: x.get()
    sum_var = 0.
    if var=='Ekin3': sum_var=np.zeros(3)
    sub_sum = map(func, cell_source.iter_dsets())
    for sub in sub_sum: sum_var += get_f(sub) 
    return sum_var
#############################################################################################
#############################################################################################
#############################################################################################
def calc_sigma_rho(readdir, outputs, writedir=None, scl=20, scale='log',var='density',weight='volume', nbins=40, overwrite=False, ncpu=1,filt_frac=0.01,cumulated=False, show=False):
    if writedir is None: writedir=readdir
    figname = writedir+'/sigma_rho_all'
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
#    unit_bin = 1.
    unit_bin = {'density': 1, 'temperature': temperature_K, 'ekin': mass_kg * vel_ms**2 / (len_cm/1.e2)**3 * 0.5}
    unit_bin = unit_bin[var]
    f_pdf = plt.figure()
#    if cumulated: f_cdf = plt.figure()
    for ioutput in outputs:
        #ro = RamsesOutput(readdir, ioutput, order='<')
        ro = RamsesOutput(readdir, ioutput)
        time = ro.info["time"]*time_Myr
        pdf_Ekin, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var='density', weight='Ekin', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac)
        pdf_Mv, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var='density', weight='vel', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac,vmin=bin_edges[0],vmax=bin_edges[-1])
       # print pdf_Ekin.shape
       # return
       # pdf_Ekin0, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var='density', weight='Ekin0', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac)
       # pdf_Ekin1, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var='density', weight='Ekin1', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac)
       # pdf_Ekin2, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var='density', weight='Ekin2', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac)
        pdf_M, bin_edges = pdf_var(ro, lmax=scl, scale=scale,var='density', weight='mass', nbins=nbins, ncpu=ncpu, filt_frac=filt_frac,vmin=bin_edges[0],vmax=bin_edges[-1])
        #sigma2_rho = pdf_Ekin/pdf_M[:,np.newaxis] - (pdf_Mv/pdf_M[:,np.newaxis])**2
        sigma2_rho = pdf_Ekin/pdf_M[:,np.newaxis] - (np.sum(pdf_Mv)/np.sum(pdf_M))**2
#        sigma2_rho0 = pdf_Ekin0/pdf_M
#        sigma2_rho1 = pdf_Ekin1/pdf_M
#        sigma2_rho2 = pdf_Ekin2/pdf_M
        bin_centers = np.sqrt(bin_edges[:-1]*bin_edges[1:])
        plt.figure(f_pdf.number)
        plt.loglog(bin_centers*unit_bin,np.sqrt(np.sum(sigma2_rho,axis=1))*vel_ms,label='%i year'%(time*1.e6))
#        plt.loglog(bin_centers*unit_bin,sigma2_rho0*vel_ms**2,label='%i year'%(time*1.e6))
#        plt.loglog(bin_centers*unit_bin,sigma2_rho1*vel_ms**2,label='%i year'%(time*1.e6))
#        plt.loglog(bin_centers*unit_bin,sigma2_rho2*vel_ms**2,label='%i year'%(time*1.e6))
#        if cumulated:
#            pdf = np.cumsum(pdf)
#            plt.figure(f_cdf.number)
#            plt.semilogx(bin_centers*unit_bin,pdf*unit_pdf,label='%i year'%(time*1.e6))
#    plt.figure(f_pdf.number)
    plt.legend(loc='best')
    plt.title(titlepath)
    if not os.path.exists(figname+'.pdf') or overwrite: save_fig(figname,tight=False)
#    if cumulated:
#        figname = figname+'_cum'
#        plt.figure(f_cdf.number)
#        plt.legend(loc='best')
#        plt.title(titlepath+' '+weight+'-'+var)
#        if not os.path.exists(figname+'.pdf') or overwrite: save_fig(figname)
    if show: plt.show()
#############################################################################################
def vel_dispersion(path,ioutput,path_out=None,gcomp=True):
    if path_out is None: path_out = path
    ro=pymses.RamsesOutput(path,ioutput)
    amr = ro.amr_source(["rho","phi","g"], grav_compat=gcomp)
#    list_max =[] #list of local maximum
#    nmax = 0
#    list_max_phi =[] #list of local maximum
#    nmax_phi = 0
#    mat_tidal=np.zeros((3,3))
    
    w_v=[]
    rho_v=[]
    x_v=[]
    y_v=[]
    z_v=[]
    dx_v=[]

    w2_v=[]
    rho2_v=[]

    ncell=0
    ncell_neg=0
    ncell_tid2_neg=0
    ncell_tid3_neg=0
    ncell2_neg=0
    ncell2_tid2_neg=0
    ncell2_tid3_neg=0

    for dset in amr.iter_dsets():
        ngrids = dset.amr_struct["ngrids"]
        ndim = dset.amr_header["ndim"]
        dset.amr_struct["neighbors"] = -np.ones((ngrids, 2*ndim), dtype='i')
        octree_compute_neighbors(dset)
        dset.get_grid_levels()
        dset.amr_struct["cell_centers"] = dset.get_cell_centers()
        dset.amr_struct["idomain_grid"] = dset.get_idomain_grid()
#        ndim = dset.amr_struct["ndim"]
        twotondim = 2**ndim
        ngridlevel = dset.amr_struct["ngridlevel"]
        ngridbound = dset.amr_struct["ngridbound"]
        ncpu = dset.amr_struct["ncpu"]
        nboundary = dset.amr_struct["nboundary"]
        son = dset.amr_struct["son_indices"]
        # Initialize array indexes
        iarr = 0

        # Loop over grids
        ilevel = 0
        igrid = 0
        ngridsremain = 0

        icpu = dset.icpu

        print 'ncell ',ncell , 'ncell_neg ',ncell_neg, 'ncell_tid2_neg ', ncell_tid2_neg , 'ncell_tid3_neg ', ncell_tid3_neg
        print 'ncell ',ncell , 'ncell2_neg ',ncell2_neg, 'ncell2_tid2_neg ', ncell2_tid2_neg , 'ncell2_tid3_neg ', ncell2_tid3_neg

        while(ilevel < dset.amr_struct["levelmax"]):
            print 'ilevel',ilevel
            dx = 0.5**ilevel

            # Skip all CPUs 0..icpu-1 for this level
            for kcpu in range(0, icpu): igrid += ngridlevel[kcpu, ilevel]

            # Loop over the CPU grids of this level
            ngridsremain = ngridlevel[icpu, ilevel]
            while(ngridsremain > 0):
                #select the oct which have only leaf cells
                find_leaf=True
                # Loop over cells to check whether all of them are leaf
                for ind in range(twotondim):
                    # Only stop at leaf cells
                    if (son[igrid, ind] >= 0): find_leaf=False

                #select the oct which have at least one leaf cell
                find_leaf=False
                # Loop over cells to check whether all of them are leaf
                for ind in range(twotondim):
                    # Only stop at leaf cells
                    if (son[igrid, ind] == -1): find_leaf=True

                if(find_leaf):
#                    print 'igrid',igrid
                    #position of the grid
                     xg = dset.amr_struct["grid_centers"][igrid,0]
                     yg = dset.amr_struct["grid_centers"][igrid,1]
                     zg = dset.amr_struct["grid_centers"][igrid,2]

                    #average over the oct and calculate the eigen values
                     g_xm = 0.25*(dset["g"][igrid,0]+dset["g"][igrid,2]+dset["g"][igrid,4]+dset["g"][igrid,6])
                     g_xp = 0.25*(dset["g"][igrid,1]+dset["g"][igrid,3]+dset["g"][igrid,5]+dset["g"][igrid,7])

                     g_ym = 0.25*(dset["g"][igrid,0]+dset["g"][igrid,1]+dset["g"][igrid,4]+dset["g"][igrid,5])
                     g_yp = 0.25*(dset["g"][igrid,2]+dset["g"][igrid,3]+dset["g"][igrid,6]+dset["g"][igrid,7])

                     g_zm = 0.25*(dset["g"][igrid,0]+dset["g"][igrid,1]+dset["g"][igrid,2]+dset["g"][igrid,3])
                     g_zp = 0.25*(dset["g"][igrid,4]+dset["g"][igrid,5]+dset["g"][igrid,6]+dset["g"][igrid,7])

                     mat_tidal[:,0] = (g_xp - g_xm) / dx
                     mat_tidal[:,1] = (g_yp - g_ym) / dx
                     mat_tidal[:,2] = (g_zp - g_zm) / dx

                     w,v=eig(mat_tidal)
                     w = np.sort(w)
                     w_v.append(w)

                     rho = np.sum(dset["rho"][igrid,:])/8.
                     rho_v.append(rho)
                     dx_v.append(dx)

                     x_v.append(xg)
                     y_v.append(yg)
                     z_v.append(zg)

                     ncell = ncell + 1.
                     if(np.sum(w) <= 0.):
                         ncell_neg = ncell_neg + 1
                     if(w[2]<=0.):
                         ncell_tid3_neg = ncell_tid3_neg + 1
                     if(w[1]<=0.):
                         ncell_tid2_neg = ncell_tid2_neg + 1

                    #average over the oct and calculate the eigen values
                     g_xm = dset["g"][igrid,0]
                     g_xp = dset["g"][igrid,1]
                     g_ym = dset["g"][igrid,0]
                     g_yp = dset["g"][igrid,2]
                     g_zm = dset["g"][igrid,0]
                     g_zp = dset["g"][igrid,4]

                     mat_tidal[:,0] = (g_xp - g_xm) / dx
                     mat_tidal[:,1] = (g_yp - g_ym) / dx
                     mat_tidal[:,2] = (g_zp - g_zm) / dx

                     w,v=eig(mat_tidal)
                     w = np.sort(w)
                     w2_v.append(w)
                     rho = dset["rho"][igrid,0]
                     rho2_v.append(rho)

                     if(np.sum(w) <= 0.):
                         ncell2_neg = ncell2_neg + 1
                     if(w[2]<=0.):
                         ncell2_tid3_neg = ncell2_tid3_neg + 1
                     if(w[1]<=0.):
                         ncell2_tid2_neg = ncell2_tid2_neg + 1
#                    for ind in range(twotondim):
#                        print dset["g"][igrid,ind]
#                mat_tidal[0,0] = dset["g"][igrid,ind][] - dset["g"][igrid,ind][]
#                if (dset["rho"][igrid,ind] > thres_dens and  dset.amr_struct["cell_levels"][igrid] >= thres_level):
#                    xg = dset.amr_struct["grid_centers"][igrid,0]
#                    yg = dset.amr_struct["grid_centers"][igrid,1]
#                    zg = dset.amr_struct["grid_centers"][igrid,2]
#                    print 'igrid: ',igrid,'ind: ',ind, dset["rho"][igrid,ind],dset["g"][igrid,ind]
#                     pos = dset.amr_struct["grid_centers"][igrid,:]
#                     pos_v.append(pos)
            # Go to next grid at this CPU and level
                ngridsremain -= 1
        # Skip all CPUs icpu+1, ... ncpu for this level
            for kcpu in range(icpu+1, ncpu):
                igrid += ngridlevel[kcpu, ilevel]
        # Skip all boundary grids for this level
            for kcpu in range(nboundary):
                igrid += ngridbound[kcpu, ilevel]
            ilevel += 1

    rho_v  = np.asarray(rho_v)
    w_v = np.asarray(w_v)
    rho2_v  = np.asarray(rho_v)
    w2_v = np.asarray(w_v)
    dx_v = np.asarray(dx_v)
    x_v = np.asarray(x_v)
    y_v = np.asarray(y_v)
    z_v = np.asarray(z_v)

#    print 'number of cells) ', len(wtot_v)
#    print 'number of cells with positive trace(tidal) ', np.sum(mask)

    name=path_out+'/v_dis_'+str(ioutput).zfill(5)+'.save'
    np.savez(name,rho_v=rho_v,w_v=w_v,x_v=x_v,y_v=y_v,z_v=z_v,w2_v=w2_v,rho2_v=rho2_v,dx_v=dx_v)
    
#############################################################################################
def display_vel_dispersion(path_in,num,path_out=None):


    if(path_out is None):
        directory_out=path_in
    else:
        directory_out=path_out

    name=path_in+'tidal_'+str(num)+'.save.npz'
    npzfile = np.load(name)

    rho_v = npzfile['rho_v']
    w_v   = npzfile['w_v']
    x_v = npzfile['x_v']
    y_v = npzfile['y_v']
    z_v = npzfile['z_v']

    dx_v = npzfile['dx_v']


    mask1 =  np.imag(w_v[:,1]) == 0.
    mask2 = w_v[:,0] < 0.
    mask3 = w_v[:,1] != 0.
    mask4 = w_v[:,2] != 0.

    mask=mask1*mask2*mask3*mask4

    rho_v = rho_v[mask]
    w1_v   = np.real(w_v[:,0][mask])
    w2_v   = np.real(w_v[:,1][mask])
    w3_v   = np.real(w_v[:,2][mask])

    dx_v = dx_v[mask]

    mass_v = ((dx_v)**3)*rho_v

    wtot_v = w1_v+w2_v+w3_v
    logwt_v = np.log10(abs(wtot_v))*np.sign(wtot_v)

    logw_w1_v = np.log10(abs(wtot_v/w1_v))*np.sign(wtot_v)

    logw1_v = np.log10(abs(w1_v))*np.sign(w1_v)
    logw2_v = np.log10(abs(w2_v))*np.sign(w2_v)
    logw3_v = np.log10(abs(w3_v))*np.sign(w3_v)


    logw13_v = np.log10(abs(w3_v/w1_v))

    logw12_v = np.log10(abs(w2_v/w1_v))


##1D histogram
    log_rho_min = 3.
    log_rho_max = 13.

    nbin=50.

    width = (log_rho_max - log_rho_min) / nbin

    P.clf()


    hist_mass , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v)

    P.bar(hist_edges[:-1],hist_mass,width)

    P.yscale('log', nonposy='clip')

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$\log(M) \, (code units)$')


    P.savefig(directory_out+'hist_mass'+str(num)+'.ps')
    P.savefig(directory_out+'hist_mass'+str(num)+'.jpeg')


    P.clf()

    hist_Tidal1 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*(w1_v))

    hist_Tidal1 = abs(hist_Tidal1) / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal1,width)
    P.yscale('log', nonposy='clip')

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$\log(\lambda1) \, (code units)$')


    P.savefig(directory_out+'hist_T1'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T1'+'_'+str(num)+'.jpeg')


    P.clf()

    hist_Tidal2 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*w2_v)

    hist_Tidal2 = abs(hist_Tidal2) / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal2,width)

    P.yscale('log', nonposy='clip')

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$\log(\lambda1) \, (code units)$')


    P.savefig(directory_out+'hist_T2'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T2'+'_'+str(num)+'.jpeg')


    P.clf()

    hist_Tidal3 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*w3_v)

    hist_Tidal3 = abs(hist_Tidal3) / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal3,width)

    P.yscale('log', nonposy='clip')

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$\log(\lambda1) \, (code units)$')


    P.savefig(directory_out+'hist_T3'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T3'+'_'+str(num)+'.jpeg')



    P.clf()
    P.bar(hist_edges[:-1],hist_Tidal2/hist_Tidal1,width)

#P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$\lambda2 / \lambda1) $')


    P.savefig(directory_out+'hist_T12'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T12'+'_'+str(num)+'.jpeg')




    P.clf()

    mask_pos = w2_v > 0.

    hist_Tidal2_pos , hist_edges = P.histogram(np.log10(rho_v[mask_pos]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_pos]*(w2_v[mask_pos]))


    mask_neg = w2_v < 0.
    hist_Tidal2_neg , hist_edges = P.histogram(np.log10(rho_v[mask_neg]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_neg]*abs(w2_v[mask_neg]))

    hist_Tidal2_rap = hist_Tidal2_pos / hist_Tidal2_neg

    P.bar(hist_edges[:-1],hist_Tidal2_rap,width)

    P.yscale('log', nonposy='clip')

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$r_2$')


    P.savefig(directory_out+'hist_T22'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T22'+'_'+str(num)+'.jpeg')




    P.clf()

    P.bar(hist_edges[:-1],hist_Tidal3/hist_Tidal1,width)

#P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$\lambda3 / \lambda1) $')


    P.savefig(directory_out+'hist_T13'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T13'+'_'+str(num)+'.jpeg')



    P.clf()

    w12_v = w1_v + w2_v
    hist_Tidal1p2 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*abs(w12_v))

    hist_Tidal1p2 = hist_Tidal1p2 / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal3/hist_Tidal1p2,width)

#P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$\lambda3 / (\lambda1+\lambda2)) $')


    P.savefig(directory_out+'hist_T1p23'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T1p23'+'_'+str(num)+'.jpeg')




    P.clf()

    mask_pos = w3_v > 0.

    hist_Tidal3_pos , hist_edges = P.histogram(np.log10(rho_v[mask_pos]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_pos]*(w3_v[mask_pos]))


    mask_neg = w3_v < 0.
    hist_Tidal3_neg , hist_edges = P.histogram(np.log10(rho_v[mask_neg]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_neg]*abs(w3_v[mask_neg]))

    hist_Tidal3_rap = hist_Tidal3_pos / hist_Tidal3_neg

    P.bar(hist_edges[:-1],hist_Tidal3_rap,width)

    P.yscale('log', nonposy='clip')

    P.xlabel(r'$\rho (cm^{-3})$')
    P.ylabel(r'$r_3$')

    P.savefig(directory_out+'hist_T33'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T33'+'_'+str(num)+'.jpeg')


#############################################################################################


#############################################################################################

