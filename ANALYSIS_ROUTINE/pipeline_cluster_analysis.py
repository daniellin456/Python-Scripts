#! /usr/bin/env python


#from module_visualization import cal_rot_mat
import argparse
from pymses import RamsesOutput
from pymses.utils import constants as C
#from pymses.utils.regions import Sphere, Cylinder
#from pymses.filters import CellsToPoints, RegionFilter
import numpy as np
#from pymses.analysis.point_sampling import PointSamplingProcessor
#import multiprocessing as mp
import os
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import module_analysis as md_ana
import module_disc as md_disc
import module_statistics as md_stat
from module_ramses_io import write_pymsesrc
#from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis
from module_core import renormalize_units
#from module_analysis import cyl_samples, cylinder_flux
from module_visualization import save_fig,make_image_zoom

#############################################################################################
def density_PDF(ro, ioutput, writedir, overwrite=True,scl=12):
    filename = writedir+'/hist_rho_'+str(ioutput).zfill(5)+'.npz'
    #vmin2 = 1.e-6/mag_gauss; vmax2 = 1.e2/mag_gauss
    vmin = 1.e2; vmax=1.e12
    nbins = 100; nbins2 = 100
    if os.path.exists(filename) and not overwrite:
         fileload = np.load(filename)
         hist_rho = fileload["hist_rho"]; rho_edges = fileload["rho_edges"]#; B_edges = fileload["B_edges"]
    else:
        hist_rho, rho_edges = md_stat.pdf_var(ro, lmax=lmax, scale='log', var='density', weight='volume', bin_edges=None, vmin=vmin, vmax=vmax, nbins=nbins, ncpu=5, filt_frac=0.01, region=None)#, var2='B', vmin2=vmin2, vmax2=vmax2, nbins2=nbins2, scale2='log',bin_edges2=None)
#        hist_B_rho, rho_edges, B_edges = md_stat.pdf_var(ro, lmax=8, scale='log', var='density', weight='mass', bin_edges=None, vmin=None, vmax=None, nbins=10, ncpu=5, filt_frac=0.01, region=None, var2='B', vmin2=None, vmax2=None, nbins2=10, scale2='log',bin_edges2=None)
        np.savez(filename, hist_rho = hist_rho, rho_edges = rho_edges)
    #print hist_B_rho
    #hist_B_rho = np.ma.masked_where((hist_B_rho==0.), hist_B_rho)
    #im = plt.imshow(np.log10(hist_B_rho.T*mass_sol), origin='lower',extent=np.log10((rho_edges[0],rho_edges[-1],B_edges[0]*mag_gauss,B_edges[-1]*mag_gauss)),aspect=1,cmap='YlOrBr')
    #plt.ylim(-5,0)
    #plt.xlabel(r'$\mathrm{log}_{10}n (\mathrm{cm}^{-3})$',fontsize=18)
    #plt.ylabel(r'$\mathrm{log}_{10}B (\mathrm{G})$',fontsize=18)
    #plt.tick_params(axis='both', which='major', labelsize=14)
    #ax=plt.gca()
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
#    c = plt.colorbar(im,fraction=0.035,pad=0.02)
    #c = plt.colorbar(im,cax=cax)
    #c.set_clim(-7,1)
    #c.ax.get_yaxis().labelpad = 18
    #c.ax.tick_params(labelsize=14)
    #c.set_label(r'$\mathrm{log}_{10}M$',rotation=270,fontsize=18)#,labelpad=2.6)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #save_fig(writedir+'/hist2d_B_rho_'+str(ioutput).zfill(5),ps=True)
#    plt.show()

#############################################################################################
def plot_B_rho(ro, ioutput, writedir, overwrite=True,scl=12):
    filename = writedir+'/hist_B_rho_'+str(ioutput).zfill(5)+'.npz'
    vmin2 = 1.e-6/mag_gauss; vmax2 = 1.e2/mag_gauss
    vmin = 1.e4; vmax=1.e12
    nbins = 100; nbins2 = 100
    if os.path.exists(filename) and not overwrite:
         fileload = np.load(filename)
         hist_B_rho = fileload["hist_B_rho"]; rho_edges = fileload["rho_edges"]; B_edges = fileload["B_edges"]
    else:
        hist_B_rho, rho_edges, B_edges = md_stat.pdf_var(ro, lmax=lmax, scale='log', var='density', weight='mass', bin_edges=None, vmin=vmin, vmax=vmax, nbins=nbins2, ncpu=5, filt_frac=0.01, region=None, var2='B', vmin2=vmin2, vmax2=vmax2, nbins2=nbins2, scale2='log',bin_edges2=None)
#        hist_B_rho, rho_edges, B_edges = md_stat.pdf_var(ro, lmax=8, scale='log', var='density', weight='mass', bin_edges=None, vmin=None, vmax=None, nbins=10, ncpu=5, filt_frac=0.01, region=None, var2='B', vmin2=None, vmax2=None, nbins2=10, scale2='log',bin_edges2=None)
        np.savez(filename, hist_B_rho = hist_B_rho, rho_edges = rho_edges, B_edges = B_edges)
    print hist_B_rho
    hist_B_rho = np.ma.masked_where((hist_B_rho==0.), hist_B_rho) 
    im = plt.imshow(np.log10(hist_B_rho.T*mass_sol), origin='lower',extent=np.log10((rho_edges[0],rho_edges[-1],B_edges[0]*mag_gauss,B_edges[-1]*mag_gauss)),aspect=1,cmap='YlOrBr')
    plt.ylim(-5,0)
    plt.xlabel(r'$\mathrm{log}_{10}n (\mathrm{cm}^{-3})$',fontsize=18)
    plt.ylabel(r'$\mathrm{log}_{10}B (\mathrm{G})$',fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=14)
    ax=plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
#    c = plt.colorbar(im,fraction=0.035,pad=0.02)
    c = plt.colorbar(im,cax=cax)
    c.set_clim(-7,1)
    c.ax.get_yaxis().labelpad = 18
    c.ax.tick_params(labelsize=14)
    c.set_label(r'$\mathrm{log}_{10}M$',rotation=270,fontsize=18)#,labelpad=2.6)
    plt.gcf().subplots_adjust(bottom=0.15)
    save_fig(writedir+'/hist2d_B_rho_'+str(ioutput).zfill(5),ps=True)    
#    plt.show()
#############################################################################################
def cluster_visu(ro, ioutput, readdir, writedir, overwrite=True,scl=12):
    #radius = 100. # in AU
    #nbins = 30
    zoom_v = 0.5 / 4**np.arange(6)
#    zoom_v = np.append(zoom_v, zoom_v[1]*0.5)
#hydro_lowturb
#    c_center = np.array([ 0.49737453,  0.5016709,   0.4951839 ])#([0.5,0.5,0.5])
#    c_axis =  np.array([ 0.21040191,  0.85262068, -0.47829804])
#pint
#    c_center = np.array([ 0.5002806,   0.50004352,  0.49998823])
#    c_axis = np.array([-0.43124987, -0.24715896, -0.86771884])
    c_center = np.array([0.5,0.5,0.5])
    c_axis = None
    center_def = None #'baricenter'
    rotate = None #'AM'
    #rotate = 'AM'
    add_B = True
    #PSCAL = range(4)
    PSCAL=[]
    ps=False
    col_dens_only = True
    overwrite=True
    scl = 15
    #for iout in range(0,nout,20):#range(9,104,10):
   #    ioutput = outputs[iout]
    center_AU = make_image_zoom(readdir,ioutput,zoom_v,path_out=writedir,plot_sinks=True,overwrite=overwrite,center_img=c_center,make_phi=False,deli=None,col_dens_only=col_dens_only,tag='',gcomp=False,map_size=400,vel_size=20,xyz_sink=2,center_def=center_def,order='<',ps=ps,ind_sink=0,all_sink=False,rotate=rotate,make_slice=False,col_dens=True,c_axis=c_axis,add_B=add_B,PSCAL=PSCAL)
#############################################################################################
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Find ellipsoids for a series of radii")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
#    parser.add_argument("write_dir", type=str, help="write directory")
#    parser.add_argument("ncpu", type=int, help="number of cpu for pooling")
#    parser.add_argument("nsample", type=int, help="number of timestep samples")
#    parser.add_argument("scl", type=int, help="smallest cell level to read")
#    parser.add_argument("Rmin", type=float, help="minimum radius in code unit")
#    parser.add_argument("Rmax", type=float, help="maximum radius in code unit")
#    parser.add_argument("nstart",type=int,nargs='?',default=1,help="first output to calculate")
#    parser.add_argument("nend",type=int,nargs='?',default=-1,help="last output to calculate")
    args = parser.parse_args()

    #readdir = '/gpfs/data1/phennebe/'+args.read_dir
    #writedir = '/gpfs/data1/ylee/'+args.read_dir
    #readdir = '/gpfs/users/ynlee/'+args.read_dir
    #writedir = '/gpfs/users/ynlee/'+args.read_dir
    readdir = '/dsm/anais/storageA/ylee/'+args.read_dir
    writedir = '/dsm/anais/storageA/ylee/'+args.read_dir


    if not os.path.exists(writedir): os.makedirs(writedir)
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(readdir)
#    if args.ncpu<=1: pool = None
#    else: pool = mp.Pool(args.ncpu) #mp.cpu_count()
    #outputs = search_ro(readdir)
    outputs = np.loadtxt(readdir+'/sampled_outputs.txt')
    outputs = [int(out) for out in outputs]
    print 'unit_mag', mag_gauss, 'mass_sol', mass_sol
#    print 'outputs', outputs, args.nstart, args.nend
#    if args.nend==-1: args.nend=outputs[-1]
#    outputs = [output for output in outputs if output >= args.nstart]
#    outputs = [output for output in outputs if output <= args.nend]
#    print 'outputs', outputs
    nout = len(outputs)
    overwrite = False
#    overwrite = True
    scl = 20
    for iout in range(nout):
        plt.clf()
        ioutput = outputs[iout]
        ro = RamsesOutput(readdir, ioutput, order='<')
        #plot_B_rho(ro, ioutput, writedir,overwrite=overwrite,scl=scl)
        #density_PDF(ro, ioutput, writedir,overwrite=overwrite,scl=scl)
        cluster_visu(ro, ioutput, readdir, writedir, overwrite=overwrite,scl=scl)
