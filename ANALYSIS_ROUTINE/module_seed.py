#! /usr/bin/env python

## The module flux contains analyses of disc flux properties

##### pymses dependences ##################################################################
#import pymses
#from pymses.filters import PointFunctionFilter
#from pymses import RamsesOutput
#from pymses.utils import constants as C
#from pymses.utils.regions import Sphere, Cylinder
#from pymses.filters import CellsToPoints, RegionFilter
#from pymses.analysis.point_sampling import PointSamplingProcessor
#from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
##### python dependences ##################################################################
#import argparse
import numpy as np
#import os
#import matplotlib.pyplot as plt
#import pickle as pck
#from scipy.interpolate import interp1d, interp2d
#from scipy.misc import derivative
#from scipy.optimize import curve_fit
#import scipy as sp
#import pdb
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import LineCollection
#from matplotlib.colors import ListedColormap, BoundaryNorm
#from scipy.signal import savgol_filter, argrelextrema
#import multiprocessing as mp
##### module dependences ##################################################################
from module_core import renormalize_units #, find_center_dmax, find_baricenter, find_axis, normalisation, read_sink_cvs
#from module_visualization import cal_rot_mat, save_fig
#from module_analysis import cylinder_flux, save_pickle, read_pickle
#from module_core import normalisation
## module usage:
## boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
Ms_sink = 1.9889e33
#############################################################################################

def write_points(readdir,nseeds, rs, rad,writedir=None,outs=[]):
    if writedir==None: writedir=readdir
    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(readdir)
    len_AU = len_pc * 2.e5
#    ind = readdir[::-1].find('/')
#    if ind>0: titlepath = path[len(path)-ind:]
#    else: titlepath = path
#np.savez(path_out+'/disc_basic_parameters', ioutput_list=ioutput_list, time_list=time_list, Msink_list=Msink_list, Mdisc_list=Mdisc_list, thres_rho_list=thres_rho_list, center_sink_list=center_sink_list, center_disc_list=center_disc_list, cvel_sink_list=cvel_sink_list, cvel_disc_list=cvel_disc_list, ax_csink_list=ax_csink_list, ax_ctot_list=ax_ctot_list, r_csink_list=r_csink_list, r_ctot_list=r_ctot_list)
    discload = np.load(writedir+'/disc_basic_parameters.npz')
    outputs = discload['ioutput_list']
    center_sink_list=np.asarray(discload['center_sink_list'])
    cvel_sink_list=np.asarray(discload['cvel_sink_list'])
    cvel_disc_list=np.asarray(discload['cvel_disc_list'])
    ax_csink = np.asarray(discload['ax_csink_list'])
    r_csink = np.asarray(discload['r_csink_list'])
    time_list = np.asarray(discload['time_list'])
    nout = len(outputs)
#    for i in range(nout):
#        print i, outputs[i], time_list[i]*1.e3
#    return

    seed_types = [1] * nseeds 
    if len(outs)>0:
        list_out = []
        for out in outs: list_out.append(np.argmin(abs(np.asarray(outputs)-out)))
    else: list_out = range(nout) 

    for iout in list_out:
        ax = ax_csink[iout,:]
        ax = ax / np.linalg.norm(ax)
        center = center_sink_list[iout,:]
        if center[0]==0. : continue
        rvec = cvel_sink_list[iout,:]-cvel_disc_list[iout,:]
        rvec = rvec - np.sum(ax*rvec) * ax 
        rvec = rvec / np.linalg.norm(rvec)
        rvec90 = np.cross(ax,rvec)
        print ax, rvec, rvec90
        rvec = rvec * np.cos(rad) + rvec90 * np.sin(rad)
        #rvec = np.array([np.cos(rad),np.sin(rad),0])
        #rvec = rvec - np.sum(ax*rvec) * ax 
        rvec = rvec / np.linalg.norm(rvec)

        print iout, outputs[iout], time_list[iout]*1.e3
        print center, ax, rvec
        f = open(writedir+"/seeds_passive_%i_%i.dat"%(outputs[iout],rad/np.pi*180),"w")
        f.write(str(nseeds)+"\n")
        f.write("  ".join(map(str, seed_types))+"\n")
        for i in range(nseeds):
            seed = center + rvec * rs[i]/len_AU
            f.write("  ".join(map(str, seed))+"\n")
        f.close()

         

#lmax = 15
#ax = [0,0,1]
#rvec = [1,0,0]
#center = [0.5,0.5,0.5]
#dx = 1./2.**lmax
#len_AU = 15312.4238652
#
#nseeds = 4
#seed_types = [1, 1, 1, 1]
#dis = [4, 6, 8, 10] # in AU
#
#ax = np.asarray(ax)
#rvec = np.asarray(rvec)
#center = np.asarray(center)
#rvec = rvec - np.sum(ax*rvec) * ax / np.sum(ax**2)
#rvec = rvec / np.linalg.norm(rvec)
#
#f = open("seed_passive.dat","w")
#f.write(str(nseeds)+"\n")
#f.write("  ".join(map(str, seed_types))+"\n")
#for i in range(nseeds):
#    seed = center + rvec * dis[i]/len_AU
#    f.write("  ".join(map(str, seed))+"\n")
#f.close()
