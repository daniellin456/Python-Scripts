#! /usr/bin/env python

##### pymses dependences ##################################################################
from pymses import RamsesOutput
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
#from pymses.utils import constants as C
#from pymses.utils.regions import Sphere, Cylinder
#from pymses.filters import CellsToPoints, RegionFilter
#from pymses.analysis.point_sampling import PointSamplingProcessor

##### python dependences ##################################################################
import argparse
import numpy as np
import os
#import matplotlib.pyplot as plt
#import multiprocessing as mp

##### module dependences ##################################################################
import module_analysis as md_ana
import module_disc as md_disc
import module_statistics as md_stat
import module_flux as md_flux
from module_local_dir import local_dir
import module_diffuse as md_df
#from module_visualization import cal_rot_mat
#from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis
#from module_core import normalisation
#from module_analysis import cyl_samples, cylinder_flux
#from module_visualization import save_fig


## module usage:
## boxlen, len_pc, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis, c_vel = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## side_xy, side_z, face_xy, face_r = cyl_samples(r_over_h,n_h,random=True)
## flux_r, flux_z_out, flux_z_in, r_shell, M_shell, face_surf, rad_surf = cylinder_flux(ro,c_center,c_axis,c_radius,c_height,n_shell=3,n_h=20,random_sample=True,scl=10)
#############################################################################################
def calculate_diffusion():
    md_df.evaluate_diffusion(readdir,pscal,overwrite=True,order='<',scl=20,plot_CDF=plot_CDF)
#############################################################################################


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Find ellipsoids for a series of radii")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()

    readdir, writedir = local_dir(args.read_dir)

    outputs = search_ro(readdir)
    nout = len(outputs)
    #outputs = outputs[1:38]
    overwrite = False
    scl = 20
    pscal = range(6)
    plot_CDF = True

    calculate_diffusion()    
#global_dynamics()
    #disc_flux()
    #write_grids()
#    detail_ana()

#    nsam = max(nout/args.nsample,1)
#    md_disc.time_series(readdir,path_out=writedir,overwrite=True,order='<',scl=scl,center_def='baricenter',rotate='AM')
#    md_disc.time_series(readdir,path_out=writedir,overwrite=True,order='<')
#    md_disc.disc_cylinder(readdir,path_out=writedir,overwrite=True,order='<')
    #md_disc.simple_disc_model(writedir)
    #md_disc.disc_characteristics(readdir,path_out=writedir,overwrite=True,order='<')
#    md_disc.disc_basic_parameters(readdir,path_out=writedir,overwrite=True,order='<',scl=20)
    #md_disc.calc_disc_properties(readdir,path_out=writedir,overwrite=True,order='<',scl=20)
    #md_disc.visulize_disc_flows(readdir,path_out=writedir,overwrite=True,order='<',scl=20)
#    md_disc.trace_radial_properties(readdir,path_out=writedir,overwrite=True,order='<',scl=20)
#    md_disc.trace_disc_movement(readdir,path_out=writedir,overwrite=True)
#    for iout in range(nout)[-100:]:
#        ioutput = outputs[iout]
#        print 'output', ioutput
#        # radius and height in AU
#        c_center, c_axis, c_vel = md_ana.plot_flux(readdir,ioutput,writedir=writedir,c_radius=60,c_hieghts=[30,20,10,5],scl=18,ns=24,nh=12,overwrite=overwrite)
#        if c_center is None: continue
#        md_disc.calc_rho_relations(readdir,ioutput,path_out=writedir,overwrite=False,order='<',scl=scl)#,center_img=c_center,c_axis=c_axis,c_vel=c_vel)
#        md_disc.mass_core_disc(readdir,ioutput,path_out=writedir,overwrite=True,order='<',scl=scl,center_img=c_center,c_axis=c_axis,c_vel=c_vel,disc_thres=1.e-13,core_thres=1.e-11)

    #outputs = outputs[-1::-nout/6]
    #md_stat.plot_hists(readdir, outputs, writedir=writedir, scl=scl, var='density', weight='volume', nbins=20,ncpu=10)
#    nbins = 24
#    ncpu = 20
#    filt_frac = 1
#    #var = 'density'
#    #weight = 'volume'
#    var = 'temperature'
#    weight = 'mass' 
#    cumulated = True
#    md_stat.plot_hists(readdir, outputs, writedir=writedir, scl=scl, var=var, weight=weight,nbins=nbins,ncpu=ncpu,filt_frac=filt_frac,cumulated=cumulated,overwrite=overwrite,show=True)
