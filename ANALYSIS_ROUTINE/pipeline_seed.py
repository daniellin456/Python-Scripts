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
#import module_analysis as md_ana
#import module_disc as md_disc
#import module_statistics as md_stat
#import module_flux as md_flux
from module_local_dir import local_dir
import module_seed as md_sd

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
def write_point_seeds(nseeds, rs, rad,outs=[]):
    md_sd.write_points(readdir,nseeds, rs, rad,writedir=writedir,outs=outs)

#############################################################################################

#############################################################################################

#############################################################################################

#############################################################################################


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Find ellipsoids for a series of radii")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()

    readdir, writedir = local_dir(args.read_dir)

    outputs = search_ro(readdir)
    nout = len(outputs)
    overwrite = True
    scl = 20

    nseeds = 6
    rs = [4,6,8,10,14,18] # in AU
    rad = 270 # angle
    outs = [600,800,1000]
    outs = [500,700]
    for rad in [0,90,180,270]:
        write_point_seeds(nseeds, rs, rad/180.*np.pi,outs=outs)
