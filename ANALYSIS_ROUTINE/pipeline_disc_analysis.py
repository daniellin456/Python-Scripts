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
import matplotlib.pyplot as plt
#import multiprocessing as mp

##### module dependences ##################################################################
import module_analysis as md_ana
import module_disc as md_disc
import module_statistics as md_stat
import module_flux as md_flux
from module_local_dir import local_dir
#from module_visualization import cal_rot_mat
#from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis
#from module_core import normalisation
#from module_analysis import cyl_samples, cylinder_flux
#from module_visualization import save_fig
from module_analysis import save_pickle, read_pickle

## module usage:
## boxlen, len_pc, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc, unit_J, time_Myr, unit_Ek, unit_Eg0, unit_Eg1, unit_Ep, unit_Eb = renormalize_units(readdir)
## c_center = find_center(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## c_axis, c_vel = find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20) 
## side_xy, side_z, face_xy, face_r = cyl_samples(r_over_h,n_h,random=True)
## flux_r, flux_z_out, flux_z_in, r_shell, M_shell, face_surf, rad_surf = cylinder_flux(ro,c_center,c_axis,c_radius,c_height,n_shell=3,n_h=20,random_sample=True,scl=10)
#############################################################################################
def threshold_disc(outputs=False):
    scl = 18
    overwrite=False
#    for iout in range(nout)[300::20]:
#        ioutput = outputs[iout]
#        c_center, c_axis, c_vel = md_ana.plot_flux(readdir,ioutput,writedir=writedir,c_radius=60,c_hieghts=[30,20,10,5],scl=scl,ns=24,nh=12,overwrite=overwrite)
#        #if c_center is None: continue
#        md_disc.calc_rho_relations(readdir,ioutput,path_out=writedir,overwrite=False,order='<',scl=scl)#,center_img=c_center,c_axis=c_axis,c_vel=c_vel)
#    outputs = [400,600,800,1200,1600,2000,2400,2800,3200]
    md_disc.disc_basic_parameters(readdir,path_out=writedir,overwrite=False,order='<',scl=scl,soutputs=outputs)  ## define disc
#    md_disc.visulize_disc_flows(readdir,path_out=writedir,overwrite=False,order='<',scl=20)
#    md_disc.trace_radial_properties(readdir,path_out=writedir,overwrite=True,order='<',scl=scl)

#############################################################################################
def global_dynamics():
    md_disc.trace_disc_movement(readdir,path_out=writedir,overwrite=True)

def disc_flux():
    md_disc.measure_flux(readdir,path_out=writedir,overwrite=True,order='<',scl=20)
    #md_disc.measure_disc_flux(readdir,path_out=writedir,overwrite=True,order='<',scl=20)
#############################################################################################
def cylinder_disc():
#    md_disc.time_series(readdir,path_out=writedir,overwrite=True,order='<',scl=scl,center_def='baricenter',rotate='AM')
#    md_disc.time_series(readdir,path_out=writedir,overwrite=True,order='<')
    #md_disc.disc_cylinder(readdir,path_out=writedir,overwrite=True,order='<')
    md_disc.calc_disc_properties(readdir,path_out=writedir,overwrite=True,order='<',scl=scl)
#    md_disc.disc_characteristics(readdir,path_out=writedir,overwrite=True,order='<')
    #md_disc.simple_disc_model(writedir)

#############################################################################################
def write_grids():
    md_disc.write_grids(readdir,path_out=writedir,overwrite=True,order='<',scl=20)

#############################################################################################
def detail_ana():
    #md_flux.detail_disc_flux(readdir,path_out=writedir,overwrite=True,order='<',scl=20,Msink_sams=[0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45])
#    md_flux.detail_disc_flux(readdir,path_out=writedir,overwrite=True,order='<',scl=20,Msink_sams=[0.45,0.5,0.525,0.55])
    md_flux.plot_disc_history(writedir)
#############################################################################################
def pretty_plots(outputs=False):
    Disk_Hist = []
    if args.read_dir == 'DC_2_fld_rest2_hres_hsink': iouts = [1612, 1622, 1632, 1642,1652,1662,1672,1682,1692,1702,1712 ]
    if args.read_dir == 'DC_2_fld_rest2_hres_hsink': iouts = [ 1632,1672,1712 ]
    if args.read_dir == 'DC_2_fld_rest2_hres_hsink': iouts = [ 1723, 1733, 1743, 1753 ]
#    if args.read_dir == 'DC_2_fld_rest2_hres_hsink': iouts = [1712 ]
    if args.read_dir == 'DC_2_fld_hres_hsink': iouts = [819,829,839,849,859,869,879,889,899,909,919,929,939,949,959,969,979,989,999,1009,1019]
#    if args.read_dir == 'DC_2_fld_hres_hsink': iouts = [859,929,1019]
    if args.read_dir == 'DC_2_fld_hres_hsink': iouts = [1019]
#    if args.read_dir == 'DC_2_fld_hres_hsink': iouts = [839]
#    if args.read_dir == 'DC_2_fld_hres_hsink': iouts = [839,929,1019]
    if args.read_dir == 'DC_2_fld_res_hsink': iouts = [731,751,771,791,812,832,852,872,892,912,932,952,992,1012,1032,1052,1072,1092,1112,1132,1152,1172]
    if args.read_dir == 'DC_2_fld_res_hsink': iouts = [751,852,952,1052,1152,1172]
    if args.read_dir == 'DC_2_fld': iouts = [122,162,202,242,282,322,362,402,442,482,562,602,642,621,722,842,882,922,962,1002,1042,1082,1122,1162,1202,1242,1282,1322,1362,1402,1442,1482,1522,1562,1602,1642,1682,1722,1762,1802,1842,1882,1922,1962,2002,2042,2082,2122,2162]
#    if args.read_dir == 'DC_2_fld': iouts = iouts[::3]
#    if args.read_dir == 'DC_2_fld': iouts = [121,561,1322,1883,2163]
#    if args.read_dir == 'DC_2_fld': iouts = [1603]
#    if args.read_dir == 'DC_2_fld': iouts = [1723]
    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink': iouts = [1750 ]
    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink': iouts = [1751,1752,1753,1754,1755,1756,1757,1758 ]
    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink': iouts = [1720, 1730,1740,1750,1758]
    if args.read_dir == 'DC_2_rest2': iouts=[1487]
    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink_lowdtAD': iouts=[1770, 1780, 1790]
    if args.read_dir == 'DC_2_fld_rest3_hres_hsink': iouts =  [400,600,800,1200,1600,2000,2400,2800,3200]
#   for time in [0.06, 0.07, 0.08, 0.09 ]:
#       md_disc.some_plots(readdir,path_out=writedir,overwrite=True,order='<',scl=20,time_sam=time)
#    for ms in [0.1,0.2, 0.21, 0.22, 0.23,0.24,0.25,0.26,0.3]: #[0.31,0.32,0.33]:#[0.2, 0.21, 0.22, 0.23,0.24,0.25,0.26]: #[0.1,0.2,0.3]:
#   for ms in [0.31,0.32,0.33,0.34]:
#    for ms in [0.1,0.2, 0.3]: #DC_2_fld
#       history = md_disc.some_plots(readdir,path_out=writedir,overwrite=True,order='<',scl=20,Msink_sam=ms)
#   for ms in [0.39]:
#    for time in [0.0852068922259]:
#    for time in [0.10098]:
#       history = md_disc.some_plots(readdir,path_out=writedir,overwrite=True,order='<',scl=20,time_sam=time)
    #for iout in [1520,1540,1560,1580,1600,1620]:
#    for iout in [729,741]: 
#    for iout in [819,829,839]:  #DC_2_fld_hres_hsink
#    for iout in [1612, 1622, 1632, 1642]:  #DC_2_fld_rest2_hres_hsink
#    for iout in [3396,3398,3400]:

    for iout in outputs:
#    for iout in [802]:
#    for iout in iouts:
       history = md_disc.some_plots(readdir,path_out=writedir,overwrite=True,order='<',scl=20,out_sam=iout)
#       print 
#       plt.show()
#       Disk_Hist.append(history)
#    save_pickle(writedir+'/disk_mr_histoty.npz', [Disk_Hist]) 
    plt.show()
#############################################################################################
def plot_MR_history():
    rdir, wdir0 = local_dir('DC_2_fld')
    path_outs = [wdir0]
    #for Dir in ['DC_2_fld_res_hsink','DC_2_fld_hres_hsink','DC_2_fld_rest2_hres_hsink']:
    #for Dir in ['DC_2_fld_rest2_vhres_hsink','DC_2_fld_rest3_hres_hsink','DC_2_fld_rest2_hres_hsink','DC_2_fld_rest2_vhres_hsink_lowdtAD']:
    for Dir in ['DC_2_fld_rest2_vhres_hsink','DC_2_fld_rest3_hres_hsink','DC_2_fld_nofeed']:
        rdir, wdir = local_dir(Dir)
        path_outs.append(wdir)
    md_disc.plot_MR_history(path_outs) 


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

    print outputs
#    soutputs = [600,800,1200,1600,1622,2000,2400,2800,3200]
#    soutputs = [1602,1612,1622,1632,1642,3378,3388,3398,3408,3418] #DC_2_fld

    #soutputs = [3393,3413,3433]  #rest3
    #soutputs = [1710,1720,1730,1740,1750,1760] # rest2_hres_hsink
#    soutputs = [1764,1774,1784,1794]  # rest2_vhres_hsink
    #soutputs = [1769,1779,1789,1799] # rest2_vhres_hsink_lowdtAD
    soutputs = [1760]
    if args.read_dir == 'DC_2_fld': soutputs = [3398]#, 1622,3398]
    if args.read_dir == 'DC_2_fld': soutputs = [600,800,1200,1600,1622,2000,2400,2800,3200,1602,1612,1622,1632,1642,3378,3388,3398,3408,3418,3500,3600,3700,3800,3900,4000]
    #if args.read_dir == 'DC_2_fld': soutputs = [2000,2400,2800,3200]
    if args.read_dir == 'DC_2_fld_rest3_hres_hsink': soutputs = [3600]
#    if args.read_dir == 'DC_2_fld_rest3_hres_hsink': soutputs = [3400,3450,3500,3550]#,3600]
    #soutputs = [2000,2400,2800,3200]
    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink': soutputs = [1794]
#    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink': soutputs = [1760,1770,1780]#,1794]
    if args.read_dir == 'DC_2_fld_rest2_hres_hsink': soutputs = [1760]
#    if args.read_dir == 'DC_2_fld_rest2_hres_hsink': soutputs = [1600,1620,1640,1660,1680,1700,1720,1740,1760]
    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink_lowdtAD': soutputs = [1798]
    if args.read_dir == 'DC_2_fld_rest2_vhres_hsink_lowdtAD': soutputs = [1770,1780,1790,1798]
    if args.read_dir == 'DC_2_fld': soutputs = [162,202,242,282,322,362,402,442,482,562,602,642,682,722,842,882,922,962,1002,1042,1082,1122,1162,1202,1242,1282,1322,1362,1402,1442,1482,1522,1562,1602,1642,1682,1722,1762,1802,1842,1882,1922,1962,2002,2042,2082,2122,2162]
    #if args.read_dir == 'DC_2_fld': soutputs = [242,282]
    #if args.read_dir == 'DC_2_fld': soutputs = [322,362,402,442,482,562,602,642,682,722,842,882,922,962,1002,1042,1082,1122,1162,1202,1242,1282,1322,1362,1402,1442,1482,1522,1562,1602,1642,1682,1722,1762,1802,1842,1882,1922,1962,2002,2042,2082,2122,2162]
    if args.read_dir == 'DC_2_fld': soutputs = [600,800,1200,1600,2000,2400,2800,3200,3300,3500,3600,3700,3800,3900,4000,140,160]
    if args.read_dir == 'DC_2_fld': soutputs = [3500,3600,3700,3800,3900,4000,140,160]
    if args.read_dir == 'DC_2_fld': soutputs = [1000,1100,1300,1400,1500,1700,1800,1900,2100,2200,2300,2500,2600,2700,2900,3000,3100,3300,3400]
    if args.read_dir == 'DC_2_fld': soutputs = [3400]
    if args.read_dir == 'DC_2_fld_nofeed': soutputs = [2520,2560,2600,2640,2680]
    if args.read_dir == 'DC_2_fld': soutputs = [1620]
    if args.read_dir == 'DC_2_fld': soutputs = [2600]
    if args.read_dir == 'DC_2_fld_nofeed': soutputs = [2600]
#    if args.read_dir == 'DC_2_fld': soutputs = [1602,1632,1642,140,160]
    #soutputs = [3398]
    #soutputs = soutputs[-4::-1]
#    threshold_disc(outputs=soutputs)
#    cylinder_disc()
    #global_dynamics()
    #disc_flux()
    #write_grids()
    #detail_ana()
    #soutputs = soutputs[2:]
    #if args.read_dir == 'DC_2_fld':     soutputs = [3398]
#
    #threshold_disc(outputs=soutputs)

    pretty_plots(outputs=soutputs)

#    plot_MR_history()
    #md_disc.flux_cube(readdir,path_out=writedir,outexist=outputs)

#    nsam = max(nout/args.nsample,1)
#    md_disc.time_series(readdir,path_out=writedir,overwrite=True,order='<',scl=scl,center_def='baricenter',rotate='AM')
#    md_disc.time_series(readdir,path_out=writedir,overwrite=True,order='<')
#    md_disc.disc_cylinder(readdir,path_out=writedir,overwrite=True,order='<')
#    md_disc.simple_disc_model(writedir)
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
