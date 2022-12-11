#! /usr/bin/env python

#This script provides modules for analysis of the core collapse 
    
from pymses import RamsesOutput
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
from pymses.utils import constants as C
from scipy.interpolate import interp1d

from pymses.utils.regions import Sphere, Cylinder
from pymses.filters import CellsToPoints, RegionFilter
import numpy as np
import multiprocessing as mp
import os
##############################################################################
def normalisation(lbox):

    mu  =  2.31
    mp =  mu * 1.660531e-24  #n gramme
    G = 6.7e-8
    kbol  =  1.38062e-16   # erg/degre
    pc=3.08e18 #cm in pc
    Ms = 2.e33 #solar mass in g
    Myr =  3.15e13 # time in Myrs
#scale_n converts particle density from user unit into number of particle per cm^-3
    scale_n = 1. 
#scale_d converts mass density from user unit into g/cc
    scale_d = scale_n*mp  #code density unit is in particle/cc
#scale_t converts time from user unit into seconds
    scale_t = 1.0/np.sqrt(G*scale_d)
#scale_l converts distance from user unit into cm
    scale_l = pc #code distance unit is in parsec
#scale_v convert velocity in user unit into cm/s
    scale_v = scale_l / scale_t
#scale_T2 converts (P/rho) in user unit into T in Kelvin
    scale_T2 = mp/kbol * scale_v**2
#scale energy
    scale_ener = scale_d * scale_v**2
#scale magnetic field
    scale_mag = np.sqrt(scale_ener*4.*np.pi)
    microG = 1.e6
#km/s
    km_s = 1.e5
#Cwnm
    Cwnm=np.sqrt(5./3.*kbol*8000./mp)
    scale_mass = mp*scale_n*(lbox*pc)**3   #mass in g
    unit_col = lbox * pc * scale_n
#lbox en pc
    lbox_pc = lbox #pour cette normalisation....
    return  pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc

#############################################################################################

def renormalize_units(readdir): #return all necessary renormalization units
    outputs = search_ro(readdir)
    ro = RamsesOutput(readdir, outputs[-1])
    boxlen = ro.info["boxlen"]
    lmax = ro.info["levelmax"]
    lmin = ro.info["levelmin"]
    fac_len = boxlen/ro.info["unit_length"].express(C.pc)
    mu = 2.31 #2.31  #ramses2.31, pymses 1/0.76
    fac_mu = 2.31 * 0.76

    ro.info["unit_length"] = ro.info["unit_length"]*fac_len
    ro.info["unit_density"] = ro.info["unit_density"]
    ro.info["unit_mass"] = ro.info["unit_mass"]*fac_len**3
    ro.info["unit_time"] = ro.info["unit_time"]
    ro.info["unit_velocity"] = ro.info["unit_velocity"]*fac_len
    ro.info["unit_pressure"] = ro.info["unit_pressure"]*fac_len**2
    ro.info["unit_mag"] = ro.info["unit_mag"]*fac_len
    ro.info["unit_temperature"] = ro.info["unit_temperature"]*fac_len**2 * mu
    ro.info["unit_mag"] = ro.info["unit_mag"] *fac_len
 
    pressure_P = ro.info["unit_pressure"].express(C.kg/C.m/C.s**2)
    temperature_K = ro.info["unit_temperature"].express(C.K)
    len_pc = ro.info["unit_length"].express(C.pc)
    len_cm = ro.info["unit_length"].express(C.cm)
    vel_ms = ro.info["unit_velocity"].express(C.m/C.s)
    mass_sol = (ro.info["unit_density"]*ro.info["unit_length"]**3).express(C.Msun)
    mass_kg = (ro.info["unit_density"]*ro.info["unit_length"]**3).express(C.kg)
    dens_gcc = ro.info["unit_density"].express(C.g_cc)
    dens_Hcc = ro.info["unit_density"].express(C.H_cc)/fac_mu
    unit_J = (ro.info["unit_density"]*ro.info["unit_length"]**4*ro.info["unit_velocity"]).express(C.kg*C.m**2/C.s)
    time_Myr = ro.info["unit_time"].express(C.Myr)
    mag_gauss = ro.info["unit_mag"].express(C.Gauss)
    unit_Ek = (ro.info["unit_density"]*ro.info["unit_length"]**3*ro.info["unit_velocity"]**2).express(C.J)
    unit_Eg0 = (ro.info["unit_density"]*ro.info["unit_velocity"]/ro.info["unit_time"]*ro.info["unit_length"]**4).express(C.J)
#    unit_Eg1 = (C.G * ro.info['unit_density']**2 * ro.info['unit_length']**5).express(C.J)
    unit_Ep = (ro.info["unit_pressure"]*ro.info["unit_length"]**3).express(C.J)
    unit_Eb = (ro.info["unit_mag"].express(C.T)**2*ro.info["unit_length"].express(C.m)**3)*1.e7/4./np.pi
    return boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb
#############################################################################################
def norm_units(readdir,ioutput=-1):
    outputs = search_ro(readdir)
    if ioutput==-1:  ro = RamsesOutput(readdir, outputs[-1])
    else: ro = RamsesOutput(readdir, ioutput)
    boxlen = ro.info["boxlen"]
    lmax = ro.info["levelmax"]
    lmin = ro.info["levelmin"]
    fac_len = boxlen/ro.info["unit_length"].express(C.pc)
    mu = 2.31 #2.31  #ramses2.31, pymses 1/0.76
    fac_mu = 2.31 * 0.76

    ro.info["unit_length"] = ro.info["unit_length"]*fac_len
    ro.info["unit_density"] = ro.info["unit_density"]
    ro.info["unit_mass"] = ro.info["unit_mass"]*fac_len**3
    ro.info["unit_time"] = ro.info["unit_time"]
    ro.info["unit_velocity"] = ro.info["unit_velocity"]*fac_len
    ro.info["unit_pressure"] = ro.info["unit_pressure"]*fac_len**2
    ro.info["unit_mag"] = ro.info["unit_mag"]*fac_len
    ro.info["unit_temperature"] = ro.info["unit_temperature"]*fac_len**2 * mu
    ro.info["unit_mag"] = ro.info["unit_mag"] *fac_len

    units = dict()
    units["boxlen"] = boxlen
    units["lmax"] = lmax
    units["lmin"] = lmin
    units["len_m"] = ro.info["unit_length"].express(C.m)
    units["len_AU"] = ro.info["unit_length"].express(C.au)
    units["len_pc"] = ro.info["unit_length"].express(C.pc)
    units["time_s"] = ro.info["unit_time"].express(C.s)
    units["time_yr"] = ro.info["unit_time"].express(C.year)
    units["time_kyr"] = units["time_yr"] * 1.e-3
    units["time_Myr"] = units["time_yr"] * 1.e-6
    units["dens_gcc"] = ro.info["unit_density"].express(C.g_cc)
    units["mass_kg"] = (ro.info["unit_density"]*ro.info["unit_length"]**3).express(C.kg)
    units["mass_sol"] = (ro.info["unit_density"]*ro.info["unit_length"]**3).express(C.Msun)
    units["vel_ms"] = ro.info["unit_velocity"].express(C.m/C.s)
    units["J_mvr"] = (ro.info["unit_density"]*ro.info["unit_length"]**4*ro.info["unit_velocity"]).express(C.kg*C.m**2/C.s)
    units["pres_P"] = ro.info["unit_pressure"].express(C.kg/C.m/C.s**2)
    units["temp_K"] = ro.info["unit_temperature"].express(C.K)
    units["mag_gauss"] = ro.info["unit_mag"].express(C.Gauss)
    units["Ek_J"] = (ro.info["unit_density"]*ro.info["unit_length"]**3*ro.info["unit_velocity"]**2).express(C.J)
    units["Eg_J"] = (ro.info["unit_density"]*ro.info["unit_velocity"]/ro.info["unit_time"]*ro.info["unit_length"]**4).express(C.J)
    units["Eb_J"] = (ro.info["unit_mag"].express(C.T)**2*ro.info["unit_length"].express(C.m)**3)*1.e7/4./np.pi
    units["Ep_J"] = (ro.info["unit_pressure"]*ro.info["unit_length"]**3).express(C.J)
#    print boxlen, units["Ek_J"], units["Eg_J"], units["Eb_J"], units["Ep_J"]
    return units
#############################################################################################

def cal_M_R(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.5,scl=20): #calculate mass as function of radius
    amr = ro.amr_source(["rho"])
    sph = Sphere(c_center, c_radius)
    amr = RegionFilter(sph,amr)
    cell_source = CellsToPoints(amr, smallest_cell_level=scl)
    cells = cell_source.flatten()
    dm = cells["rho"]*cells.get_sizes()**3
    r = np.sqrt( np.sum((cells.points-c_center)**2,axis=1) )
    sort = np.argsort(r)
    return interp1d(np.append(0,r[sort]),np.append(0,np.cumsum(dm[sort])))

#############################################################################################

def find_center_dmax(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20): #find the amr cell of the highest density as the center of the cluster
    amr = ro.amr_source(["rho"])
    sph = Sphere(c_center, c_radius)
    amr = RegionFilter(sph,amr)
    cell_source = CellsToPoints(amr, smallest_cell_level=scl)
    cells = cell_source.flatten()
    c_center = cells.points[cells["rho"]==max(cells["rho"])][0]
    return c_center

#############################################################################################

def find_baricenter(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20): #find the amr cell of the highest density as the center of the cluster
    c_center = np.asarray(c_center)
    amr = ro.amr_source(["rho","vel"])
    dis = 1.
    scl2 = scl/2 
    print 'radius', c_radius, 'scl', scl
    #sph = Sphere(c_center, c_radius)
    #amr_sph = RegionFilter(sph,amr)
    while dis>0.5**(scl+1) and scl2<scl+2:
        print 'read level', scl2
        sph = Sphere(c_center, c_radius)
        amr_sph = RegionFilter(sph,amr)
        cell_source = CellsToPoints(amr_sph, smallest_cell_level=scl2)
        cells = cell_source.flatten()
        dm = cells["rho"]*cells.get_sizes()**3
        if len(dm)==0: 
            scl2 += 1
            continue
        dx2 = np.sum( (cells.points-c_center)**2,axis=1 )
        b_center = np.sum((dm)[:,np.newaxis]*cells.points,axis=0)/np.sum(dm)
        dis = np.linalg.norm(b_center-c_center)        
        print 'new baricenter', b_center
        print 'dis =', dis, 0.5**(scl2+1),0.5**(scl+1)
        scl2 = np.minimum(np.ceil(np.log(dis)/np.log(0.5)),scl)
        c_center = b_center
    return c_center

#############################################################################################

def find_axis(ro, c_center=np.array([0.5,0.5,0.5]),c_radius=0.05,scl=20): #find the axis of rotation
    c_center = np.asarray(c_center)
    amr = ro.amr_source(["rho","vel"])
    sph = Sphere(c_center, c_radius)
    amr = RegionFilter(sph,amr)
    cell_source = CellsToPoints(amr, smallest_cell_level=scl)
    cells = cell_source.flatten()
    dm = cells["rho"]*cells.get_sizes()**3
    c_vel = np.sum(dm[:,np.newaxis]*cells["vel"],axis=0)/np.sum(dm)
    c_axis = np.cross((cells.points-c_center),(cells["vel"]-c_vel))
    c_axis = np.sum(dm[:,np.newaxis]*c_axis,axis=0)
    c_axis = c_axis / np.linalg.norm(c_axis)
    print 'rotation axis', c_axis
    return c_axis, c_vel

#############################################################################################
## rotation matrix for given axis
def cal_rot_mat(c_axis):
    c_axis = c_axis/np.linalg.norm(c_axis)
    if c_axis[2]<1.:
        cos_th = c_axis[2]
        sin_th = np.sqrt(1.-cos_th**2)
        cos_ph = c_axis[0]/sin_th
        sin_ph = c_axis[1]/sin_th
    else: cos_th=1.; sin_th=0.; cos_ph=1.; sin_ph=0.
    rot_x = np.array([cos_th*cos_ph, cos_th*sin_ph, -sin_ph])
    rot_y = np.array([-sin_ph, cos_ph, 0])
    mat_Rot = np.vstack((rot_x, rot_y, c_axis))
    return mat_Rot

##############################################################################
## read sink particle properties from the csv file
def read_sink_cvs(num,directory,deli=''):
    name = directory+'/output_'+str(num).zfill(5)+'/sink_'+str(num).zfill(5)+'.csv'
    print 'read sinks from ', name
    sinks = np.loadtxt(name,ndmin=2,delimiter=deli)#,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12))
    #if len(sinks)==0: sinks = np.empty([0, 0]) 
    #elif len(sinks.shape)==1: sinks=sinks[np.newaxis,:]    
    return sinks
########################################################################################
##make images
##amr and ro pymses stuff
##center 3D array for coordinates center
##radius the size of the region to be plotted (between 0 and 1)
##num output number
##directory : the directory in which the images are being written 
##x_v, y_v and z_v coordinates of the points to be overplotted
##i_im image number
##force test whether the first of the image series exists already. If yes the return
##make_phi draw gravitational potential
########################################################################################


                                                                                               
##############################################################################
### do a series of zoom
### path : the path
### num the output number
### zoom_v an arrays which contains the zoom
### sinks particle can be overplotted and used to place the zoom
##############################################################################
def make_image_zoom(path,num,zoom_v,path_out=None,sinks=False,force=True,center_img=[0.5,0.5,0.5],make_phi=False,deli=',',col_dens_only=False,tag='',gcomp=True,mag_im=True,vel_red=20,map_size=512,im_pos=0.,conv_sink=1,center_def='None',order='<',ps=False,ind_sink=0,all_sink=True,rotate=False):
    ro=pymses.RamsesOutput(path,num,order=order)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)
    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)
    if(not make_phi): amr = ro.amr_source(["rho","vel","P","Bl","Br"])
    else: amr = ro.amr_source(["rho","vel","P","phi","g"], grav_compat=gcomp)
    if center_def=='dmax': center_img = find_center_dmax(ro, c_center=center_img,c_radius=0.005,scl=20)  
    elif center_def=='baricenter': center_img = find_baricenter(ro, c_center=center_img,c_radius=0.05,scl=20)

    mass_v = None
    pos_v = [None, None, None]
    if sinks: #read sink file
        sinks = read_sink_cvs(num,path,deli=deli)
        if len(sinks)>0:
            mass_v = sinks[:,1] * scale_mass / Ms / lbox**3
            pos_v = sinks[:,1+conv_sink:4+cin_sink]
#position of the most massive sink
#            ind_max = np.argmax(sinks[:,1])
            args = np.argsort(sinks[:,1])[::-1]
            ind_sink = np.min((ind_sink,len(sinks)))
            ind_max = args[ind_sink]
            if center_def=='sink':
                center_img=pos_v[ind_max,:]/lbox
                tag= 'ns'+str(ind_sink)+'_'
    for i_im in ragne(len(zoom_v)):
        rad = zoom_v[i_im]
        center = center_img.copy()
        center[center<rad]=rad; center[center>(1-rad)]=1-rad #re-center big box
        print 'center ', center, 'rad ', rad

        if make_phi: make_image_phi(amr,ro,center,rad,num,path,x_v=x_v,y_v=y_v,z_v=z_v,i_im=i_im,force=force,path_out=path_out)
        else:
            make_image(amr,ro,center,rad,num,path,x_v=x_v,y_v=y_v,z_v=z_v,mass_v=mass_v,i_im=i_im,force=force,path_out=path_out,col_dens_only=col_dens_only,tag=tag,vel_red=vel_red,map_size=map_size,im_pos=im_pos,ps=ps,all_sink=all_sink)
            if (not col_dens_only and mag_im):  make_image_B(amr,ro,center,rad,num,path,x_v=None,y_v=None,z_v=None,mass_v=None,i_im=i_im,force=force,path_out=path_out,tag=tag,ps=ps)
    return
##############################################################################

