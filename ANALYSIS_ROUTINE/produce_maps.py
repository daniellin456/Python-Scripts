#! /usr/bin/env python
# This module contains functions for imaging the ramses output
import matplotlib
matplotlib.use('Agg')
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
#from module_core import renormalize_units, find_center_dmax, find_baricenter, find_axis, normalisation, read_sink_cvs, cal_M_R
#from module_analysis import cylinder_flux, save_pickle, read_pickle
#############################################################################################


def local_dir(dir):
    #ldirs = ['/gpfs/users/', '/drf/projets/alfven-data/','/dsm/anais/storageA/','/drf/projets/capucine/']
    ldirs = ['/data/']
    for ldir in ldirs:
        readdir = dir
        writedir = dir + '/figures'
        #if ldir == ldirs[3]:  writedir = ldir[2]+'ylee/'+dir
        if os.path.exists(readdir):
            if not os.path.exists(writedir):
                os.makedirs(writedir)
            return readdir, writedir
        else:
            if os.path.exists(writedir):
                readdir = writedir
                return readdir, writedir
    return dir, dir

###########################################################################


def save_fig(figname, ps=False, tight=False):
    if tight:
        plt.tight_layout()
    plt.savefig(figname+'.pdf')
    if ps:
        plt.savefig(figname+'.eps')
#############################################################################################


def norm_units(readdir, ioutput=-1):
    outputs = search_ro(readdir)
    if ioutput == -1:
        ro = RamsesOutput(readdir, outputs[-1])
    else:
        ro = RamsesOutput(readdir, ioutput)
    boxlen = ro.info["boxlen"]
    lmax = ro.info["levelmax"]
    lmin = ro.info["levelmin"]
    fac_len = boxlen/ro.info["unit_length"].express(C.pc)
    mu = 2.31  # 2.31  #ramses2.31, pymses 1/0.76
    fac_mu = 2.31 * 0.76
    print fac_len

    ro.info["unit_length"] = ro.info["unit_length"]*fac_len
    ro.info["unit_density"] = ro.info["unit_density"]
    ro.info["unit_mass"] = ro.info["unit_mass"]*fac_len**3
    ro.info["unit_time"] = ro.info["unit_time"]
    ro.info["unit_velocity"] = ro.info["unit_velocity"]*fac_len
    ro.info["unit_pressure"] = ro.info["unit_pressure"]*fac_len**2
    ro.info["unit_mag"] = ro.info["unit_mag"]*fac_len
    ro.info["unit_temperature"] = ro.info["unit_temperature"]*fac_len**2 * mu
    ro.info["unit_mag"] = ro.info["unit_mag"] * fac_len

    units = dict()
    units["boxlen"] = boxlen
    units["lmax"] = lmax
    units["lmin"] = lmin
    units["len_m"] = ro.info["unit_length"].express(C.m)
    units["len_cm"] = ro.info["unit_length"].express(C.cm)
    units["len_AU"] = ro.info["unit_length"].express(C.au)
    units["len_pc"] = ro.info["unit_length"].express(C.pc)
    units["time_s"] = ro.info["unit_time"].express(C.s)
    units["time_yr"] = ro.info["unit_time"].express(C.year)
    units["time_kyr"] = units["time_yr"] * 1.e-3
    units["time_Myr"] = units["time_yr"] * 1.e-6
    units["dens_gcc"] = ro.info["unit_density"].express(C.g_cc)
    units["mass_kg"] = (ro.info["unit_density"] *
                        ro.info["unit_length"]**3).express(C.kg)
    units["mass_sol"] = (ro.info["unit_density"] *
                         ro.info["unit_length"]**3).express(C.Msun)
    units["vel_ms"] = ro.info["unit_velocity"].express(C.m/C.s)
    units["J_mvr"] = (ro.info["unit_density"]*ro.info["unit_length"]
                      ** 4*ro.info["unit_velocity"]).express(C.kg*C.m**2/C.s)
    units["pres_P"] = ro.info["unit_pressure"].express(C.kg/C.m/C.s**2)
    units["temp_K"] = ro.info["unit_temperature"].express(C.K)
    units["mag_gauss"] = ro.info["unit_mag"].express(C.Gauss)
    return units
###########################################################################
# produce maps for various properties


def colden_map(ro, amr, camera):  # column density map
    rho_op = ScalarOperator(lambda dset: dset["rho"], ro.info["unit_density"])
    rt = raytracing.RayTracer(amr, ro.info, rho_op)
    map_colden = rt.process(camera, surf_qty=True)
    map_colden = map_colden.map
    return map_colden.T


def rho_map(ro, amr, camera, dis=0.):  # density slice
    rho_op = ScalarOperator(lambda dset: dset["rho"], ro.info["unit_density"])
    dmap_rho = slicing.SliceMap(amr, camera, rho_op, z=dis)  # ,verbose=False)
    map_rho = dmap_rho.map
    return map_rho.T


def vel_map(ro, amr, camera, ax1, ax2, dis=0.):  # pos velocity slice map
    vx_op = ScalarOperator(lambda dset: np.sum(
        dset["vel"]*ax1[np.newaxis, :], axis=1), ro.info["unit_velocity"])
    dmap_vx = slicing.SliceMap(amr, camera, vx_op, z=dis)  # ,verbose=False)
    map_vx = dmap_vx.map
    vy_op = ScalarOperator(lambda dset: np.sum(
        dset["vel"]*ax2[np.newaxis, :], axis=1), ro.info["unit_velocity"])
    dmap_vy = slicing.SliceMap(amr, camera, vy_op, z=dis)  # ,verbose=False)
    map_vy = dmap_vy.map
    return map_vx.T, map_vy.T


def velz_map(ro, amr, camera, ax, dis=0.):  # los velocity slice map
    vz_op = ScalarOperator(lambda dset: np.sum(
        dset["vel"]*ax[np.newaxis, :], axis=1), ro.info["unit_velocity"])
    dmap_vz = slicing.SliceMap(amr, camera, vz_op, z=dis)
    map_vz = dmap_vz.map
    return map_vz.T


def P_map(ro, amr, camera):  # pressure slice map
    P_op = ScalarOperator(lambda dset: dset["P"],  ro.info["unit_pressure"])
    dmap_P = slicing.SliceMap(amr, camera, P_op, z=0.)  # ,verbose=False)
    map_P = dmap_P.map
    return map_P.T


def T_map(ro, amr, camera):  # temperature slice map
    T_op = ScalarOperator(lambda dset: dset["T"],  ro.info["unit_temperature"])
    dmap_T = slicing.SliceMap(amr, camera, T_op, z=0.)  # ,verbose=False)
    map_T = dmap_T.map
    return map_T.T


# B rho integrated map, divide by column density for mass-weighted B
def Bxy_colden_map(ro, amr, camera, ax1, ax2):
    cells = CellsToPoints(amr).flatten()
    B_op = ScalarOperator(lambda dset: 0.5*(np.sum(
        (dset["Br"]+dset["Bl"])*ax1[np.newaxis, :], axis=-1))*dset["rho"], ro.info["unit_mag"])
    rt = raytracing.RayTracer(amr, ro.info, B_op)
    map_Bcolden = rt.process(camera, surf_qty=True)
    map_Bxcolden = map_Bcolden.map
    B_op = ScalarOperator(lambda dset: 0.5*(np.sum(
        (dset["Br"]+dset["Bl"])*ax2[np.newaxis, :], axis=-1))*dset["rho"], ro.info["unit_mag"])
    rt = raytracing.RayTracer(amr, ro.info, B_op)
    map_Bcolden = rt.process(camera, surf_qty=True)
    map_Bycolden = map_Bcolden.map
    return map_Bxcolden.T, map_Bycolden.T


def Bxy_slice_map(ro, amr, camera, ax1, ax2, dis=0.):  # B slice map
    B_op = ScalarOperator(
        lambda dset: 0.5*(np.sum((dset["Br"]+dset["Bl"])*ax1[np.newaxis, :], axis=-1)), ro.info["unit_mag"])
    dmap_B = slicing.SliceMap(amr, camera, B_op, z=dis)  # ,verbose=False)
    map_Bx = dmap_B.map
    B_op = ScalarOperator(
        lambda dset: 0.5*(np.sum((dset["Br"]+dset["Bl"])*ax2[np.newaxis, :], axis=-1)), ro.info["unit_mag"])
    dmap_B = slicing.SliceMap(amr, camera, B_op, z=dis)  # ,verbose=False)
    map_By = dmap_B.map
    return map_Bx.T, map_By.T


def Pscal_map(ro, amr, camera, i):  # passive scalar map
    Pscal_op = ScalarOperator(
        lambda dset: dset["rho"]*dset["passive"][:, :, i], ro.info["unit_density"])
    rt = raytracing.RayTracer(amr, ro.info, Pscal_op)
    dmap_Pscal = rt.process(camera, surf_qty=True)
    map_Pscal = dmap_Pscal.map
    return map_Pscal.T


def MaxLev_map(ro, amr, camera, read=30):  # max level map
    level_op = MaxLevelOperator()
    amr.set_read_levelmax(read)
    rt = raytracing.RayTracer(amr, ro.info, level_op)
    datamap = rt.process(camera, surf_qty=True)
    map_level = datamap.map
    return map_level.T


def cpu_map(ro, amr, camera):  # cpu map
    cpu_op = ScalarOperator(lambda dset: dset.icpu *
                            (np.ones_like(dset["P"])), ro.info["unit_pressure"])
    rt = raytracing.RayTracer(amr, ro.info, cpu_op)
    datamap = rt.process(camera, surf_qty=True)
    map_cpu = datamap.map
    return map_cpu.T


###########################################################################
def define_cameras(center=[0.5, 0.5, 0.5], radius=0.5, map_size=512, ax1=None, ax2=None, ax3=None):
    '''define 3 cameras with given center, size, resolution, orientation, axes must be perpendicular if asigned '''
    if ax1 is not None:
        ax_z = np.asarray(ax1)/np.linalg.norm(ax1)
        if ax2 is not None:
            ax_x = np.asarray(ax2)
        else:
            if ax_z[2] >= 1:
                ax_x = np.array([1, 0, 0])
            else:
                cos_sin = ax_z[2]/np.sqrt(ax_z[0]**2+ax_z[1]**2)
                ax_x = ax_z.copy()
                ax_x[0:2] *= -cos_sin
                ax_x[2] /= cos_sin
        ax_x = ax_x / np.linalg.norm(ax_x)
        if ax3 is not None:
            ax_y = ax3
        else:
            ax_y = np.cross(ax_z, ax_x)
    else:
        ax_x = 'x'
        ax_y = 'y'
        ax_z = 'z'
    cam_x = Camera(center=center, line_of_sight_axis=ax_x, region_size=[
                   2.*radius, 2.*radius], distance=radius, far_cut_depth=radius, up_vector=ax_z, map_max_size=map_size)
    cam_y = Camera(center=center, line_of_sight_axis=ax_y, region_size=[
                   2.*radius, 2.*radius], distance=radius, far_cut_depth=radius, up_vector=ax_x, map_max_size=map_size)
    cam_z = Camera(center=center, line_of_sight_axis=ax_z, region_size=[
                   2.*radius, 2.*radius], distance=radius, far_cut_depth=radius, up_vector=ax_y, map_max_size=map_size)
    return [cam_x, cam_y, cam_z]

###########################################################################


def make_image(amr, ro, center, radius, num, path, pos_sink=None, mass_sink=None, i_im=0, overwrite=False, path_out=None, col_dens_only=False, map_size=20, vel_size=20, tag='', ps=False, all_sink=True, ax1=None, ax2=None, PSCAL=[], add_B=False, col_dens = False, histogram=False, pressure=False, temperature=False, rho_and_vel=False, amr_level=False):
    # ax1: new z axis, ax2: new x axis
    if (path_out is not None):
        directory = path_out
    else:
        directory = path

    #boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    units = norm_units(path, ioutput=1)

    len_AU = units["len_AU"]
    time_Myr = units["time_Myr"]
    dens_gcc = units["dens_gcc"]
    len_cm = units["len_cm"]
    vel_ms = units["vel_ms"]
    pressure_P = units["pres_P"]
    mag_gauss = units["mag_gauss"]
    lmax = units["lmax"]
    dx = 0.5**lmax
    print 'len_AU', len_AU
    time = ro.info['time'] * time_Myr
    ind = path[::-1].find('/')
    if ind > 0:
        titlepath = path[len(path)-ind:]
    else:
        titlepath = path
    titre = 'Time = %.3f kyr' % (time*1.e3)

   # name = directory+'/coldens_z'+'_'+str(i_im)+'_'+format(num,'05')+'.pdf'
   # if (os.path.exists(name) and not overwrite): return

    loss = ['x', 'y', 'z']
    x_axis = np.array([1, 0, 0])
    y_axis = np.array([0, 1, 0])
    z_axis = np.array([0, 0, 1])
    axes = np.vstack((x_axis, y_axis, z_axis))

    cameras = define_cameras(center=center, radius=radius,
                             map_size=map_size, ax1=z_axis, ax2=x_axis)
    cameras_vel = define_cameras(
        center=center, radius=radius, map_size=vel_size, ax1=z_axis, ax2=x_axis)
    ind_selects = [[1, 2, 0], [2, 0, 1], [0, 1, 2]]

    if len(PSCAL) > 0:
        if PSCAL[-1] == 100:
            PSCAL_only = True
            PSCAL = PSCAL[:-1]
        else:
            PSCAL_only = False
    for los, cam, cam_v, ind in zip(loss, cameras, cameras_vel, ind_selects):
        print 'center', center, 'ind', ind
        cen1 = center[ind[0]]
        cen2 = center[ind[1]]
        cen3 = center[ind[2]]
        fig_lim = [(cen1-radius-0.5)*len_AU, (cen1+radius-0.5)*len_AU,
                   (cen2-radius-0.5)*len_AU, (cen2+radius-0.5)*len_AU]
        #fig_lim = [-radius*len_AU,radius*len_AU,-radius*len_AU,radius*len_AU]
        if len(PSCAL) > 0:
            if PSCAL[0] == -1:  # plot all passive scalar on same figure
                figname_sum = directory+'/sum_pscal_'+los + \
                    '_'+tag+str(i_im)+'_'+format(num, '05')
                if not os.path.exists(figname_sum+'.pdf') or overwrite:
                    comp = ['passive_'+str(p) for p in PSCAL[1:]]
                    for i in PSCAL[1:]:
                        figname = directory+'/' + \
                            comp[i]+'_'+los+'_'+tag + \
                            str(i_im)+'_'+format(num, '05')
                        map_Pscal = Pscal_map(ro, amr, cam, i)
                        if i == PSCAL[1]:
                            map_Pscal_sum = map_Pscal
                        else:
                            map_Pscal_sum += map_Pscal
                        plt.clf()
                        print 'sum pscal', np.sum(map_Pscal)
                        im = plt.imshow(np.log10(map_Pscal*dens_gcc*len_cm)+1., extent=fig_lim,
                                        origin='lower', interpolation='none', zorder=1, vmin=-4, vmax=5)
                        plt.title(titre)
                        plt.xlabel(loss[ind[0]]+' (AU)')
                        plt.ylabel(loss[ind[1]]+' (AU)')
                        cbar = plt.colorbar(im, extend='both')
                        cbar.set_label(
                            r'$\mathrm{log}_{10}\rho (\mathrm{kg}\/\mathrm{m}^{-2})$')
                        save_fig(figname, ps=ps)
                    plt.clf()
                    im = plt.imshow(np.log10(map_Pscal_sum*dens_gcc*len_cm)+1., extent=fig_lim,
                                    origin='lower', interpolation='none', zorder=1, vmin=-4, vmax=5)
                    plt.scatter((cen1-0.5)*len_AU, (cen2-0.5) *
                                len_AU, color='y', marker='*', zorder=2)
                    plt.title(titre)
                    plt.xlabel(loss[ind[0]]+' (AU)')
                    plt.ylabel(loss[ind[1]]+' (AU)')
                    cbar = plt.colorbar(im, extend='both')
                    cbar.set_label(
                        r'$\mathrm{log}_{10}\rho (\mathrm{kg}\/\mathrm{m}^{-2})$')
                    save_fig(directory+'/sum_pscal_'+los+'_'+tag +
                             str(i_im)+'_'+format(num, '05'), ps=ps)

            else:
                comp = ['passive_'+str(p) for p in PSCAL]
                for i in range(len(PSCAL)):
                    figname = directory+'/' + \
                        comp[i]+'_'+los+'_'+tag+str(i_im)+'_'+format(num, '05')
                    if not os.path.exists(figname+'.pdf') or overwrite:
                        #if len(PSCAL)==3: comp = ['dust','gas','CAI']
                        #if len(PSCAL)==2: comp = ['1500K','1650K']
                        plt.clf()
                        map_Pscal = Pscal_map(ro, amr, cam, i)
                        print 'sum pscal', np.sum(map_Pscal)
                        im = plt.imshow(np.log10(map_Pscal*dens_gcc*len_cm)+1., extent=fig_lim,
                                        origin='lower', interpolation='none', zorder=1, vmin=-4, vmax=2)
                        plt.title(titre)
                        plt.xlabel(loss[ind[0]]+' (AU)')
                        plt.ylabel(loss[ind[1]]+' (AU)')
                        cbar = plt.colorbar(im, extend='both')
                        cbar.set_label(
                            r'$\mathrm{log}_{10}\rho (\mathrm{kg}\/\mathrm{m}^{-2})$')
                        save_fig(figname, ps=ps)
            if PSCAL_only:
                print "plot only pscal, continue!"
            if PSCAL_only:
                continue

        if col_dens:
            figname = directory+'/coldens_'+los+'_' + \
                tag+str(i_im)+'_'+format(num, '05')
            if not os.path.exists(figname+'.pdf') or overwrite:
                map_colden = colden_map(ro, amr, cam)  # column density map
                plt.clf()
                im = plt.imshow(np.log10(map_colden*dens_gcc*len_cm)+1., extent=fig_lim,
                                origin='lower', interpolation='none', vmin=-3, vmax=13, cmap='plasma')
                plt.title(titre)
                plt.xlabel(loss[ind[0]]+' (AU)')
                plt.ylabel(loss[ind[1]]+' (AU)')
                cbar = plt.colorbar(im, extend='both')
                cbar.set_label(
                    r'$\mathrm{log}_{10}\Sigma (\mathrm{kg}\/\mathrm{m}^{-2})$')
                save_fig(figname, ps=ps)

                if mass_sink is not None:  # add sinks
                    mask = np.where((abs(pos_sink[:, ind[0]]-cen1) <= radius)
                                    & (abs(pos_sink[:, ind[1]]-cen2) <= radius) == True)
                    if not all_sink:
                        mask = np.where((abs(pos_sink[:, ind[0]]-cen1) <= radius) & (abs(
                            pos_sink[:, ind[1]]-cen2) <= radius) & (abs(pos_sink[:, ind[2]]-cen3) <= radius) == True)
                    siz_sym = np.sqrt(mass_sink[mask])*50.
                    plt.scatter((pos_sink[mask, ind[0]]-0.5)*len_AU, (pos_sink[mask,
                                ind[1]]-0.5)*len_AU, marker='o', s=siz_sym, color='r', zorder=2)
                    circle = plt.Circle(((cen1-0.5)*len_AU, (cen2-0.5)*len_AU), dx *
                                        len_AU*4, fill=False, edgecolor='w', linestyle=':', zorder=2)
                    ax = plt.gca()
                    ax.add_artist(circle)
                    print 'SINK POSITION', (cen1-0.5)*len_AU, (cen2-0.5)*len_AU, dx*len_AU*10
                    plt.xlim(fig_lim[:2])
                    plt.ylim(fig_lim[2:])
                    # plt.show()
                    save_fig(figname+'_sink', ps=ps)

        #if col_dens_only:
        #    continue

        if rho_and_vel:
            figname = directory+'/rho_'+los+'_' + \
                tag+str(i_im)+'_'+format(num, '05')
            if not os.path.exists(figname+'.pdf') or overwrite:
                cen1 = center[ind[0]]
                cen2 = center[ind[1]]
                cen3 = center[ind[2]]
                map_rho = rho_map(ro, amr, cam)  # density slice map + velocity
                map_vx, map_vy = vel_map(
                    ro, amr, cam_v, axes[ind[0], :], axes[ind[1], :])
                plt.clf()
                pos_vel = np.linspace(-radius, radius, vel_size+1)
                pos_vel = (pos_vel[:-1] + radius/vel_size) * len_AU
                xx, yy = np.meshgrid(pos_vel+(cen1-0.5)*len_AU,
                                     pos_vel+(cen2-0.5)*len_AU)
                q = plt.quiver(xx, yy, map_vx*vel_ms, map_vy*vel_ms,
                               color='w', pivot='middle', zorder=2)
                qk = plt.quiverkey(q, 0.8, -0.08, 1.e4,
                                   '10 km/s', color='k', labelpos='W')
                im = plt.imshow(np.log10(map_rho*dens_gcc), extent=fig_lim, origin='lower',
                                interpolation='none', zorder=1, cmap='plasma', vmin=-25, vmax=-13)
                plt.title(titre)
                plt.xlabel(loss[ind[0]]+' (AU)')
                plt.ylabel(loss[ind[1]]+' (AU)')
                cbar = plt.colorbar(im, extend='both')
                cbar.set_label(
                    r'$\mathrm{log}_{10}\rho (\mathrm{g}\/\mathrm{cm}^{-3})$')
                save_fig(figname, ps=ps)

        if add_B:
            figname_slice = directory+'/B_slice_'+los + \
                '_'+tag+str(i_im)+'_'+format(num, '05')
            figname_int = directory+'/B_int_'+los + \
                '_'+tag+str(i_im)+'_'+format(num, '05')
            if not os.path.exists(figname+'.pdf') or overwrite:
                map_Bx_slice, map_By_slice = Bxy_slice_map(
                    ro, amr, cam, axes[ind[0], :], axes[ind[1], :])
                map_Bx_colden, map_By_colden = Bxy_colden_map(
                    ro, amr, cam, axes[ind[0], :], axes[ind[1], :])
                map_Bx_int = map_Bx_colden/map_colden
                map_By_int = map_By_colden/map_colden
                nvec = 20
                B_size = map_size / nvec
                plt.clf()
                pos_B = np.linspace(-radius, radius, B_size+1)
                pos_B = (pos_B[:-1] + radius/B_size) * len_AU
                xx, yy = np.meshgrid(
                    pos_B+(cen1-0.5)*len_AU, pos_B+(cen2-0.5)*len_AU)
                print xx.shape, map_Bx_slice[nvec/2::nvec, nvec/2::nvec].shape, 'shapes!'
                q = plt.quiver(xx, yy, map_Bx_slice[nvec/2::nvec, nvec/2::nvec]*mag_gauss,
                               map_By_slice[nvec/2::nvec, nvec/2::nvec]*mag_gauss, color='w', pivot='middle', zorder=2)
           # qk = plt.quiverkey(q, 0.8, -0.08, 1.e4, '10 km/s',color='k', labelpos='W')
                im = plt.imshow(np.log10(np.sqrt(map_Bx_slice**2+map_By_slice**2)*mag_gauss), extent=fig_lim,
                                origin='lower', interpolation='none', zorder=1, cmap='summer', vmin=-3, vmax=0)
                plt.title(titre)
                plt.xlabel(loss[ind[0]]+' (AU)')
                plt.ylabel(loss[ind[1]]+' (AU)')
                cbar = plt.colorbar(im, extend='both')
                cbar.set_label(r'$\mathrm{log}_{10}B (\mathrm{G})$')
                save_fig(figname_slice, ps=ps)

                plt.clf()
                q = plt.quiver(xx, yy, map_Bx_int[nvec/2::nvec, nvec/2::nvec]*mag_gauss,
                               map_By_int[nvec/2::nvec, nvec/2::nvec]*mag_gauss, color='w', pivot='middle', zorder=2)
           # qk = plt.quiverkey(q, 0.8, -0.08, 1.e4, '10 km/s',color='k', labelpos='W')
                im = plt.imshow(np.log10(np.sqrt(map_Bx_int**2+map_By_int**2)*mag_gauss), extent=fig_lim,
                                origin='lower', interpolation='none', zorder=1, cmap='summer', vmin=-3, vmax=0)
                plt.title(titre)
                plt.xlabel(loss[ind[0]]+' (AU)')
                plt.ylabel(loss[ind[1]]+' (AU)')
                cbar = plt.colorbar(im, extend='both')
                cbar.set_label(r'$\mathrm{log}_{10}B (\mathrm{G})$')
                save_fig(figname_int, ps=ps)

        if pressure:
            plt.clf()
            figname = directory+'/P_'+los+'_' + \
                tag+str(i_im)+'_'+format(num, '05')
            if not os.path.exists(figname+'.pdf') or overwrite:
                map_P = P_map(ro, amr, cam)
                im = plt.imshow(np.log10(map_P*pressure_P), extent=fig_lim, origin='lower',
                                cmap='YlGnBu', interpolation='none', zorder=1, vmin=-16, vmax=-10)
                plt.title(titre)
                plt.xlabel(loss[ind[0]]+' (AU)')
                plt.ylabel(loss[ind[1]]+' (AU)')
                cbar = plt.colorbar(im, extend='both')
                cbar.set_label(r'$\mathrm{log}_{10}P$ (pascal)')
                # plt.show()
                save_fig(figname, ps=ps)

        if temperature:
            plt.clf()
            figname = directory+'/T_'+los+'_' + \
                tag+str(i_im)+'_'+format(num, '05')
            if not os.path.exists(figname+'.pdf') or overwrite:
                map_T = T_map(ro, amr, cam)
                im = plt.imshow(np.log10(map_T), extent=fig_lim, origin='lower',
                                cmap='BuPu_r', interpolation='none', zorder=1, vmin=1, vmax=5)
                plt.title(titre)
                plt.xlabel(loss[ind[0]]+' (AU)')
                plt.ylabel(loss[ind[1]]+' (AU)')
                cbar = plt.colorbar(im, extend='both')
                cbar.set_label(r'$\mathrm{log}_{10}T$ (K)')
                save_fig(figname, ps=ps)

        if histogram:
            plt.clf()
            figname = directory+'/2dHist_'+los+'_' + \
                tag+str(i_im)+'_'+format(num, '05')
            if not os.path.exists(figname+'.pdf') or overwrite:
                map_T = T_map(ro, amr, cam)
                map_rho = rho_map(ro, amr, cam)
                
                plt.title(r"$\rho \; v.s \; T$")
                plt.xlabel(r"$log_{10}\rho$")
                plt.ylabel(r"$log_{10} T$")
                H, xedges, yedges = np.histogram2d(
                    np.log10(map_rho*dens_gcc).flatten(), np.log10(map_T).flatten(), bins=20)
                plt.imshow(H, origin='lower', cmap='BuGn', extent=[
                           xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm())
                plt.colorbar()
                save_fig(directory+'/2dHist_'+los+'_'+tag +
                         str(i_im)+'_'+format(num, '05'), ps=ps)

        if amr_level:
            plt.clf()
            figname = directory+'/amr_'+los+'_' + tag+str(i_im)+'_'+format(num, '05')
            if not os.path.exists(figname+'.pdf') or overwrite:
            
                cmap = plt.get_cmap('rainbow', 7)
                #norm = BoundaryNorm(np.arange(0.5, 16, 1), cmap.N)
                map_level = MaxLev_map(ro,amr,cam,read=30)
                im = plt.imshow( map_level, extent=fig_lim, origin='lower',
                                cmap=cmap, interpolation='none',  vmin=7.5, vmax=14.5)
                bounds = np.linspace(8, 14, 7)
                ticks = np.linspace(8, 14, 7)
                plt.title(titre)
                plt.xlabel(loss[ind[0]]+' (AU)')
                plt.ylabel(loss[ind[1]]+' (AU)')
                cbar = plt.colorbar(im, ticks=np.arange(8, 15) )
                #cbar = plt.colorbar(im,  spacing='proportional', ticks=ticks, boundaries=bounds)
                cbar.set_label('AMR level')
                save_fig(figname, ps=ps)

#        map_P = P_map(ro,amr,cam)
#        map_level = MaxLev_map(ro,amr,cam,read=30) ## max refinement level
#        map_cpu = cpu_map(ro,amr,cam) ## cpu domain map
##############################################################################
# do a series of zoom
# path : the path
# num the output number
# zoom_v an arrays which contains the zoom
# sinks particle can be overplotted and used to place the zoom
##############################################################################


def make_image_zoom(path, num, zoom_v, i_ims, path_out=None, plot_sinks=False, overwrite=True, center_img=[0.5, 0.5, 0.5], make_phi=False, deli=',', col_dens_only=False, tag='', gcomp=True, map_size=512, vel_size=20, xyz_sink=3, center_def='None', order='<', ps=False, ind_sink=0, all_sink=True, rotate=None, make_slice=True, col_dens=True, c_axis=None, PSCAL=[], add_B=False, histogram=False, pressure=False, temperature=False, rho_and_vel=False, amr_level=False):

    ro = pymses.RamsesOutput(path, num, order=order)
#    boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units(path)
    units = norm_units(path, ioutput=1)
    len_AU = units["len_AU"]
    lmax = units["lmax"]
    dx = 0.5**lmax
    vars = ["rho", "vel", "P", "T"]
    if make_phi:
        vars = vars + ["phi", "g"]
    if len(PSCAL) > 0:
        vars = vars + ["passive"]
    if add_B:
        vars = vars + ["Br", "Bl"]
    if not make_phi:
        amr = ro.amr_source(vars)
    else:
        amr = ro.amr_source(vars, grav_compat=gcomp)

    mass_sink = None
    pos_sink = None
    if plot_sinks:  # read sink file
        sinks = read_sink_cvs(num, path, deli=deli)
        if len(sinks) > 0:
            mass_sink = sinks[:, 1] * mass_sol
            pos_sink = sinks[:, xyz_sink:xyz_sink+3]/boxlen
#            ind_max = np.argmax(sinks[:,1])
            args = np.argsort(mass_sink)[::-1]
            ind_sink = np.min((ind_sink, len(sinks)))
            ind_max = args[ind_sink]  # position of the most massive sink
            if center_def == 'sink':
                center_img = pos_sink[ind_max, :]
                tag = 'ns'+str(ind_sink)+'_'
    # for i_im in range(len(zoom_v))[1:2]:
    for i_im, rad in zip(i_ims, zoom_v):
        #        rad = zoom_v[i_im]
        center = center_img.copy()
        center[center < rad] = rad
        center[center > (1-rad)] = 1-rad  # re-center big box
        print 'center ', center, 'rad ', rad

        #if col_dens:
        make_image(amr, ro, center, rad, num, path, pos_sink=pos_sink, mass_sink=mass_sink, i_im=i_im, overwrite=overwrite, path_out=path_out,
                       col_dens_only=col_dens_only, tag=tag, map_size=map_size, vel_size=vel_size, ps=ps, all_sink=all_sink, ax1=c_axis, PSCAL=PSCAL, add_B=add_B, histogram=histogram, pressure=pressure, temperature=temperature, col_dens=col_dens, rho_and_vel=rho_and_vel, amr_level=amr_level)
    return center
##############################################################################

###########################################################################


##############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="core visualization with axis alignement")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()
    #readdir = args.read_dir
    #writedir = args.read_dir

    readdir, writedir = local_dir(args.read_dir)
    #outputs = search_ro(readdir)
    outputs = [1,6,12]
    nout = len(outputs)

    overwrite = True
    radius = 10.  # in AU
    nbins = 30
    c_center = np.array([0.5, 0.5, 0.5])
    c_axis = None
    # center_def = 'sink' #'baricenter'
    center_def = None
    rotate = None
    #rotate = 'AM'
    add_B = False
    #PSCAL = range(6)
    #PSCAL=[-1] + range(6)
    #PSCAL += [100] # PSCAL_only
    PSCAL = []
    ps = False  # save .eps file
    col_dens_only = False
    map_size = 100
    plot_sinks = False
    rho_and_vel = True
    col_dens = True
    histogram = False
    pressure = True
    temperature = True
    amr_level = True
    
    #zoom_v = 0.5 / 2**np.arange(3)
    #i_ims = [0, 1, 2]

    zoom_v = 0.5 / 2**np.arange(1)
    i_ims = [0]

    for iout in range(0, nout, 1):

        ioutput = outputs[iout]
        center_AU = make_image_zoom(readdir, ioutput, zoom_v, i_ims, path_out=writedir, plot_sinks=plot_sinks, overwrite=overwrite, center_img=c_center, make_phi=False, deli=',', col_dens_only=col_dens_only, tag='',
                                    gcomp=True, map_size=map_size, vel_size=20, xyz_sink=3, center_def=center_def, order='<', ps=ps, ind_sink=0, all_sink=True, rotate=rotate, make_slice=False, col_dens=col_dens, c_axis=c_axis, add_B=add_B, PSCAL=PSCAL, histogram=histogram, pressure=pressure, temperature=temperature, rho_and_vel=rho_and_vel, amr_level=amr_level)
