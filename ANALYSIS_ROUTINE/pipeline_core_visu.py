#! /usr/bin/env python

##### python dependences ##################################################################
import os
import argparse
import pymses
import numpy as np
#from pymses.utils import constants as C
#from pymses.utils.regions import Sphere
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
#from pymses import RamsesOutput
#from pymses.filters import CellsToPoints
#from pymses.sources.ramses.tree_utils import octree_compute_neighbors
#from pymses.analysis import Camera, raytracing, slicing
#from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
#from pymses.analysis.point_sampling import PointSamplingProcessor
#from pymses.utils.regions import Sphere
#from pymses.utils.regions import Box
#from pymses.filters import RegionFilter
#from pymses.filters import PointFunctionFilter
#from matplotlib.colors import LogNorm

##### module dependences ##################################################################
import  module_visualization as md_visu
from module_local_dir import local_dir

##############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="core visualization with axis alignement")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()
    #readdir = args.read_dir
    #writedir = args.read_dir

    readdir, writedir = local_dir(args.read_dir)
    outputs = search_ro(readdir)
    #outputs = range(1019,1021)
    #outputs = range(822,830,5)
    #outputs = [818,820]
    #outputs = [400,600,800,1000,1500,2000,2500,3000,3500]
    #outputs = range(3433,3419,-3)
    #outputs = [3399]
#    outputs = [1000,1005,1010]
    nout = len(outputs)

    overwrite = True
    radius = 100. # in AU
    nbins = 30
    zoom_v = 0.5 / 4**np.arange(6)
#    zoom_v = [0.0033,0.0016]#  / np.arange(1,4) / 2**3
#    i_ims = [1,2]
#    zoom_v = np.append(zoom_v, zoom_v[1]*0.5)
#hydro_lowturb
#    c_center = np.array([ 0.49737453,  0.5016709,   0.4951839 ])#([0.5,0.5,0.5])
#    c_axis =  np.array([ 0.21040191,  0.85262068, -0.47829804])
#pint
#    c_center = np.array([ 0.5002806,   0.50004352,  0.49998823])
#    c_axis = np.array([-0.43124987, -0.24715896, -0.86771884])
    c_center = np.array([0.5,0.5,0.5])
    c_axis = None
    #center_def = 'sink' #'baricenter'
    center_def = None
    rotate = None 
    #rotate = 'AM'
    add_B = True
    #PSCAL = range(6)
    #PSCAL=[-1] + range(6)
    #PSCAL += [100] # PSCAL_only
    PSCAL = []
    ps=True
    col_dens_only = False
    map_size = 400
    plot_sinks = False
#    zoom_v = zoom_v[1:]#  / np.arange(1,4) / 2**3
#    i_ims = i_ims[1:]
    #zoom_v = [0.0080,0.0033,0.0016]
    zoom_v = 0.5 ** 5 / 4**np.arange(3)
    i_ims = [3,4,5]
    zoom_v = 0.5 ** 3 / 4**np.arange(3)
    i_ims = [0,1,2]

    for iout in range(10,nout,10):#range(9,104,10):
        ioutput = outputs[iout]
        center_AU = md_visu.make_image_zoom(readdir,ioutput,zoom_v,i_ims,path_out=writedir,plot_sinks=plot_sinks,overwrite=overwrite,center_img=c_center,make_phi=False,deli=',',col_dens_only=col_dens_only,tag='',gcomp=True,map_size=map_size,vel_size=20,xyz_sink=3,center_def=center_def,order='<',ps=ps,ind_sink=0,all_sink=True,rotate=rotate,make_slice=False,col_dens=True,c_axis=c_axis,add_B=add_B,PSCAL=PSCAL)
#        print 'center', center_AU
        #md_visu.velocity_ana(readdir,ioutput,zoom_v,path_out=writedir,overwrite=overwrite,tag='',gcomp=True,center_def='baricenter',order='<',ps=False,radius=radius, nbins=nbins,rotate='AM',i_im=1)



