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
import module_extract as md_ex
#from module_core import find_baricenter, find_axis, renormalize_units, cal_M_R

##############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="core visualization with axis alignement")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()
    #readdir = args.read_dir
    #writedir = args.read_dir

    readdir, writedir = local_dir(args.read_dir)
    outputs = search_ro(readdir)
    nout = len(outputs)

    overwrite = True
    radius = 100. # in AU
    nbins = 30
    zoom_v = 0.5 / 4**np.arange(6)
    zoom_v = 0.5 / np.arange(1,4) / 2**3
#    zoom_v = np.append(zoom_v, zoom_v[1]*0.5)
#hydro_lowturb
#    c_center = np.array([ 0.49737453,  0.5016709,   0.4951839 ])#([0.5,0.5,0.5])
#    c_axis =  np.array([ 0.21040191,  0.85262068, -0.47829804])
#pint
#    c_center = np.array([ 0.5002806,   0.50004352,  0.49998823])
#    c_axis = np.array([-0.43124987, -0.24715896, -0.86771884])
    c_center = np.array([0.5,0.5,0.5])
    c_axis = None
    center_def = 'sink' #'baricenter'
    rotate = None #'AM'
    #rotate = 'AM'
    add_B = False
    PSCAL = range(4)
    #PSCAL=[]
    ps=False
    col_dens_only = False
    map_size = 400
    map_size = 2048
    md_ex.look2(readdir,227)
#        print 'center', center_AU
        #md_visu.velocity_ana(readdir,ioutput,zoom_v,path_out=writedir,overwrite=overwrite,tag='',gcomp=True,center_def='baricenter',order='<',ps=False,radius=radius, nbins=nbins,rotate='AM',i_im=1)



