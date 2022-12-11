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

#import module_analysis as md_ana
#import module_disc as md_disc
#import module_statistics as md_stat
from module_local_dir import local_dir

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Find ellipsoids for a series of radii")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()
    readdir, writedir = local_dir(args.read_dir)
    outputs = search_ro(readdir)
    nout = len(outputs)
    ro = RamsesOutput(readdir,outputs[0])#,order='<') 
    ncpu = ro.info["ncpu"]
    clump_all = []
    boxlen = ro.info['boxlen']**3
    mass_sol = ro.info['unit_mass'].express(C.Msun)/boxlen
    len_pc = ro.info['unit_length'].express(C.pc)
##index  lev   parent      ncell    peak_x             peak_y             peak_z             rho-               rho+               rho_av             mass_cl            relevance 
#    f1 = plt.figure(1)
    #f2 = plt.figure(2)
    for iout in [1,3,5,15,30,40]:
        clumpsave = readdir+'/clump_list'+str(outputs[iout]).zfill(5)+'.txt'
        if os.path.exists(clumpsave): clump_all = np.loadtxt(clumpsave)
        else:
          for icpu in range(1,ncpu+1):
            clumpfile = readdir+'/output_'+str(outputs[iout]).zfill(5)+'/clump_'+str(outputs[iout]).zfill(5)+'.txt'+str(icpu).zfill(5)
            if not os.path.exists(clumpfile): continue
            clumps = np.loadtxt(clumpfile,skiprows=1)
            print clumps.shape
            if len(clumps)>0:
                if len(clump_all)==0: clump_all = clumps
                else: clump_all = np.vstack((clump_all,clumps))
          np.savetxt(clumpsave,clump_all)
        print clump_all.shape
        print clump_all
        clump_all[:,10] = clump_all[:,10] * mass_sol
        print clump_all[:,4:7]
        print boxlen, len_pc
        #plt.scatter(clump_all[:,4],clump_all[:,5])
        #plt.show()
        #clump_all[:,4:7] = (clump_all[:,4:7] / boxlen) - 0.5
        clump_all[:,4:7] = clump_all[:,4:7]
        clump_cen = np.sum(clump_all[:,4:7]*clump_all[:,10,np.newaxis],axis=0) / np.sum(clump_all[:,10])
        print 'clump_cen', clump_cen
#        hist_mass, logmass_edges = np.histogram(np.log10(clump_all[:,10]),bins=16)
#        plt.figure(1)
#        plt.plot(np.repeat(logmass_edges,2), np.hstack((0,np.repeat(np.log10(hist_mass),2),0)))
        distance = np.linalg.norm(clump_all[:,4:7]-clump_cen[np.newaxis,:],axis=1)
        hist_mass_pos, logmass_edges, dist_edges = np.histogram2d(np.log10(clump_all[:,10]),distance, bins=(16,5))
        plt.figure()
        for idis in range(len(dist_edges)-1):
            plt.plot(np.repeat(logmass_edges,2), np.hstack((0,np.repeat(np.log10(hist_mass_pos[:,idis]),2),0)),label='%.3f pc'%(dist_edges[idis+1]))
            if idis==0:
                plt.plot(logmass_edges[-8:],logmass_edges[-8:]*(-1))
                plt.plot(logmass_edges[-8:],logmass_edges[-8:]*(-1.3))
        plt.plot(np.repeat(logmass_edges,2), np.hstack((0,np.repeat(np.log10(np.sum(hist_mass_pos,axis=1)),2),0)),'k')
        plt.legend()
#    plt.figure(1)
#    plt.plot(logmass_edges[-8:]+np.log10(mass_sol),logmass_edges[-8:]*(-1)+3)
#    plt.plot(logmass_edges[-8:]+np.log10(mass_sol),logmass_edges[-8:]*(-1.3)+3.5)
    plt.show()
