#! /usr/bin/env python

import argparse
import os
import pymses
import pylab as plt
import numpy as np
from pymses.utils import constants as C
from pymses.utils.regions import Sphere, Box
from pymses.filters import RegionFilter
#from pymses.analysis import sample_points
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
#from pymses.analysis.visualization import Camera, SliceMap, ScalarOperator, raytracing, FractionOperator, ScalarOperator
#from pymses.analysis.visualization import fft_projection
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib.colors import LogNorm
import pickle as pc
from scipy.spatial.distance import cdist
import bisect
from scipy import stats
##### This module uses the output series of one simulation to properties time series for all sinks.
##### All quantities are stored in code units. File: sink_property_list.save
def newlist(oldlist, idsnew):
    return [oldlist[idnew] for idnew in idsnew]

def read_sink_csv(filename,deli=None):
    if deli is None: sinks = np.loadtxt(filename)
    else: sinks = np.loadtxt(filename,delimiter=deli)
    if len(sinks) == 0: return 0
    if len(sinks.shape)==1: sinks = sinks[np.newaxis,:]
    return sinks
# Id     M[Msol]    dmf[Msol]      x           y           z        vx       vy       vz     rot_period[y] lx/|l|  ly/|l|  lz/|l| acc_rate[Msol/y] acc_lum[Lsol]  age[y]  int_lum[Lsol]   Teff[K] 
# generate history for each sink
def t_list(outputs, ids, masses, positions, velocities, momentums, ages, birth, ioutputs,deli=None):
    noutputs = len(outputs)
    for i in range(noutputs)[:]:
        iout = outputs[i]
        ro = pymses.RamsesOutput(readdir, iout)
        t = ro.info["time"]*ro.info["unit_time"].express(C.year)
        sinks = read_sink_csv(readdir+'/output_'+str(iout).zfill(5)+'/sink_'+str(iout).zfill(5)+'.csv',deli)
        if sinks is 0: continue
#        sinks[:,1] = sinks[:,1]/boxlen**3
        sinks[:,2:5] = sinks[:,2:5]/boxlen
        nsinks = sinks.shape[0]
        for i in range(nsinks):
            if not sinks[i,0] in ids:#append new sinks
                ids.append(sinks[i,0])
                masses.append([sinks[i,1]])
                positions.append([sinks[i,3:6]])
                velocities.append([sinks[i,6:9]])
                momentums.append([sinks[i,10:13]])
                ages.append([sinks[i,15]])
                birth.append(t)
                ioutputs.append([iout])
            else:#follow sink history
                indi = ids.index(sinks[i,0])
                masses[indi].append(sinks[i,1])
                positions[indi].append(sinks[i,3:6])
                velocities[indi].append(sinks[i,6:9])
                momentums[indi].append(sinks[i,10:13])
                ages[indi].append(sinks[i,15])
                ioutputs[indi].append(iout)
        #remove merged sinks
        idsnew = [int(ids.index(id)) for id in sinks[:,0]]
        ids = newlist(ids,idsnew)
        masses = newlist(masses,idsnew)
        positions = newlist(positions,idsnew)
        velocities = newlist(velocities,idsnew)
        momentums = newlist(momentums,idsnew)
        ages = newlist(ages,idsnew)
        birth = newlist(birth,idsnew)
        ioutputs = newlist(ioutputs,idsnew)
    return ids, masses, positions, velocities, momentums, ages, birth, ioutputs

def save_sink_list(readdir,writedir,deli=None,overwrite=False):
    filename=writedir+'/sink_property_list.save'
    if os.path.exists(filename) and not overwrite: return
    outputs = search_ro(readdir)
    ids, masses, positions, velocities, momentums, ages, birth, ioutputs = [], [], [], [], [], [], [], []
    ids, masses, positions, velocities, momentums, ages, birth, ioutputs = t_list(outputs, ids, masses, positions, velocities, momentums, ages, birth, ioutputs, deli=deli)
    ids = [int(id) for id in ids]
    filename=writedir+'/sink_property_list.save'
    f = open(filename,'w')
    for list in [ids, masses, positions, velocities, momentums, ages, birth, ioutputs]:
        pc.dump(list,f)
    f.close()

#### Read sink property file for one simulation (all in code units) 
#### File: sink_property_list.save
def read_pickle(filename):
    f = open(filename,'rb')
    vars = []
    for i in range(8):
        vars.append(pc.load(f))
    f.close()
    return vars

#### This module uses the output series of one simulation to mass evolution for all sinks.
#### Calculate time for accreting a fraction of mass (0.6, 0.9)
#### File: sink_M_t.eps
aj = np.pi**(2.5)/6.
Cs = 2.e4 #cm/s
G = 6.67e-8 #cgs
pc_cm = 3.e18 #cm
Ms = 2.e33 #g
eta = 0.4
Myr = 86400.*365.*1.e6 #s
yr = 86400.*365.
def calc_M_of_R(R,Mach2):
    M = aj**(2./3.) * Cs**2*(R*pc_cm)/G*(1+ Mach2/3.*R**(2*eta)) #g
    rho = M*3./4./np.pi/(R*pc_cm)**3 #g/cm3
    t_ff = np.sqrt(3*np.pi/32./G/rho)
    return M, rho, t_ff

def t_accret(ages, masses,frac_list):
    nfrac = len(frac_list)
    t_acc = [[] for i in range(nfrac)]
    for age, mass in zip(ages,masses):
        for i in range(nfrac):
            t_acc[i].append(age[bisect.bisect(mass,mass[-1]*frac_list[i])]*time_yr)
    return t_acc

def plot_M_t(readdir,writedir,new_tacc=True):
    fig,ax2 = plt.subplots(1, figsize=(6,3.2))
    filename=writedir+'/sink_property_list.save'
    [ids, masses, positions, velocities, momentums, ages, birth, ioutputs] = read_pickle(filename)
    for i in range(len(ids)):
        ax2.plot(np.asarray(ages[i])+birth[i], masses[i])
    plt.show()

def sink_M_t(readdir,writedir,new_tacc=True):
   # fig, [ax1,ax2] = plt.subplots(2,sharex=True, figsize=(6,7))
    fig,ax2 = plt.subplots(1, figsize=(6,3.2))
    filename=writedir+'/sink_property_list.save'
    col=['b','g','r','c','m','y']
    hat=['/','|','-','o','O']
    filename2=writedir+'/sink_tacc.save'
    frac_list=[0.9,0.6]
    if new_tacc:
        [ids, masses, positions, velocities, momentums, ages, birth, ioutputs] = read_pickle(filename)
        t_acc = t_accret(ages,masses,frac_list)
        t_acc9 = t_acc[0]
        t_acc6 = t_acc[1]
        M_end = [m[-1]* mass_sun for m in masses]
        f = open(filename2,'w')
        for list in [t_acc9,t_acc6,M_end]:
            pc.dump(list,f)
        f.close()
    else: 
        [t_acc9,t_acc6,M_end] = read_pickle(filename2)
        t_acc = [t_acc9,t_acc6]
    ifrac=1
    ntbins = 20
    ## histogram of age at fraction (0.9, 0.6) of mass accreted
    bins = np.logspace(np.log10(min(t_acc[ifrac])),np.log10(max(t_acc[ifrac])),ntbins)
    hist, bins = np.histogram(t_acc9,bins=bins)
#    ax1.semilogx(np.repeat(bins,2),np.hstack((0,np.repeat(hist,2),0)),'k:',lw=1,label=r'$t_\mathrm{acc,90}$')
    hist, bins = np.histogram(t_acc6,bins=bins)
#    ax1.semilogx(np.repeat(bins,2),np.hstack((0,np.repeat(hist,2),0)),'k',lw=1,label=r'$t_\mathrm{acc,60}$')

    ## histogram of accretion age (0.6) with mass bins by dex
    ifrac=0
    nbins = 5
    abin_edges = np.logspace(np.log10(min(t_acc[ifrac])),np.log10(max(t_acc[ifrac])),ntbins)
    mbin_edges = np.array([0,0.01,0.1,1,10,100])#np.logspace(-2.-nbins/2,-1+nbins/2,nbins+1)
    hist, abins, mbins = np.histogram2d(t_acc[ifrac],M_end,bins=[abin_edges,mbin_edges])
    bars = []
#    for im in range(len(mbin_edges)-2):
#        ax1.semilogx(np.repeat(abins,2),np.hstack((0,np.repeat(hist[:,im],2),0)),c=col[im],lw=2,label=str(mbin_edges[im+1])+' $M_\odot$')
##    plt.xlabel(r'$t_\mathrm{acc} (\mathrm{yr})$',fontsize=20)
#    ax1.set_xscale('log')
#    ax1.set_ylabel('Number of sinks',fontsize=20)
#    ax1.tick_params(axis='both',which='major',labelsize=15)
#    ax1.legend(loc='upper right',frameon=False,fontsize=14,labelspacing=0.1)
##    ax2 = ax1.twinx()
    ## Fractional  mass (0.6) versus accretion age
    ax2.loglog(t_acc[ifrac],np.asarray(M_end),'o',c='r',alpha=0.25,markeredgecolor='none')
#    t_acc_avg = np.mean(t_acc9)
#    t_acc_sig = np.std(t_acc9)
#    xpos = np.max(t_acc9)*0.7
#    ypos = np.min(M_end)
#    plt.text(xpos, ypos*10, '90% mass at '+str(int(t_acc_avg))+'$\pm$'+str(int(t_acc_sig))+' yr')
#    t_acc_avg = np.mean(t_acc6)
#    t_acc_sig = np.std(t_acc6)
#    plt.text(xpos, ypos*5, '60% mass at '+str(int(t_acc_avg))+'$\pm$'+str(int(t_acc_sig))+' yr')
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(t_acc[ifrac]),np.log(M_end))
    for run, R in zip(['isothermal','lesscompact','compact','extreme'],[0.75,0.33,0.084,0.042]):
        if readdir.find(run)>=0: Rc = R*1.3; break
    Mc = 1000 #Ms
    R = np.logspace(-5,-1,50) #pc
    Mach2 = 3.*G* (Mc*Ms)/3./(Rc*pc_cm)/ Cs**2 /2.
    M, rho, t_ff = calc_M_of_R(R,Mach2)
    ax2.loglog(t_ff*0.6/yr, M/Ms, label = r'$R=%.2f\mathrm{ pc}$'%(Rc))
    ax2.loglog(t_acc[ifrac],np.exp(intercept)*t_acc[ifrac]**slope,'k',zorder=0)
    ax2.text(2.e2,1.e1,r'$M_\mathrm{f} \propto t_\mathrm{acc,60}^{%.1f}$'%(slope),fontsize=16)
    ax2.set_xlim([1.e2,1.e6])
    ax2.set_ylim([1.e-4,1.e2])
    #plt.text(xpos, ypos*2.5, [round(slope,2), round(intercept,2), round(r_value,2), round(p_value,2)])
    ax2.set_xlabel(r'$t_\mathrm{acc} (\mathrm{yr})$',fontsize=20)
    ax2.set_ylabel(r'$M_\mathrm{f} (M_\odot)$',fontsize=20)
    plt.setp(ax2.get_yticklabels()[-2:], visible=False) 
    ax2.tick_params(axis='both',which='major',labelsize=15)
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.savefig(writedir+'/sink_M_t_2.pdf')
    plt.show()

#### This module calculates the distance nearest neighbor of each sink at the moment of birth 
#### File: nearest_neighbor.eps

def nearest_neighbor(readdir,writedir):
    filename=writedir+'/sink_property_list.save'
    [ids, masses, positions, velocities, momentums, ages, birth, ioutputs] = read_pickle(filename)
    nsinks = len(ids)
    near_neig = []
    dis_err = []
    for isink in range(1,nsinks):
        id = ids[isink]
        out = ioutputs[isink][0]
        old_sink_list = [i for i, e in enumerate(ioutputs) if e[0] <= out] #list of existing sinks
        old_sink_list.remove(isink)
        cur_pos = [positions[i][ioutputs[i].index(out)] for i in old_sink_list]
        cur_vel = [velocities[i][ioutputs[i].index(out)] for i in old_sink_list]
        diss = cdist(positions[isink][0][np.newaxis,:],np.asarray(cur_pos),'euclidean')[0]
        diss0 = cdist((positions[isink][0]-ages[isink][0]*velocities[isink][0])[np.newaxis,:],np.asarray(cur_pos)-ages[isink][0]*np.asarray(cur_vel),'euclidean')[0]
        sig = abs(diss-diss0) 
        print 'id',id,'out',out,'min dis', min(diss0)*len_au, sig[np.argmin(diss0)]*len_au
        near_neig.append(min(diss0)*len_au)
        dis_err.append(sig[np.argmin(diss0)]*len_au)
    dismin = min(near_neig)
    dismax = max(near_neig)
    nbin = 15
    bin_edges = np.logspace(np.log10(dismin),np.log10(dismax),nbin+1)
    hist, _ = np.histogram(near_neig,bins = bin_edges)
    err, _ = np.histogram(near_neig,bins = bin_edges,weights=np.asarray(dis_err)**2)
    hist_err = np.sqrt(err/hist)
    hist_cum = np.cumsum(hist)
    err_cum = np.cumsum(err)
    hist_cum_err = np.sqrt(err_cum/hist_cum)
    bin_c = np.sqrt(bin_edges[:-1]*bin_edges[1:])
#    plt.step(bin_c,hist,where='mid',c='b')
    plt.errorbar(bin_c,hist,xerr=hist_err,fmt='.',c='b')
#    plt.step(bin_c,hist_cum,where='mid',c='g')
    plt.errorbar(bin_c,hist_cum,xerr=hist_cum_err,fmt='.',c='g')
    plt.xscale('log')
    plt.xlim([dismin,dismax])
    plt.title(args.out_dir)
    plt.xlabel('distance to nearest sink (AU)')
    plt.savefig(writedir+'/nearest_neighbor.eps')
    plt.show()

#### This module calculates the distance of nearby newly formed stars 
#### File: nearby_formation.eps
def nearby_formation(readdir,writedir):
    filename=writedir+'/sink_property_list.save'
    [ids, masses, positions, velocities, momentums, ages, birth, ioutputs] = read_pickle(filename)
    nsinks = len(ids)
    neig_dis = []
    neig_dis0 = []
    age_cen = []
    dis_err = []
    mass_cen = []
    for isink in range(nsinks-1):
        id = ids[isink]
        out = ioutputs[isink][0]
        young_sink_list = [i for i, e in enumerate(birth) if e >= birth[isink]] #list of newly-formed sinks
        young_sink_list.remove(isink)
        #print 'young sink list', len(young_sink_list), young_sink_list
        #print 'positions[i][0]', len([positions[i][0] for i in young_sink_list]), [positions[i][0] for i in young_sink_list] 
        form_pos = np.asarray([positions[i][0] for i in young_sink_list])
        form_vel = np.asarray([velocities[i][0] for i in young_sink_list])
        inst_pos = np.asarray([positions[isink][ioutputs[isink].index(ioutputs[i][0])] for i in young_sink_list])
        inst_vel = np.asarray([velocities[isink][ioutputs[isink].index(ioutputs[i][0])] for i in young_sink_list])
        

        form_age = np.asarray([ages[i][0] for i in young_sink_list])
        form_pos0 = form_pos - form_age[:,np.newaxis]*form_vel
        inst_pos0 = inst_pos - form_age[:,np.newaxis]*inst_vel
        #print 'outputs dist', [ioutputs[i][0]-ioutputs[isink][0] for i in young_sink_list]
        inst_age = np.asarray([ages[isink][ioutputs[isink].index(ioutputs[i][0])]-ages[i][0] for i in young_sink_list])
        inst_mass = np.asarray([masses[isink][ioutputs[isink].index(ioutputs[i][0])]-ages[i][0] for i in young_sink_list]) 
        #print 'form_pos', form_pos
        diss = cdist( form_pos, inst_pos ,'euclidean' )[0] #distance calculated directly at output time
        diss0 = cdist( form_pos0, inst_pos0 ,'euclidean' )[0] #distance estimated at formation time
        sig = abs(diss-diss0)
        
        #print 'id',id,'out',out,'min dis', min(diss0)*len_au, sig[np.argmin(diss0)]*len_au
        age_cen.append(inst_age*time_yr)
        mass_cen.append(inst_mass*mass_sun)
        neig_dis.append(diss*len_au)
        neig_dis0.append(diss0*len_au)
        dis_err.append(sig*len_au)

    ineig = 9
    #neig_dis0 = neig_dis

    for isink in range(nsinks-1-ineig):
        #ord = np.argsort(age_cen[isink])
        #plt.plot(neig_dis[isink][ord], age_cen[isink][ord])
        ord = np.argsort(neig_dis0[isink])
        #print mass_cen[isink][ord[:10]],neig_dis0[isink][ord[:10]]
        plt.plot(mass_cen[isink][ord[0:11:9]],neig_dis0[isink][ord[0:11:9]],'-',zorder=isink)#,marker='o',linestyle='-')

    plt.gca().set_color_cycle(None)
    for isink in range(nsinks-1-ineig):
        ord = np.argsort(neig_dis0[isink])
        plt.scatter(mass_cen[isink][ord[9]],neig_dis0[isink][ord[9]],s=40,edgecolor='none',zorder=isink)
    plt.gca().set_color_cycle(None)
    for isink in range(nsinks-1-ineig):
        ord = np.argsort(neig_dis0[isink])
        plt.scatter(mass_cen[isink][ord[0]],neig_dis0[isink][ord[0]],s=20,edgecolor='none',zorder=isink)
        #plt.scatter(masses[isink][-1]*mass_sun,np.mean(np.sort(neig_dis0[isink])[:3]))
        #print neig_dis0[isink].shape #np.sort(neig_dis0[isink])
        #plt.scatter(masses[isink][-1]*mass_sun,np.sort(neig_dis0[isink])[ineig])
    plt.xlabel('Sink mass (M$_\odot$)',fontsize=18)
    plt.ylabel('Neighbor distance (AU)',fontsize=18)
    #plt.title('1st and 10th neighbors')
#    dismin = min(near_neig)
#    dismax = max(near_neig)
#    nbin = 15
#    bin_edges = np.logspace(np.log10(dismin),np.log10(dismax),nbin+1)
#    hist, _ = np.histogram(near_neig,bins = bin_edges)
#    err, _ = np.histogram(near_neig,bins = bin_edges,weights=np.asarray(dis_err)**2)
#    hist_err = np.sqrt(err/hist)
#    hist_cum = np.cumsum(hist)
#    err_cum = np.cumsum(err)
#    hist_cum_err = np.sqrt(err_cum/hist_cum)
#    bin_c = np.sqrt(bin_edges[:-1]*bin_edges[1:])
#    plt.errorbar(bin_c,hist,xerr=hist_err,fmt='.',c='b')
#    plt.errorbar(bin_c,hist_cum,xerr=hist_cum_err,fmt='.',c='g')
    plt.yscale('log')
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlim([1.e-3,1.e2])
    plt.ylim([1.e1,1.e5])
#    plt.title(args.out_dir)
#    plt.xlabel('distance to newly-formed sink (AU)')
#    plt.ylabel('age (yr)')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tight_layout()
    plt.savefig(writedir+'/neighbor_distance.pdf')
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate profiles")
    parser.add_argument("out_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()

    readdir = '/gpfs/data1/phennebe/'+args.out_dir
    writedir = '/gpfs/data1/ylee/'+args.out_dir
    if not os.path.exists(writedir): os.makedirs(writedir)

    ro = pymses.RamsesOutput(readdir,search_ro(readdir)[0])
    len_pc = ro.info["unit_length"].express(C.pc)
    boxlen = ro.info["boxlen"]
    dens = ro.info["unit_density"].express(C.H_cc)
    dens_par = (ro.info["unit_density"]*ro.info["unit_length"]).express(C.H_cc*C.cm)
    vel = ro.info["unit_velocity"].express(C.m/C.s)
    e_dens = (ro.info["unit_density"] * ro.info["unit_velocity"]**2).express(C.J/C.m**3)
    mass_sun = (ro.info["unit_mass"].express(C.Msun))
    temp_K = ro.info["unit_temperature"].express(C.K)
    lmax = ro.info["levelmax"]
    mag_gauss = ro.info["unit_mag"].express(C.Gauss)*0.5 * 1.e6
    time_yr = ro.info["unit_time"].express(C.year)
    len_au = ro.info["unit_length"].express(C.au)

    save_sink_list(readdir,writedir,deli=',',overwrite=True)
    plot_M_t(readdir,writedir)
#    sink_M_t(readdir,writedir,new_tacc=True)
    #plt.clf()
    #nearest_neighbor(readdir,writedir)
#    nearby_formation(readdir,writedir)
