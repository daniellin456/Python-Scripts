import matplotlib
matplotlib.use('Agg') #_tkinter.TclError: no display name and no $DISPLAY environment variable #due to python2
import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import myfunction as mf
from numpy.core.fromnumeric import reshape
import numpy as np
import os
import pymses #running python2 rather than python3 due to pymses
from pymses.analysis import bin_spherical
from pymses.analysis import Camera, raytracing, slicing, ScalarOperator, FractionOperator, MaxLevelOperator
from pymses.filters import CellsToPoints
from pymses.filters import RegionFilter
from pymses.utils import constants as C
from pymses.utils.regions import Sphere
import scipy as sp
import scipy.special as s
import shutil
import sys #use sys.exit()

pwd = "/home/cmcheng/result/" #save results in this directory
inputfile = "alpha0.2" #ramses output directory name
start = 1 #for loop start (including this value)
end = 1654 #for loop end (excluding this value)
maplength = 30. #map length with au unit

##delete folder with old data
if os.path.isdir(pwd+inputfile+"/property/"):
	shutil.rmtree(pwd+inputfile+"/property/")

##make new empty folder
if not os.path.isdir(pwd+inputfile+"/property/"):
    os.makedirs(pwd+inputfile+"/property/")

##get resolution of op.info["levelmax"]=16
for i in range(1, 2): 
    op = pymses.RamsesOutput("/data/cmcheng/"+inputfile, i)
    unit_l_au=op.info["unit_length"].express(C.au)
    dr_bin = unit_l_au*1./2**16

oparray = range(start, end)
time = []
mapDarray = []
mapTarray = []
maxDarray = []
minDarray = []
maxTarray = []
minTarray = []
D_r = np.array([]) #average radial profile
T_r = np.array([]) #average radial profile
v_r = np.array([]) #average radial profile

for i in range(start, end): 
    op = pymses.RamsesOutput("/data/cmcheng/"+inputfile, i) #print(op.info) shows all parameters.
    
    ##unit	
    unit_M=op.info["unit_density"]*(op.info["unit_length"]**3)
    unit_M_kg=unit_M.express(C.kg)
    unit_M_g=unit_M.express(C.g)
    unit_M_sun=unit_M.express(C.Msun)
	
    unit_l_m=op.info["unit_length"].express(C.m)
    unit_l_cm=op.info["unit_length"].express(C.cm)
    unit_l_au=op.info["unit_length"].express(C.au)	

    unit_t_s=op.info["unit_time"].express(C.s)
    unit_t_y=op.info["unit_time"].express(C.year)
	
    unit_P_Pa=op.info["unit_pressure"].express(C.N/C.m**2)
	
    unit_D_SI=op.info["unit_density"].express(C.kg/C.m**3)
    unit_D_cgs=op.info["unit_density"].express(C.g/C.cm**3)
    unit_D_Hcc=op.info["unit_density"].express(C.H_cc)
	
    unit_acceler_SI=unit_l_m/unit_t_s**2
    unit_acceler_cgs=unit_l_cm/unit_t_s**2

    unit_vel_SI=unit_l_m/unit_t_s
    unit_vel_cgs=unit_l_cm/unit_t_s
    
    op_amr = op.amr_source(["rho","T", "P", "g", "vel"]) #It needs pymsesrc file.

    time.append(op.info["time"]*unit_t_y)
    
    ####map
    cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', region_size=[maplength/unit_l_au, maplength/unit_l_au], up_vector='y', map_max_size=512)
    density = ScalarOperator(lambda dset: dset["rho"], op.info["unit_density"])
    map_d = slicing.SliceMap(op_amr, cam, density, z=0.)
    map_density = map_d.map
    temperature = ScalarOperator(lambda dset: dset["T"], op.info["unit_temperature"])
    map_T = slicing.SliceMap(op_amr, cam, temperature, z=0.)
    map_temperature = map_T.map

    mapDarray.append(np.log10(map_density*unit_D_cgs))   
    mapTarray.append(map_temperature)
    
    plt.title("No."+str(i)+" Time: "+str(round(op.info["time"]*unit_t_y, 1))+" (year)")
    snapshot = plt.imshow(np.log10(map_density*unit_D_cgs), cmap='winter', extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength])
    plt.xlabel("Length (au)")
    plt.ylabel("Length (au)")
    plt.colorbar(snapshot, extend='max', label="Log Scale Density "+r"$\rm(g/cm^3)$")
    plt.savefig(pwd+inputfile+"/property/D Map"+str(i)+".png")
    plt.clf()
    
    plt.title("No."+str(i)+" Time: "+str(round(op.info["time"]*unit_t_y, 1))+" (year)")
    snapshot = plt.imshow(map_temperature, cmap='cool', extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength])
    plt.xlabel("Length (au)")
    plt.ylabel("Length (au)")
    plt.colorbar(snapshot, extend='max', label="Temperature (K)")
    plt.savefig(pwd+inputfile+"/property/T Map"+str(i)+".png")
    plt.clf()
    
    ####radial profile
    source = CellsToPoints(op_amr)
    cell = source.flatten()
    
	##physical quantities of cells
    cood = cell.points
    dx = cell.get_sizes()
    V = dx**3
    D = cell["rho"]
    M = V*D
    T = cell["T"]
    P = cell["P"]
    v = cell["vel"]
    cood_D_max = cood[np.argmax(cell["rho"])][:]
    x = cood[:,0]-cood_D_max[0]
    y = cood[:,1]-cood_D_max[1]
    z = cood[:,2]-cood_D_max[2]
    r = np.linalg.norm(cood-cood_D_max[np.newaxis,:],axis=1)
    vc = v[np.argmin(r),:]
    vr = np.sum((v-vc)*(cood-cood_D_max) ,axis=1)/r
    vr[r==0.]=0. #invalid value encountered in divide

    maxDarray.append(np.amax(D)*unit_D_cgs)
    minDarray.append(np.amin(D)*unit_D_cgs)
    maxTarray.append(np.amax(T))
    minTarray.append(np.amin(T))

    ##radius bin 
    r_1 = np.arange(0., 20.00001/unit_l_au, dr_bin/unit_l_au)
    r_2 = np.arange(21./unit_l_au, 50.00001/unit_l_au, 1./unit_l_au)
    r_3 = np.arange(100./unit_l_au, 10000.00001/unit_l_au, 100./unit_l_au)
    r_12 = np.append(r_1, r_2)
    r_bin = np.append(r_12, r_3)
    r_cent = 0.5*(r_bin[:-1] + r_bin[1:])

    ##mass & volume & density
    M_hist,_ = np.histogram(r, bins=r_bin, weights=M)
    V_hist,_ = np.histogram(r, bins=r_bin, weights=V)
    D_hist = M_hist/V_hist

    ##temperature
    TM_hist,_ = np.histogram(r, bins=r_bin, weights=T*M)
    T_hist = TM_hist/M_hist

    ##radial velocity
    vrM_hist,_ = np.histogram(r, bins=r_bin, weights=vr*M)	
    vr_hist = vrM_hist/M_hist

    D_r = np.append(D_r, D_hist*unit_D_cgs)
    T_r = np.append(T_r, T_hist)
    v_r = np.append(v_r, vr_hist*unit_vel_SI)

D_r = np.reshape(D_r, (len(time), -1))
T_r = np.reshape(T_r, (len(time), -1))
v_r = np.reshape(v_r, (len(time), -1))
np.savetxt(pwd+inputfile+"/property/alltime.txt", time)
np.savetxt(pwd+inputfile+"/property/r.txt", r_cent*unit_l_au)
np.savetxt(pwd+inputfile+"/property/D(r).txt", D_r)
np.savetxt(pwd+inputfile+"/property/T(r).txt", T_r)
np.savetxt(pwd+inputfile+"/property/vr(r).txt", v_r)
np.savetxt(pwd+inputfile+"/property/maxD.txt", [np.amax(maxDarray)])
np.savetxt(pwd+inputfile+"/property/minD.txt", [np.amin(minDarray)])
np.savetxt(pwd+inputfile+"/property/maxT.txt", [np.amax(maxTarray)])
np.savetxt(pwd+inputfile+"/property/minT.txt", [np.amin(minTarray)])

####make central temperature-density figure
nptime = np.array(time)
timeabove = nptime[T_r[:,0]>=11]
timebelow = nptime[T_r[:,0]<11]
fig = plt.figure(figsize=(12,8))
norm = colors.Normalize(vmin=np.amin(timeabove), vmax=np.amax(timeabove))
sm = cm.ScalarMappable(cmap=cm.rainbow, norm=norm)
sm.set_array([])
plt.colorbar(sm, aspect=30, pad=0.02, label="Time (year)")
for i in range(len(time)):
    if T_r[i,0] >= 11:
        plt.plot(D_r[i,0], T_r[i,0], "o", markersize=2, color=cm.rainbow(norm(timeabove[i-len(timebelow)])))   
    else:
        plt.plot(D_r[i,0], T_r[i,0], "o", markersize=2, color="black")
plt.xlabel("Central Density "+r"$\rm(g/cm^3)$")
plt.ylabel("Central Temperature (K)")
plt.xscale('log')
plt.yscale('log')
plt.xlim(10**-18, 1) 
plt.ylim(10**0.5, 10**5)
plt.tick_params(which='major',width=1, length=5) 
plt.tick_params(which='minor',width=1, length=3)
#plt.axvspan(xmin=10**-13, xmax=5.78e-09, facecolor='silver')
fig.savefig(pwd+inputfile+"/property/central.png", bbox_inches='tight')
plt.clf()

####calculate slope of temperature-density diagram
dlogD = np.diff(np.log10(D_r), axis=0)
dlogT = np.diff(np.log10(T_r), axis=0)
np.savetxt(pwd+inputfile+"/property/slope.txt", dlogT/dlogD)

####make profile figure
for i in range(len(time)):
    plt.figure(figsize=(20,6))
    plt.suptitle("No."+str(oparray[i])+" Time: "+str(round(time[i], 1))+" (year)")
    plt.subplot(131)
    plt.plot(r_cent*unit_l_au, D_r[i], "-", markersize=2, color="black")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(dr_bin/2, 10000) 
    plt.ylim(10**-19, 1) 
    plt.xlabel("Radius (au)")
    plt.ylabel("Density "+r"$\rm(g/cm^3)$")
    plt.tick_params(which='major',width=1, length=5) 
    plt.tick_params(which='minor',width=1, length=3)
    plt.subplot(132)
    plt.plot(r_cent*unit_l_au, T_r[i], "-", markersize=2, color="black")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(dr_bin/2, 10000) 
    plt.ylim(5, 10**5) 
    plt.xlabel("Radius (au)")
    plt.ylabel("Temperature "+r"$\rm(K)$")
    plt.tick_params(which='major',width=1, length=5) 
    plt.tick_params(which='minor',width=1, length=3)
    plt.subplot(133)
    plt.plot(r_cent*unit_l_au, v_r[i], "-", markersize=2, color="black")
    plt.xscale('log')
    plt.xlim(dr_bin/2, 10000)
    plt.ylim(np.nanmin(v_r), np.nanmax(v_r))
    plt.axhline(y=0., c="r", ls="--", lw=1)
    plt.xlabel("Radius (au)")
    plt.ylabel("Radial Velocity "+r"$\rm(m/s)$")
    plt.tick_params(which='major',width=1, length=5) 
    plt.tick_params(which='minor',width=1, length=3)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(pwd+inputfile+"/property/radius"+str(oparray[i])+".png", bbox_inches='tight')
    plt.clf()

####make map animation
fig2 = plt.figure()
def animate2(i):
    fig2.clear() #clear colorbar
    plt.title("No."+str(oparray[i])+" Time: "+str(round(time[i], 1))+" (year)")
    s = plt.imshow(mapDarray[i], cmap = 'winter', vmin=np.amin(mapDarray), vmax=np.amax(mapDarray), extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength]) #fixed colorbar coordinate 
    #s = plt.imshow(densityarray[i], cmap = 'winter', extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength]) #floated colorbar coordinate
    plt.xlabel("Length (au)")
    plt.ylabel("Length (au)")
    fig2.colorbar(s, extend='max', label="Log Scale Density "+r"$\rm(g/cm^3)$")

if len(oparray)-1<600:
    frameInterval = 600//(len(oparray)-1)
else:
    frameInterval = 1
anim = animation.FuncAnimation(fig2, animate2, interval=frameInterval, frames=range(end-start)) 
anim.save(pwd+inputfile+"/property/D Map.gif", writer="imagemagick")

"""
fig4 = plt.figure()
def animate4(i):
    fig4.clear() #clear colorbar
    plt.title("No."+str(oparray[i])+" Time: "+str(round(time[i], 1))+" (year)")
    s = plt.imshow(mapTarray[i], cmap = 'cool', vmin=np.amin(mapTarray), vmax=np.amax(mapTarray), extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength]) #fixed colorbar coordinate 
    #s = plt.imshow(temperaturearray[i], cmap = 'cool', extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength]) #floated colorbar coordinate
    plt.xlabel("Length (au)")
    plt.ylabel("Length (au)")
    fig4.colorbar(s, extend='max', label="Temperature (K)")

if len(oparray)-1<600:
    frameInterval = 600//(len(oparray)-1)
else:
    frameInterval = 1
anim = animation.FuncAnimation(fig4, animate4, interval=frameInterval, frames=range(end-start)) 
anim.save(pwd+inputfile+"/property/T Map.gif", writer="imagemagick")

####make map-profile animation
fig = plt.figure(figsize=(12,12))
def animate(i):
    fig.clear() #clear colorbar
    plt.suptitle("No."+str(oparray[i])+" Time: "+str(round(time[i], 1))+" (year)")
    plt.subplot(221)
    dmap = plt.imshow(mapDarray[i], cmap = 'winter', vmin=np.amin(mapDarray), vmax=np.amax(mapDarray), extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength])
    plt.xlabel("Length (au)")
    plt.ylabel("Length (au)")
    plt.colorbar(dmap, extend='max', label="Log Scale Density "+r"$\rm(g/cm^3)$")
    plt.subplot(222)
    plt.plot(r_cent*unit_l_au, D_r[i], "-", markersize=2, color="black")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(dr_bin/2, 10000) 
    plt.ylim(10**-19, 1) 
    plt.xlabel("Radius (au)")
    plt.ylabel("Density "+r"$\rm(g/cm^3)$")
    plt.tick_params(which='major',width=1, length=5) 
    plt.tick_params(which='minor',width=1, length=3)
    plt.subplot(223)
    tmap = plt.imshow(mapTarray[i], cmap = 'cool', vmin=np.amin(mapTarray), vmax=np.amax(mapTarray), extent=[-0.5*maplength, 0.5*maplength, -0.5*maplength, 0.5*maplength])
    plt.xlabel("Length (au)")
    plt.ylabel("Length (au)")
    plt.colorbar(tmap, extend='max', label="Temperature (K)")
    plt.subplot(224)
    plt.plot(r_cent*unit_l_au, v_r[i], "-", markersize=2, color="black")
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlim(dr_bin/2, 10000) 
    plt.ylim(np.nanmin(v_r), np.nanmax(v_r))
    plt.xlabel("Radius (au)")
    plt.ylabel("Radial Velocity"+r"$\rm(m/s)$")
    plt.tick_params(which='major',width=1, length=5) 
    plt.tick_params(which='minor',width=1, length=3)
    fig.tight_layout(rect=[0, 0, 1, 0.95])	

if len(oparray)-1<600:
    frameInterval = 600//(len(oparray)-1)
else:
    frameInterval = 1
anim = animation.FuncAnimation(fig, animate, interval=frameInterval, frames=range(end-start)) 
anim.save(pwd+inputfile+"/property/radius.gif", writer="imagemagick")
plt.clf()
"""

####M(rho, T) histogram
##delete folder with old data
if os.path.isdir(pwd+inputfile+"/histogram2D/"):
	shutil.rmtree(pwd+inputfile+"/histogram2D/")

####2D histogram
##make new empty folder
if not os.path.isdir(pwd+inputfile+"/histogram2D/"):
    os.makedirs(pwd+inputfile+"/histogram2D/")

maxD = np.loadtxt(pwd+inputfile+"/property/maxD.txt")
minD = np.loadtxt(pwd+inputfile+"/property/minD.txt")
maxT = np.loadtxt(pwd+inputfile+"/property/maxT.txt")
minT = np.loadtxt(pwd+inputfile+"/property/minT.txt")
time2 = []
Mzarray = []

for i in range(start, end): 
    op = pymses.RamsesOutput("/data/cmcheng/"+inputfile, i)
    unit_M=op.info["unit_density"]*(op.info["unit_length"]**3)
    unit_M_sun=unit_M.express(C.Msun)
    unit_D_cgs=op.info["unit_density"].express(C.g/C.cm**3)
    unit_t_y=op.info["unit_time"].express(C.year)
    op_amr = op.amr_source(["rho","T", "P", "g", "vel"])
    time2.append(op.info["time"]*unit_t_y)
    source = CellsToPoints(op_amr)
    cell = source.flatten()
    dx = cell.get_sizes()
    V = dx**3
    D = cell["rho"]
    M = V*D
    T = cell["T"]

    ##M(rho, T)
    D_bin_log = np.linspace(np.log10(minD/unit_D_cgs), np.log10(maxD/unit_D_cgs), 50)
    #D_bin_log = np.linspace(np.log10(10**-15/unit_D_cgs), np.log10(10**-7/unit_D_cgs), 50)
    D_bin = 10**D_bin_log

    T_bin_log = np.linspace(np.log10(minT), np.log10(maxT), 50)
    #T_bin_log = np.linspace(np.log10(1), np.log10(10000), 50)
    T_bin = 10**T_bin_log

    M_2Dhist, D_, T_ = np.histogram2d(D, T, bins=(D_bin, T_bin), weights=M)
    M2D = np.where(M_2Dhist*unit_M_sun>0, M_2Dhist*unit_M_sun, np.nan) # To mask zero for log scale colormap
    Mzarray.append(M2D)
    np.savetxt(pwd+inputfile+"/histogram2D/D_axis.txt", D_*unit_D_cgs)
    np.savetxt(pwd+inputfile+"/histogram2D/T_axis.txt", T_)

    fig = plt.figure()
    plt.title("No."+str(i)+" Time: "+str(round(op.info["time"]*unit_t_y, 1))+" (year)")
    snapshot = plt.pcolormesh(D_*unit_D_cgs, T_, M2D, norm = colors.LogNorm(vmin=np.nanmin(M2D), vmax=np.nanmax(M2D)), cmap='rainbow')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Density "+r"$\rm(g/cm^3)$")
    plt.ylabel("Temperature "+r"$\rm(K)$")
    plt.tick_params(which='major',width=1, length=5) 
    plt.tick_params(which='minor',width=1, length=3)
    plt.colorbar(snapshot, extend='max', label="Mass "+r"$\rm(M_{\odot})$")
    plt.tight_layout()
    plt.savefig(pwd+inputfile+"/histogram2D/M(rho,T)"+str(i)+".png")
    plt.clf()

np.save(pwd+inputfile+"/histogram2D/M.npy", Mzarray)

####M(rho, T) animation
Dx = np.loadtxt(pwd+inputfile+"/histogram2D/D_axis.txt")
Ty = np.loadtxt(pwd+inputfile+"/histogram2D/T_axis.txt")

fig = plt.figure()
def animate(i):
    fig.clear() #clear colorbar
    plt.title("No."+str(oparray[i])+" Time: "+str(round(time2[i], 1))+" (year)") 
    Mz = plt.pcolormesh(Dx, Ty, Mzarray[i], norm=colors.LogNorm(vmin=np.nanmin(Mzarray), vmax=np.nanmax(Mzarray)), cmap='rainbow')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Density "+r"$\rm(g/cm^3)$")
    plt.ylabel("Temperature "+r"$\rm(K)$")
    plt.tick_params(which='major',width=1, length=5) 
    plt.tick_params(which='minor',width=1, length=3)
    plt.colorbar(Mz, extend='max',label="Mass "+r"$\rm(M_{\odot})$")
    plt.tight_layout()

if len(oparray)-1<600:
    frameInterval = 600//(len(oparray)-1)
else:
    frameInterval = 1
anim = animation.FuncAnimation(fig, animate, interval=frameInterval, frames=range(end-start)) 
anim.save(pwd+inputfile+"/histogram2D/M(rho,T).gif", writer="imagemagick")
plt.clf()

"""
####resolution analysis to decide angle bin
pwd = "/home/cmcheng/result/" #save results in this directory
inputfile = "alpha0.4" #ramses output directory name
start = 208 #for loop start (including this value) alpha0.4=208 alpha0.6=269 alpha0.8=467
end = 318 #for loop end (excluding this value) alpha0.4=318 alpha0.6=452 alpha0.8=682
dr_bin = 1.0 #unit: au

firstCoreTime = []
NumPerShell = []
for i in range(start, end): 
    op = pymses.RamsesOutput("/data/cmcheng/"+inputfile, i) #print(op.info) shows all parameters. 
    unit_l_au=op.info["unit_length"].express(C.au)
    unit_t_y=op.info["unit_time"].express(C.year)
    op_amr = op.amr_source(["rho","T", "P", "g", "vel"]) #It needs pymsesrc file.		
    firstCoreTime.append(op.info["time"]*unit_t_y)
    
    ##region filter
    center = [0.5, 0.5, 0.5]
    radius = 15./unit_l_au
    region = Sphere(center, radius)
    op_fil = RegionFilter(region, op_amr)
    source = CellsToPoints(op_fil)
    cell = source.flatten()

	##physical quantities of cells
    cood = cell.points
    cood_D_max = cood[np.argmax(cell["rho"])][:]
    r = np.linalg.norm(cood-cood_D_max[np.newaxis,:],axis=1)
    v = cell["vel"]
    vc = v[np.argmin(r),:]
    vr = np.sum((v-vc)*(cood-cood_D_max) ,axis=1)/r
    vr[r==0.]=0. #invalid value encountered in divide

    ##radius bin
    r_bin = np.arange(0., 9.00001/unit_l_au, dr_bin/unit_l_au)
    r_cent = 0.5*(r_bin[:-1] + r_bin[1:])
    #making condition list [0<=r<r_bin1, r_bin1<=r<r_bin2,..., r_binn-1<=r<=r_binn] 
    condit =[]
    for n in np.arange(len(r_cent)):
        if n != np.max(np.arange(len(r_cent))):
            condit.append((dr_bin*n<=r*unit_l_au)&(r*unit_l_au<dr_bin*(n+1)))
        else:
            condit.append((dr_bin*n<=r*unit_l_au)&(r*unit_l_au<=dr_bin*(n+1)))
        #print(dr_bin*n, dr_bin*(n+1)) #check condition
	#number of cells in each cell
    for i in condit:
        NumPerShell.append(len(vr[i]))

NumPerShell = np.reshape(np.array(NumPerShell),(len(firstCoreTime), -1))
np.savetxt(pwd+inputfile+"/property/NumberOfCell.txt", NumPerShell, fmt="%i")
resol16 = unit_l_au*1./2**16
thetaResol = 2*np.arctan(0.5*resol16/(r_bin[1:]*unit_l_au))
np.savetxt(pwd+inputfile+"/property/degResolution.txt", np.rad2deg(thetaResol))
"""
