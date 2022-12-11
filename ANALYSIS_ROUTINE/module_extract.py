#module to extract structure and visualise simulations
import sys
import numpy as np
import os 

import pymses
from pymses.filters import CellsToPoints
from pymses.sources.ramses.tree_utils import octree_compute_neighbors
from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.utils.regions import Sphere
from pymses.utils.regions import Box
from pymses.filters import RegionFilter
from pymses.filters import PointFunctionFilter


#from matplotlib import pyplot as P
import pylab as P
import glob as glob
import pdb
import pickle as pickle
from struct import pack

from numpy.linalg import eigh,eig


from pymses.sources.hop.file_formats import *

##from pymses.sources.hop.hop import HOP                      
#from hop_camille.hop import HOP

#import Olivier's package
#sys.path.append('../../clump/tests/')
#import clump_list as FL

#import plotting module
#import module_plot_clump as mpc



#######################################################################################
#######################################################################################
#this routine use the FOF algorithm developed by Olivier Iffrig to extract clumps defined 
#from a density threshold 
#as an output it provides a list of cell position ordered by the group to which they belong
#path_in is the path where the data tobe read are located
#path_out is the path of teh directory where resulting files must be written
#num output number
#name a string which is used to write the names of the various files
#thres_dens density threshold above which cells are considered
#thres_level level threshold above which cells are considered
#order is to treat little endian
#######################################################################################
#######################################################################################
def make_clump_fof(path_in,num,name,thres_dens,path_out=None,force=False,gcomp=True,order='<'):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    name_txt = path_in + name+'.txt'

    #check whether hop entry files have been created (not test is done on .txt only
    if ( len(glob.glob(name_txt)) == 0 or force ):
        ro=pymses.RamsesOutput(path_in,num,order=order)

        amr = ro.amr_source(["rho"], grav_compat=gcomp) #density only is used

        d_op=lambda dset: dset["rho"]

        flg, xx_v, yy_v, zz_v = FL.test_fof_print(amr, d_op, thres_dens)


        ngroups = flg.max()


        npart_v = np.zeros(ngroups,dtype=int)

        for i in range(ngroups):
            mask = flg==i+1
            npart_v[i]=np.sum(mask)

# write the sorted cells
        name_save_clump=directory_out+name + '.save'

        np.savez(name_save_clump,ngroups=ngroups,npart_v=npart_v,xx_v=xx_v,yy_v=yy_v,zz_v=zz_v,vect_id_group=flg,num=num,name=name,thres_dens=thres_dens)

#        pdb.set_trace() # debug mode     
    return

########################################################################################
########################################################################################
## make various PDF
########################################################################################
########################################################################################                 
def get_pdf(path_in,num,path_out=None,force=False,ps=False,xmin=0,xmax=1.,ymin=0,ymax=1.,zmin=0.,zmax=1.,tag=''):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in

    if(tag != ''):
        tag = tag+'_'

    name_save_pdf=directory_out+'pdf_' + tag + str(num).zfill(5) + '.save'





    if ( len(glob.glob(name_save_pdf+'*')) == 0 or force ):
        ro=pymses.RamsesOutput(path_in,num)

        lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)  

        amr = ro.amr_source(["rho","vel","P","Bl","Br"])

        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        dx = cells.get_sizes()


        mask1 = cells.points[:,0] > xmin
        mask2 = cells.points[:,0] < xmax
        mask3 = cells.points[:,1] > ymin
        mask4 = cells.points[:,1] < ymax
        mask5 = cells.points[:,2] > zmin
        mask6 = cells.points[:,2] < zmax
        mask = mask1*mask2*mask3*mask4*mask5*mask6


        pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

        vol  = dx**3
        mass = cells["rho"]*vol * scale_mass / Ms

        log_rho = np.log10(cells['rho'])
        log_min=min(log_rho)
        log_max=max(log_rho)
        nbin=50.
        width=(log_max-log_min)/nbin

        hist_logrho , histrho_edges = P.histogram(log_rho[mask],bins=nbin,range=(log_min,log_max),weights=mass[mask])
        pdf_logrho , histrhomv_edges = P.histogram(log_rho[mask],bins=nbin,range=(log_min,log_max),weights=vol[mask])


        B = 0.5 * (cells["Bl"] + cells["Br"])
        Bnorm = np.sqrt( B[:,0]**2 + B[:,1]**2 + B[:,2]**2 ) * scale_mag * microG
        pdf_B , histrho_edges = P.histogram(log_rho[mask],bins=nbin,range=(log_min,log_max),weights=vol[mask]*Bnorm[mask])
        pdf_B = pdf_B / pdf_logrho


        Temp = cells["P"] / cells["rho"] * scale_T2
        pdf_T , histrho_edges = P.histogram(log_rho[mask],bins=nbin,range=(log_min,log_max),weights=mass[mask]*Temp[mask])
        pdf_T = pdf_T / hist_logrho

        np.savez(name_save_pdf,hist_logrho=hist_logrho,pdf_logrho=pdf_logrho,histrho_edges=histrho_edges,histrhomv_edges=histrhomv_edges,pdf_B=pdf_B,pdf_T=pdf_T,width=width)


    else:
        name_save = glob.glob(name_save_pdf+'*')

        print 'restore file ',name_save

        npzfile = np.load(name_save[0])
        hist_logrho = npzfile['hist_logrho']
        pdf_logrho  = npzfile['pdf_logrho']
        histrho_edges  = npzfile['histrho_edges']
        histrhomv_edges  = npzfile['histrhomv_edges']
        pdf_B  = npzfile['pdf_B']
        pdf_T  = npzfile['pdf_T']
        width  = npzfile['width']


    #mass weighted histogram of density 
    P.clf()
    P.bar(histrho_edges[:-1],np.log10(hist_logrho),width)
#    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(\rho)$ (cm$^{-3}$)')
    P.ylabel(r'$log(M)$ (M$_\odot$)')

    if(ps):
        P.savefig(directory_out+'pdfmw_logrho_'+tag+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+'pdfmw_logrho_'+tag+str(num).zfill(5)+'.png')


    #volume weighted histogram of density 
    P.clf()
    P.bar(histrho_edges[:-1],np.log10(pdf_logrho),width)
#    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(\rho)$ (cm$^{-3}$)')
    P.ylabel(r'$PDF$')
    if(ps):
        P.savefig(directory_out+'pdfvw_logrho_'+tag+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+'pdfvw_logrho_'+tag+str(num).zfill(5)+'.png')



    #volume weighted histogram of magnetic intensities 
    P.clf()
    P.bar(histrho_edges[:-1],np.log10(pdf_B),width)
#    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(\rho)$ (cm$^{-3}$)')
    P.ylabel(r'$log(B) \, (\mu G)$')
    if(ps):
        P.savefig(directory_out+'pdf_B_'+tag+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+'pdf_B_'+tag+str(num).zfill(5)+'.png')


    #volume weighted histogram of temperature 
    P.clf()
    P.bar(histrho_edges[:-1],np.log10(pdf_T),width)
#    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(\rho)$ (cm$^{-3}$)')
    P.ylabel(r'$log(T) \, (K)$')
    if(ps):
        P.savefig(directory_out+'pdf_T_'+tag+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+'pdf_T_'+tag+str(num).zfill(5)+'.png')


#    pdb.set_trace() # debug mode     

    return


#######################################################################################
#######################################################################################
#this routine use the HOP algorithm to extract clumps defined from their peaks
#as an output it provides a list of cell position ordered by the group to which they belong
#path_in is the path where the data tobe read are located
#path_out is the path of teh directory where resulting files must be written
#num output number
#name a string which is used to write the names of the various files
#thres_dens density threshold above which cells are considered
#thres_level level threshold above which cells are considered
#pos_zoom the center of the zoom coordinates
#size_zoom the 3 zoom extension (x, y and z)
#######################################################################################
#######################################################################################
def make_clump_hop(path_in,num,name,thres_dens,thres_level,pos_zoom,size_zoom,path_out=None,path_hop='',force=False,gcomp=True):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    name_txt = name+'.txt'

    #check whether hop entry files have been created (not test is done on .txt only
    if ( len(glob.glob(name_txt)) == 0 or force ):
        ro=pymses.RamsesOutput(path_in,num)
        amr = ro.amr_source(["rho"], grav_compat=gcomp) #density only is used
        center=pos_zoom
        radius=size_zoom
        min_coords=np.zeros(3)
        max_coords=np.zeros(3)
        min_coords[:]=center[:]-radius/2.
        max_coords[:]=center[:]+radius/2.

        region = Box((min_coords,max_coords))

    #region = Sphere(center,radius)
        filt_amr = RegionFilter(region,amr)
        cell_source = CellsToPoints(filt_amr,)

    #selection of the cells of interest
        def cell_selec_func(dset):
            mask1 = dset["rho"] >= thres_dens
            dx=dset.get_sizes()
            mask2 = dx <= 0.5**thres_level
            return mask1*mask2

    #begin cell_selec
        cells_selec = PointFunctionFilter(cell_selec_func,cell_source).flatten()

        dx = cells_selec.get_sizes()

        ncells = cells_selec.npoints

     # fill the matrice with ID, x,y,z and masses of particles
        val_mat = np.zeros((ncells,5))

        val_mat[:,0]=np.arange(ncells)

        val_mat[:,1]=cells_selec.points[:,0]

        val_mat[:,2]=cells_selec.points[:,1]

        val_mat[:,3]=cells_selec.points[:,2]

        val_mat[:,4]=cells_selec["rho"]*(dx**3)

      # write name.txt
        head = str(ncells)

        np.savetxt(name_txt,val_mat,fmt='%10d %.10e %.10e %.10e %.10e',header=head,delimiter=' ',comments=' ')

#########################
#end of creation name.txt
#########################

#########################
#creation name.den
#########################
        f = open(name+'.den','wb')
        #f.write(ncells)                                             
        #rho32 = bytearray(cells_selec["rho"].astype("f"))                                                  
        #f.write(rho32)                                                                                         
        f.write(pack('I',ncells))
        cells_selec["rho"].astype("f").tofile(f)
        f.close()

        print name+'.den created'
#########################
#end of creation name.den
#########################

#########################
# HOP Algorithm
#########################
    print 'creation of .hop and .gbound du to hop'
#path= '../hop_v1.1/'
    fname = path_hop+name+'.txt'
    print 'look for hop in ',fname

    h=HOP(fname, path_hop)
    h.process_hop()                                                                                                                       
    print 'hop grouping is finished'


#########################
# end of HOP Algorithm 
#########################

    idpart=val_mat[:,0]
    X=val_mat[:,1]
    Y=val_mat[:,2]
    Z=val_mat[:,3]
    mass=val_mat[:,4]


    ##read the gbound file to get list of particle numbers within groups
    f = open(name+'.gbound','r')
    aline = f.readline()
    ngroups = int(aline)
    npart_v = np.zeros(ngroups,dtype=int)
    for i in range(10):
        aline = f.readline()
    for i in range(ngroups):
        aline = f.readline()
        vec = aline.split()
        igroup = int(vec[0])
        npart_v[igroup]=int(vec[1])
    f.close()

##get the igroup array
    group_ids = h.get_group_ids()


##sort it and apply the sorting to the coordinates
## this means that the particules of group 1 are written first then of group 2 etc... 
    ind_sort = np.argsort(group_ids)
    xx_v = X[ind_sort]
    yy_v = Y[ind_sort]
    zz_v = Z[ind_sort]
    vect_id_group = group_ids[ind_sort] #not so useful


# write the sorted cells
    name_save_clump=directory_out+name+'.save'

    np.savez(name_save_clump,ngroups=ngroups,npart_v=npart_v,xx_v=xx_v,yy_v=yy_v,zz_v=zz_v,vect_id_group=vect_id_group,num=num,name=name,thres_dens=thres_dens,thres_level=thres_level,pos_zoom=pos_zoom,size_zoom=size_zoom)

    return name_save_clump

#######################################################################################
#######################################################################################
## calculate the properties of clumps
#######################################################################################
#######################################################################################
def clump_properties(name,path_in,num,path_out=None,gcomp=True,nograv=False,order='<',log_min=None,log_max=None):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in

    name_save_clump=directory_out+name+'.save.npz'

    npzfile = np.load(name_save_clump)

    ngroups=npzfile['ngroups']
    npart_v=npzfile['npart_v']
    xx_v = npzfile['xx_v']
    yy_v = npzfile['yy_v']
    zz_v = npzfile['zz_v']
    thres_dens = npzfile['thres_dens']


#    pdb.set_trace() # debug mode     

    #cumulatif
    npart_vcumul=np.cumsum(npart_v)

    name_prop_clump = extract_clump_properties(path_in,num,name,xx_v,yy_v,zz_v,npart_vcumul,thres_dens,path_out=path_out,gcomp=gcomp,nograv=nograv,order=order)

#    mpc.plot_clump_properties(name_prop_clump,name,directory_out,num,path_out=path_out,log_min=log_min,log_max=log_max)

    return name_prop_clump

#######################################################################################
#######################################################################################
## make images of the cores
## they must have been extracted and their properties must have been calculated
## which means that extract_clump_properties must have been running first
#######################################################################################
#######################################################################################
def clump_images_series(name,path_in,num,path_out=None,gcomp=True,coldens=True,order='<',mass_thres=0.,force=False,im_core='IMAGE_CORE',im_core2='IMAGE_CORE2'):

    if( path_out is not None):
        directory_out = path_out 
    else:
        directory_out = path_in


    directory_im  = directory_out + im_core

    directory_im2 = directory_out + im_core2


    #check image directory exists
    if ( len(glob.glob(directory_im2)) == 0  ):
        str_mkdir = 'mkdir ' + directory_im2
        os.system(str_mkdir)
    directory_im2=directory_im2+'/'


    #check image directory exists
    if ( len(glob.glob(directory_im)) == 0  ):
        str_mkdir = 'mkdir ' + directory_im
        os.system(str_mkdir)
    directory_im=directory_im+'/'



    ro=pymses.RamsesOutput(path_in,num,order=order)
#    amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"], grav_compat=gcomp)
    amr = ro.amr_source(["rho","vel","P","Bl","Br"]) #, grav_compat=gcomp)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)  


    #restore particle list
#    name_save_clump=path_in+name+'.save.npz'
    #usually the particle list is saved in path_out
    name_save_clump=directory_out+name+'.save.npz'
    npzfile = np.load(name_save_clump)

    ngroups=npzfile['ngroups']
    npart_v=npzfile['npart_v']
    xx_v = npzfile['xx_v']
    yy_v = npzfile['yy_v']
    zz_v = npzfile['zz_v']
    thres_dens = npzfile['thres_dens']


    #cumulatif
    npart_vcumul=np.cumsum(npart_v)


    #restore core properties
    name_prop_clump=name+'_prop_struct_'+str(num)+'.save'
    f = open(directory_out+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()


    #now loop over cores and make some images
    for i in range(len(npart_vcumul)-1):
        if( prop_clump['ncell'][i] > 100 and prop_clump['mass'][i] > mass_thres):
            print 'treat clump ',i
            #center with respect to the peak density
            center=prop_clump['pos_n_max'][i,:] / lbox
            #spatial extension equal the double of the largest distance within the core
            radius=prop_clump['size'][i] / lbox #/ 2.

            selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
            #make_image_core expects to get the coordinate with respect to center
            x_v=xx_v[selec]-center[0]
            y_v=yy_v[selec]-center[1]

#            pdb.set_trace() # debug mode     

            if(prop_clump['mu'][i] > 1 and prop_clump['rho_mean'][i] < 1.e6 and prop_clump['n_max'][i] > 1.e5 and prop_clump['mass'][i] > 1.):



                print 'mu', prop_clump['mu'][i], 'rho_mean ', prop_clump['rho_mean'][i],' mass ',prop_clump['mass'][i] 


                make_image_core(amr,ro,center,radius,x_v,y_v,num,directory_im2,i_im=i,force=force,coldens=coldens)


                make_image(amr,ro,center,radius,num,path_in,x_v=None,y_v=None,z_v=None,mass_v=None,i_im=i,force=force,col_dens_only=False,tag='',ps=False,path_out=directory_im)


                make_image_B(amr,ro,center,radius,num,path_in,x_v=None,y_v=None,z_v=None,mass_v=None,i_im=i,force=force,tag='',ps=False,path_out=directory_im)


    return
#######################################################################################
#######################################################################################
## make images of a list of cores contains in core_list 
## they must have been extracted and their properties must have been calculated
## which means that extract_clump_properties must have been running first
#######################################################################################
#######################################################################################
def clump_images_list(name,path_in,num,core_list,path_out=None,force=False,gcomp=True):

    if( path_out is not None):
        directory_out = path_out 
    else:
        directory_out = path_in


    ro=pymses.RamsesOutput(path_in,num)
    amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"], grav_compat=gcomp)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)  



    #restore core properties
    name_prop_clump='prop_struct_'+str(num)+'.save'
    f = open(path_in+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()


    #now loop over cores and make some images
    for i in core_list:

        center=prop_clump['pos_n_max'][i,:] / lbox
        #spatial extension equal the double of the largest distance within the core
        radius=prop_clump['size'][i] / lbox 

#        make_image(amr,ro,center,radius,num,path_in,i_im=i,force=force,path_out=path_out)
        
        make_image_B(amr,ro,center,radius,num,path_in,i_im=i,force=force,path_out=path_out)

        


    return
########################################################################################
########################################################################################
## make a column density map for all clumps/cores using make_hierach_map
## only the cells belonging to the clumps/cores are displayed
## path_in : where to find the data
## num : output number
## name : the name of the numpy save file, which contains the coordinates of the clump/core cells 
## path_out : where to write the images
## nograv : gravity output format
## order : small or big endian 
## im_hier : another specific directory to put the images 
## ps : do also ps files
########################################################################################
########################################################################################                 
def make_clump_map_hierarch_series(name,path_in,num,path_out=None,gcomp=True,nograv=False,order='<',im_hier='IMAGE_HIERAR',ps=False):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in

    directory_im = directory_out + im_hier


    name_save_clump=directory_out+name+'.save.npz'

    npzfile = np.load(name_save_clump)

    ngroups=npzfile['ngroups']
    npart_v=npzfile['npart_v']
    xx_v = npzfile['xx_v']
    yy_v = npzfile['yy_v']
    zz_v = npzfile['zz_v']
    thres_dens = npzfile['thres_dens']


#    pdb.set_trace() # debug mode     

    #cumulatif
    npart_vcumul=np.cumsum(npart_v)


    #check image directory exists
    if ( len(glob.glob(directory_im)) == 0  ):
        str_mkdir = 'mkdir ' + directory_im
        os.system(str_mkdir)
    directory_im=directory_im+'/'



    ro=pymses.RamsesOutput(path_in,num,order=order)

    if (not nograv):
        amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"], grav_compat=gcomp)
    else:
        amr = ro.amr_source(["rho","vel","P","Bl","Br"])

    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)  


    #now query the cells using their positions
    #due to memory issue in the cpu communication split it in two (enough for the present purpose)
    ntot = len(xx_v)
    ntot_2=ntot/2


    points = np.concatenate((xx_v[0:ntot_2],yy_v[0:ntot_2],zz_v[0:ntot_2])).reshape(3,len(xx_v[0:ntot_2])).T
    psamp = PointSamplingProcessor(amr)
    cells1 = psamp.process(points,add_level=True)


    points = np.concatenate((xx_v[ntot_2:ntot],yy_v[ntot_2:ntot],zz_v[ntot_2:ntot])).reshape(3,len(xx_v[ntot_2:ntot])).T
    psamp = PointSamplingProcessor(amr)
    cells2 = psamp.process(points,add_level=True)


#    pdb.set_trace() # debug mode     

    #now merge the two structures
    cells=cells1.concatenate([cells1,cells2])

    dx = 0.5**cells["level"]


    #call for normalisation 
    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)


    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])

        rho_v = cells['rho'][selec]
        x_v = xx_v[selec]
        y_v = yy_v[selec]
        z_v = zz_v[selec]
        dx_v = dx[selec]

        xmax = np.max(x_v)
        xmin = np.min(x_v)

        ymax = np.max(y_v)
        ymin = np.min(y_v)

        zmax = np.max(z_v)
        zmin = np.min(z_v)

        ind_max=np.argmax(cells["rho"][selec])

        xref = x_v[ind_max]
        yref = y_v[ind_max]
        zref = z_v[ind_max]

        center_max=[xref,yref,zref]

        size=np.max([xmax-xref,xref-xmin,ymax-yref,yref-ymin,zmax-zref,zref-zmin])

        nzoom = np.fix(np.log(size)/np.log(0.5))

        eps=0.5**nzoom
        map_coldens , map_w13, xedges, yedges = make_hierarch_map(x_v,y_v,z_v,dx_v,rho_v,rho_v,eps,center=center_max)

        map_coldens = map_coldens.T 

        xedges = xedges * lbox
        yedges = yedges * lbox

        mask = map_coldens != 0.
#        map_w13[mask] = map_w13[mask] / map_coldens[mask]

        map_coldens = map_coldens* unit_col

        log_map_coldens = map_coldens
        log_map_coldens[mask] = np.log10(log_map_coldens[mask])


        P.clf()
        xmin = (center_max[0]-eps)*lbox
        xmax = (center_max[0]+eps)*lbox
        ymin = (center_max[1]-eps)*lbox
        ymax = (center_max[1]+eps)*lbox

        im = P.imshow(log_map_coldens,origin='lower',extent=[xmin,xmax,ymin,ymax] )

        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(N) (cm^{-2})$')

        P.xlabel(r'$x$ (pc)')     
        P.ylabel(r'$y$ (pc)')

        if(ps):
            P.savefig(directory_im+'/coldens_z_hm_'+'_'+str(num)+'_'+str(i)+'.ps')
        P.savefig(directory_im+'/coldens_z_hm'+'_'+str(num)+'_'+str(i)+'.png')


    return
########################################################################################
########################################################################################
## extract the properties of the clumps/cores
########################################################################################
########################################################################################                 
def extract_clump_properties(path_in,num,name,xx_v,yy_v,zz_v,npart_vcumul,thres_dens,path_out=None,gcomp=True,nograv=False,order='<'):




    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    ro=pymses.RamsesOutput(path_in,num,order=order)

    if (not nograv):
        amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"], grav_compat=gcomp)
    else:
        amr = ro.amr_source(["rho","vel","P","Bl","Br"])

    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)  





    #now query the cells using their positions
    #due to memory issue in the cpu communication split it in two (enough for the present purpose)
    ntot = len(xx_v)
    ntot_2=ntot/2


#    pdb.set_trace() # debug mode     

    points = np.concatenate((xx_v[0:ntot_2],yy_v[0:ntot_2],zz_v[0:ntot_2])).reshape(3,len(xx_v[0:ntot_2])).T
    psamp = PointSamplingProcessor(amr)
    cells1 = psamp.process(points,add_level=True)

#    pdb.set_trace() # debug mode     

    points = np.concatenate((xx_v[ntot_2:ntot],yy_v[ntot_2:ntot],zz_v[ntot_2:ntot])).reshape(3,len(xx_v[ntot_2:ntot])).T
    psamp = PointSamplingProcessor(amr)
    cells2 = psamp.process(points,add_level=True)


#    pdb.set_trace() # debug mode     

    #now merge the two structures
    cells=cells1.concatenate([cells1,cells2])


#    pdb.set_trace() # debug mode     

    dx = 0.5**cells["level"]
    dV = dx**3

    res_v=np.zeros(len(npart_vcumul)-1)

    res3_v=np.zeros((len(npart_vcumul)-1,3))


    res6_v=np.zeros((len(npart_vcumul)-1,3,3))


    #initiate dictionaries
    #dictionnary to contain the structure properties
    prop_clump=dict()
    #dictionnary to describe the structure properties (description, units, query index for database, show index for database)
    #query= 0 do not query on this parameter, 1 do it. 
    #show = 0 low priority, 2 high priority
    head_clump=dict()

    #dictionnary for COAST database. Contains both data and header information
#    prop_clump_CDB=dict()


    #initialise the inertia matrix
    mat_iner=np.zeros((3,3))


    #call for normalisation 
    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)
    G=6.7e-8 # G is used to calculate the mass to flux ratio

    #number of cells within the cores
    describ='number of cells within the cores'
    units='none'
    query=1
    show =1
    for i in range(len(npart_vcumul)-1):
        res_v[i] = npart_vcumul[i+1] - npart_vcumul[i]

    prop_clump.update({'ncell':np.copy(res_v)})
    head_clump.update({'ncell':[describ,units,query,show]})

#    field='ncell'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mass of the cores
    describ='mass of the core'
    units='Msun'
    query=1
    show =2
    for i in range(len(npart_vcumul)-1):
        mass_code = np.sum(cells["rho"][npart_vcumul[i]:npart_vcumul[i+1]]*dV[npart_vcumul[i]:npart_vcumul[i+1]])
        res_v[i]=mass_code

    mass_v = np.copy(res_v)
    res_v=res_v * (scale_l*lbox)**3 * scale_d / Ms 
    prop_clump.update({'mass':np.copy(res_v)})
    head_clump.update({'mass':[describ,units,query,show]})

#    field='mass'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the volume of the cores (pc^-3)
    describ='volume of the core'
    units='pc^3'
    query=0
    show =0
    for i in range(len(npart_vcumul)-1):
        vol_code = np.sum(dV[npart_vcumul[i]:npart_vcumul[i+1]])
        res_v[i]=vol_code

    vol_v = np.copy(res_v)
    res_v=res_v * (lbox*scale_l/pc)**3 
    prop_clump.update({'vol':np.copy(res_v)})
    head_clump.update({'vol':[describ,units,query,show]})

#    field='vol'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    #mean density of the cores
    describ='mean density of the core'
    units='cm^-3'
    query=1
    show =1
    res_v = mass_v / vol_v
    prop_clump.update({'rho_mean':np.copy(res_v)})
    head_clump.update({'rho_mean':[describ,units,query,show]})    

#    field='rho_mean'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mass above various density thresholds
    units='Msun'
    query=1
    show =1
    i_thr=0
    thres_v=[3.e3,1.e4,3.e4,1.e5,3.e5,1.e6,3.e6,1.e7]
    for thres in thres_v:
        describ='mass above density threshold of '+ str(thres)
        for i in range(len(npart_vcumul)-1):
            mask = cells["rho"][npart_vcumul[i]:npart_vcumul[i+1]] > thres
            selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])[mask]
            mass_code = np.sum(cells["rho"][selec]*dV[selec])
            res_v[i]=mass_code * (scale_l*lbox)**3 * scale_d / Ms 


        mass_thr = 'mass_'+str(i_thr)
        prop_clump.update({mass_thr:np.copy(res_v)})
        head_clump.update({mass_thr:[describ,units,query,show]})

#        field=mass_thr
#        prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

        i_thr=i_thr+1


    #compute the peak density of the cores and its position (pc)
    for i in range(len(npart_vcumul)-1):
        res_v[i]=np.max(cells["rho"][npart_vcumul[i]:npart_vcumul[i+1]])

        ind_max=np.argmax(cells["rho"][npart_vcumul[i]:npart_vcumul[i+1]])

        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])[ind_max] 
        res3_v[i,0] = xx_v[selec]
        res3_v[i,1] = yy_v[selec]
        res3_v[i,2] = zz_v[selec]
        
    res3_v = res3_v*lbox*scale_l/pc
    prop_clump.update({'n_max':np.copy(res_v)})
    describ='maximum density'
    units='cm^-3'
    query=1
    show =1
    head_clump.update({'n_max':[describ,units,query,show]})

#    field='n_max'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    prop_clump.update({'pos_n_max':np.copy(res3_v)})
    describ='position of maximum density'
    units='pc'
    head_clump.update({'pos_n_max':[describ,units,query,show]})

#    field='pos_n_max'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mean velocity of the cores (km/s)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        vel_code = np.sum(cells["rho"][selec]*cells["vel"][:,0][selec]*dV[selec])
        res3_v[i,0] = vel_code 
        vel_code = np.sum(cells["rho"][selec]*cells["vel"][:,1][selec]*dV[selec])
        res3_v[i,1] = vel_code 
        vel_code = np.sum(cells["rho"][selec]*cells["vel"][:,2][selec]*dV[selec])
        res3_v[i,2] = vel_code 

    res3_v[:,0] = res3_v[:,0] / mass_v
    res3_v[:,1] = res3_v[:,1] / mass_v
    res3_v[:,2] = res3_v[:,2] / mass_v
    vG = np.copy(res3_v)

    res3_v = res3_v * scale_v / km_s 
    prop_clump.update({'vel':np.copy(res3_v)})
    describ='mean velocity of the core'
    units='km/s'
    query=0
    show =0
    head_clump.update({'vel':[describ,units,query,show]})

#    field='vel'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the specific anuglar momentum of the cores (pc km/s)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        vel_code = np.sum(cells["rho"][selec]*(yy_v[selec]*(cells["vel"][:,2][selec]-vG[i,2])-zz_v[selec]*(cells["vel"][:,1][selec]-vG[i,1]))*dV[selec])
        res3_v[i,0] = vel_code 
        vel_code = np.sum(cells["rho"][selec]*(zz_v[selec]*(cells["vel"][:,0][selec]-vG[i,0])-xx_v[selec]*(cells["vel"][:,2][selec]-vG[i,2]))*dV[selec])
        res3_v[i,1] = vel_code 
        vel_code = np.sum(cells["rho"][selec]*(xx_v[selec]*(cells["vel"][:,1][selec]-vG[i,1])-yy_v[selec]*(cells["vel"][:,0][selec]-vG[i,0]))*dV[selec])
        res3_v[i,2] = vel_code 

    res3_v[:,0] = res3_v[:,0] / mass_v
    res3_v[:,1] = res3_v[:,1] / mass_v
    res3_v[:,2] = res3_v[:,2] / mass_v


    res3_v = res3_v * scale_v / km_s * lbox * scale_l / pc
    prop_clump.update({'j_v':np.copy(res3_v)})
    describ='3D specific angular momentum of the core'
    units='pc km/s'
    query=0
    show =0
    head_clump.update({'j_v':[describ,units,query,show]})

#    field='j_v'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    res_v = np.sqrt(res3_v[:,0]**2+res3_v[:,1]**2+res3_v[:,2]**2)
    prop_clump.update({'j':np.copy(res_v)})
    describ='norm of specific angular momentum of the core'
    units='pc km/s'
    query=1
    show =1
    head_clump.update({'j':[describ,units,query,show]})

#    field='j'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}



    #compute the velocity dispersion of the cores with respect to mean velocity (km/s)
    #compute also the total inner kinetic energy (code units)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        one_v = np.ones(npart_vcumul[i+1]-npart_vcumul[i])
        vel_code = np.sum(cells["rho"][selec]*((cells["vel"][:,0][selec]-vG[i,0]*one_v)**2)*dV[selec])
        res3_v[i,0] = vel_code 
        vel_code = np.sum(cells["rho"][selec]*((cells["vel"][:,1][selec]-vG[i,1]*one_v)**2)*dV[selec])
        res3_v[i,1] = vel_code 
        vel_code = np.sum(cells["rho"][selec]*((cells["vel"][:,2][selec]-vG[i,2]*one_v)**2)*dV[selec])
        res3_v[i,2] = vel_code 



    res_v = ( res3_v[:,0] + res3_v[:,1] + res3_v[:,2] ) / 2.
    prop_clump.update({'ener_kin':np.copy(res_v)})
    describ='internal kinetic energy of the core'
    units='code units'
    query=0
    show =0
    head_clump.update({'ener_kin':[describ,units,query,show]})

#    field='ener_kin'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    res3_v[:,0] = np.sqrt( res3_v[:,0] / mass_v )
    res3_v[:,1] = np.sqrt( res3_v[:,1] / mass_v )
    res3_v[:,2] = np.sqrt( res3_v[:,2] / mass_v )
    res3_v = res3_v * scale_v / km_s 

    res_v = np.sqrt( res3_v[:,0]**2 + res3_v[:,1]**2 + res3_v[:,2]**2 ) 
    prop_clump.update({'sig_v':np.copy(res3_v)})
    describ='3D internal kinetic velocity  of the core'
    units='km/s'
    query=0
    show =0
    head_clump.update({'sig_v':[describ,units,query,show]})

#    field='sig_v'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    prop_clump.update({'sig':np.copy(res_v)})
    describ='internal kinetic velocity  of the core'
    units='km/s'
    query=1
    show =2
    head_clump.update({'sig':[describ,units,query,show]})

#    field='sig'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mean position of the cores (pc)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        pos_code = np.sum(cells["rho"][selec]*xx_v[selec]*dV[selec])
        res3_v[i,0] = pos_code 
        pos_code = np.sum(cells["rho"][selec]*yy_v[selec]*dV[selec])
        res3_v[i,1] = pos_code 
        pos_code = np.sum(cells["rho"][selec]*zz_v[selec]*dV[selec])
        res3_v[i,2] = pos_code 

    res3_v[:,0] = res3_v[:,0] / mass_v
    res3_v[:,1] = res3_v[:,1] / mass_v
    res3_v[:,2] = res3_v[:,2] / mass_v
    posG = res3_v
    res3_v = res3_v * lbox * scale_l / pc
    prop_clump.update({'pos':np.copy(res3_v)})
    describ='position of barycenter'
    units='pc'
    query=0
    show =0
    head_clump.update({'pos':[describ,units,query,show]})

#    field='pos'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    #compute the length from inertia momentum (pc)
    #compute and store the eigenvectors
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        Ixx = np.sum(cells["rho"][selec]*(xx_v[selec]-posG[i,0])**2*dV[selec])
        Iyy = np.sum(cells["rho"][selec]*(yy_v[selec]-posG[i,1])**2*dV[selec])
        Izz = np.sum(cells["rho"][selec]*(zz_v[selec]-posG[i,2])**2*dV[selec])
        mat_iner[0,0] = Iyy+Izz
        mat_iner[1,1] = Ixx+Izz
        mat_iner[2,2] = Ixx+Iyy

        mat_iner[0,1] = -np.sum(cells["rho"][selec]*(xx_v[selec]-posG[i,0])*(yy_v[selec]-posG[i,1])*dV[selec])
        mat_iner[0,2] = -np.sum(cells["rho"][selec]*(xx_v[selec]-posG[i,0])*(zz_v[selec]-posG[i,2])*dV[selec])
        mat_iner[1,2] = -np.sum(cells["rho"][selec]*(yy_v[selec]-posG[i,1])*(zz_v[selec]-posG[i,2])*dV[selec])

        mat_iner[1,0] = mat_iner[0,1]
        mat_iner[2,0] = mat_iner[0,2]
        mat_iner[2,1] = mat_iner[1,2]

        if (mass_v[i] > 0. and prop_clump['ncell'][i]>10):
            w,v=eigh(mat_iner)
            w = np.sort(w)

#            res3_v[i,0] = np.sqrt(np.abs(w[0]) / mass_v[i])
#            res3_v[i,1] = np.sqrt(np.abs(w[1]) / mass_v[i])
#            res3_v[i,2] = np.sqrt(np.abs(w[2]) / mass_v[i])

            res3_v[i,0] = np.sqrt(w[0] / mass_v[i])
            res3_v[i,1] = np.sqrt(w[1] / mass_v[i])
            res3_v[i,2] = np.sqrt(w[2] / mass_v[i])

            w1 = 0.5*np.sum(w) - w
            w1 = np.sort(w1)
            if(w[0] < 0. or w1[0] < 0.):
               print 'negative eigen value', 'ncell ', prop_clump['ncell'][i], w[0],w[1],w[2]

#               pdb.set_trace() # debug mode     


            res6_v[i,:,:] = v
        else:
            res3_v[i,0] = dx[i]
            res3_v[i,1] = dx[i]
            res3_v[i,2] = dx[i]

    res3_v = res3_v * lbox * scale_l / pc
    prop_clump.update({'size_iner':np.copy(res3_v)})
    describ='3 sizes from inertia momentum tensor'
    units='pc'
    query=0
    show =0
    head_clump.update({'size_iner':[describ,units,query,show]})

#    field='size_iner'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    prop_clump.update({'dir_iner':np.copy(res6_v)})
    describ='3 eigen arrays of inertia momentum tensor (sqrt(I/M))'
    units='code units'
    query=0
    show =0
    head_clump.update({'dir_iner':[describ,units,query,show]})

#    field='dir_iner'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    res3_v = res3_v**2
    tot = 0.5*(res3_v[:,0]+res3_v[:,1]+res3_v[:,2])
    res3_0_v = np.copy(res3_v)
    #take abs because in some rare occasions this can be negative
    res3_v[:,0] = np.sqrt(abs(tot - res3_0_v[:,2]))
    res3_v[:,1] = np.sqrt(abs(tot - res3_0_v[:,1]))
    res3_v[:,2] = np.sqrt(abs(tot - res3_0_v[:,0]))
    prop_clump.update({'size_iner2':np.copy(res3_v)})
    describ='3 sizes from inertia momentum tensor (sqrt( (I1+I2+I3)/2 -Ii)'
    units='pc'
    query=1
    show =1
    head_clump.update({'size_iner2':[describ,units,query,show]})

#    field='size_iner'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}





    #compute the maximum spatial extension of the cores (pc)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        lx = np.max(xx_v[selec]) - np.min(xx_v[selec])
        ly = np.max(yy_v[selec]) - np.min(yy_v[selec])
        lz = np.max(zz_v[selec]) - np.min(zz_v[selec])

        res_v[i] = np.max([lx,ly,lz]) * lbox * scale_l / pc

    prop_clump.update({'size':np.copy(res_v)})
    describ='maximum extension of cores'
    units='pc'
    query=0
    show =0
    head_clump.update({'size':[describ,units,query,show]})

#    field='size'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mean magnetic field of the cores (given in microGauss)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        pos_code = np.sum(cells["Br"][:,0][selec]*dV[selec])
        res3_v[i,0] = pos_code 
        pos_code = np.sum(cells["Br"][:,1][selec]*dV[selec])
        res3_v[i,1] = pos_code 
        pos_code = np.sum(cells["Br"][:,2][selec]*dV[selec])
        res3_v[i,2] = pos_code 

        pos_code = np.sum(((cells["Br"][:,0][selec])**2)*dV[selec])
        res_v[i] = pos_code 
        pos_code = np.sum(((cells["Br"][:,1][selec])**2)*dV[selec])
        res_v[i] += pos_code 
        pos_code = np.sum(((cells["Br"][:,2][selec])**2)*dV[selec])
        res_v[i] += pos_code 


    res3_v[:,0] = res3_v[:,0] / vol_v
    res3_v[:,1] = res3_v[:,1] / vol_v
    res3_v[:,2] = res3_v[:,2] / vol_v
    res3_v = res3_v * scale_mag * microG
    prop_clump.update({'B':np.copy(res3_v)})
    describ='mean magnetic field'
    units='uG'
    query=1
    show =1
    head_clump.update({'B':[describ,units,query,show]})

#    field='B'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    res_v = res_v / vol_v * (scale_d*scale_v**2) / 2.
    prop_clump.update({'PB':np.copy(res_v)})
    describ='mean magnetic pressure'
    units='erg cm^-3'
    query=0
    show =0
    head_clump.update({'PB':[describ,units,query,show]})

#    field='PB'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mean volume weighted Alfven speed
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        pos_code = np.sum(np.sqrt((cells["Br"][:,0][selec]**2+cells["Br"][:,1][selec]**2+cells["Br"][:,2][selec]**2)/cells["rho"][selec])*dV[selec])
        res_v[i]=pos_code

    res_v = res_v / vol_v *scale_v / km_s
    prop_clump.update({'Va':np.copy(res_v)})
    describ='mean volume weighted Alfven speed'
    units='km/s'
    query=0
    show =0
    head_clump.update({'Va':[describ,units,query,show]})

#    field='Va'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mean volume weighted  Alfvenic Mach number
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        one_v = np.ones(npart_vcumul[i+1]-npart_vcumul[i])
        sig2 = (cells["vel"][:,0][selec]-vG[i,0]*one_v)**2+(cells["vel"][:,1][selec]-vG[i,1]*one_v)**2+(cells["vel"][:,2][selec]-vG[i,2]*one_v)**2
        pos_code = np.sum(np.sqrt(sig2*cells["rho"][selec]/(cells["Br"][:,0][selec]**2+cells["Br"][:,1][selec]**2+cells["Br"][:,2][selec]**2))*dV[selec])
        res_v[i]=pos_code

    res_v = res_v / vol_v
    prop_clump.update({'Mach_alfv':np.copy(res_v)})
    describ='mean volume weighted Alfvenic Mach number'
    units='none'
    query=1
    show =1
    head_clump.update({'Mach_alfv':[describ,units,query,show]})

#    field='Mach_alfv'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mean volume weighted  Mach number
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        one_v = np.ones(npart_vcumul[i+1]-npart_vcumul[i])
        sig2 = (cells["vel"][:,0][selec]-vG[i,0]*one_v)**2+(cells["vel"][:,1][selec]-vG[i,1]*one_v)**2+(cells["vel"][:,2][selec]-vG[i,2]*one_v)**2
        pos_code = np.sum(np.sqrt(cells["rho"][selec]*sig2/cells["P"][selec])*dV[selec])
        res_v[i]=pos_code

    res_v = res_v / vol_v
    prop_clump.update({'Mach':np.copy(res_v)})
    describ='mean volume weighted Mach number'
    units='none'
    query=1
    show =2
    head_clump.update({'Mach':[describ,units,query,show]})

#    field='Mach'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mean mass weighted  Mach number
    #the rho^3 is there to i) produce temperature (P/rho) then to weight for the mass (^2 beacuse inside sqrt)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        one_v = np.ones(npart_vcumul[i+1]-npart_vcumul[i])
        sig2 = (cells["vel"][:,0][selec]-vG[i,0]*one_v)**2+(cells["vel"][:,1][selec]-vG[i,1]*one_v)**2+(cells["vel"][:,2][selec]-vG[i,2]*one_v)**2
        pos_code = np.sum(np.sqrt((cells["rho"][selec]**3)*sig2/cells["P"][selec])*dV[selec])
        res_v[i]=pos_code

    res_v = res_v / mass_v
    prop_clump.update({'Mach_mw':np.copy(res_v)})
    describ='mean mass weighted Mach number'
    units='none'
    query=0
    show =0
    head_clump.update({'Mach_mw':[describ,units,query,show]})

#    field='Mach_mw'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the mass to magnetic flux ratio
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        x_max = prop_clump['pos_n_max'][i,0]/lbox
        y_max = prop_clump['pos_n_max'][i,1]/lbox
        z_max = prop_clump['pos_n_max'][i,2]/lbox

        mask = abs(xx_v[selec]-x_max) < 0.51*dx[selec]
        flux_x = abs(np.sum(cells["Br"][:,0][selec][mask]*(dx[selec][mask])**2))

        mask = abs(yy_v[selec]-y_max) < 0.51*dx[selec]
        flux_y = abs(np.sum(cells["Br"][:,1][selec][mask]*(dx[selec][mask])**2))

        mask = abs(zz_v[selec]-z_max) < 0.51*dx[selec]
        flux_z = abs(np.sum(cells["Br"][:,2][selec][mask]*(dx[selec][mask])**2))

        res_v[i] = 1./np.max([flux_x,flux_y,flux_z]) 

    res_v = res_v * mass_v
    #  / (0.13/np.sqrt(G)) normalisation according to Mouschovias & Spitzer
    res_v = res_v / scale_mag * (pc*lbox*scale_d)  / (0.13/np.sqrt(G))
    prop_clump.update({'mu':np.copy(res_v)})
    describ='mass to flux over critical mass to flux'
    units='none'
    query=1
    show =2
    head_clump.update({'mu':[describ,units,query,show]})

#    field='mu'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the alpha virial parameter based on radius obtained from total volume
    mass_vlog = np.log10(prop_clump['mass'])
    rad = (3./4./np.pi*prop_clump['vol'])**(0.3333) #an estimate of the radius
    res_v = 5./3.*(prop_clump['sig']*km_s)**2 / (G*prop_clump['mass']*Ms/(rad*pc))
    prop_clump.update({'alpha_vir':np.copy(res_v)})
    describ='alpha virial parameter based on radius obtained from total volume'
    units='none'
    query=0
    show =0
    head_clump.update({'alpha_vir':[describ,units,query,show,query,show]})

#    field='alpha_vir'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the alpha virial parameter based on radius obtained from inertia momentum
    mass_vlog = np.log10(prop_clump['mass'])
    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    res_v = 5./3.*(prop_clump['sig']*km_s)**2 / (G*prop_clump['mass']*Ms/(rad*pc))
    prop_clump.update({'alpha_vir1':np.copy(res_v)})
    describ='alpha virial parameter based on radius obtained from inertia momentum'
    units='none'
    query=0
    show =1
    head_clump.update({'alpha_vir1':[describ,units,query,show]})

#    field='alpha_vir1'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #compute the gravitational virial term (code units)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        if( not nograv):
            grav_vir = np.sum( cells["rho"][selec]*((xx_v[selec]-posG[i,0])*cells["g"][:,0][selec]+(yy_v[selec]-posG[i,1])*cells["g"][:,1][selec]+(zz_v[selec]-posG[i,2])*cells["g"][:,2][selec])*dV[selec])
        #we must multiply it by lbox because the gradient contains 1/lbox in ramses output
        else:
            grav_vir=0.
        res_v[i]=grav_vir*lbox

    prop_clump.update({'grav_vir':np.copy(res_v)})
    describ='gravitational virial parameter'
    units='code units'
    query=1
    show =2
    head_clump.update({'grav_vir':[describ,units,query,show]})

#    field='grav_vir'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    #compute the thermal energy virial term (code units)
    for i in range(len(npart_vcumul)-1):
        selec = np.arange(npart_vcumul[i],npart_vcumul[i+1])
        therm_vir = np.sum(cells["P"][selec]*dV[selec])*3./2.
        res_v[i]=therm_vir


    prop_clump.update({'therm_vir':np.copy(res_v)})
    describ='thermal virial parameter'
    units='code units'
    query=0
    show =0
    head_clump.update({'therm_vir':[describ,units,query,show]})

#    field='therm_vir'
#    prop_clump_CDB[field] = {'data': prop_clump[field], 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}


    #save the dictionary
    name_prop_clump=name+'_prop_struct_'+str(num)+'.save'

    print 'write ',name_prop_clump, ' in ',directory_out

    f = open(directory_out+name_prop_clump,'w')
    pickle.dump(prop_clump,f)
    pickle.dump(thres_v,f)
    pickle.dump(lbox,f)
    pickle.dump(head_clump,f)
    f.close()


    #save the COASTDB dictionary
#    name_prop_clump=name+'_prop_struct_CDB'+str(num)+'.save'

#    print 'write ',name_prop_clump, ' in ',directory_out

#    f = open(directory_out+name_prop_clump,'w')
#    pickle.dump(prop_clump_CDB,f)
#    pickle.dump(thres_v,f)
#    pickle.dump(lbox,f)
#    f.close()




#    pdb.set_trace() # debug mode     


    return name_prop_clump


########################################################################################
########################################################################################
## modify the shape of the clump catalog
########################################################################################
########################################################################################                 
def reshape_prop_clump(name_prop_clump,name,path):

    #load the catalog
    f = open(path+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    lbox = pickle.load(f)
    head_clump = pickle.load(f)
    f.close()


    #dictionnary for COAST database. Contains both data and header information
    prop_clump_CDB=dict()


    for key in prop_clump:

         data = prop_clump[key]
         describ = head_clump[key][0]
         units = head_clump[key][1]
         query = head_clump[key][2]
         show = head_clump[key][3]


         prop_clump_CDB[key] =  {'data':data, 'descr': describ, 'filter_flag': query, 'sort_flag': show, 'unit': units}

    name_prop_clump=name+'_prop_struct_CDB.save'

    print 'write ',name_prop_clump, ' in ',path

    f = open(path+name_prop_clump,'w')
    pickle.dump(prop_clump_CDB,f)
    pickle.dump(thres_v,f)
    pickle.dump(lbox,f)
    f.close()

    return

########################################################################################
########################################################################################
## establish  clump correlation
########################################################################################
########################################################################################                 
def correl_clump(name_prop_clump,name,path_in,num,path_out=None):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


#    if(subname != ''):
#        subname='_'+subname+'_'


    f = open(directory_out+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()

    ##dictionnary with clump correlation properties
    correl_clump=dict()



            #ordered the clumps with respect to x-variable in order to speed up the  
#        for key in prop_clump.keys():
#            prop_clump[key] = prop_clump[key][mask]



    #loop over clumps and calculate various statistics
    Nstneig = 10
    n_neig = np.zeros(Nstneig,dtype=np.int)
        
    n_neig[0:5] = np.arange(5)+3
    n_neig[5]   = 10
    n_neig[6]   = 20
    n_neig[7]   = 50
    n_neig[8]   = 100
    n_neig[9]   = 1000

    nclump=len(prop_clump['mass'])

    ind_first_neig = np.zeros((nclump,Nstneig),dtype=np.int)
    M_neig = np.zeros((nclump,Nstneig))
    rad_neig = np.zeros((nclump,Nstneig))
    AR_neig  = np.zeros((nclump,Nstneig))
    AR2_neig = np.zeros((nclump,Nstneig))

    ind_first_sgrav_neig = np.zeros(nclump,dtype=np.int)


    #initialise the inertia matrix
    mat_iner=np.zeros((3,3))

    for i in range(nclump):

        pos_rel = prop_clump['pos_n_max']
        pos_rel[:,0] = pos_rel[:,0] - prop_clump['pos_n_max'][i,0]
        pos_rel[:,1] = pos_rel[:,1] - prop_clump['pos_n_max'][i,1]
        pos_rel[:,2] = pos_rel[:,2] - prop_clump['pos_n_max'][i,2]

        rad = np.sqrt(pos_rel[:,0]**2+pos_rel[:,1]**2+pos_rel[:,2]**2) 

        ind_sort = np.argsort(rad)

        ind_first_neig[i,:] = ind_sort[1:Nstneig+1]


        for iineig in range(6):
            ineig = iineig + 4
            nneig= n_neig[ineig]
            selec = ind_sort[1:nneig+1]
            M_neig[i,ineig] = np.sum(prop_clump['mass'][selec])
            rad_neig[i,ineig] = np.sum(rad[selec] * prop_clump['mass'][selec])
            rad_neig[i,ineig] = rad_neig[i,ineig] / M_neig[i,ineig]


            Ixx = np.sum(prop_clump['mass'][selec]*(pos_rel[:,0][selec])**2)
            Iyy = np.sum(prop_clump['mass'][selec]*(pos_rel[:,1][selec])**2)
            Izz = np.sum(prop_clump['mass'][selec]*(pos_rel[:,2][selec])**2)

            mat_iner[0,0] = Iyy+Izz
            mat_iner[1,1] = Ixx+Izz
            mat_iner[2,2] = Ixx+Iyy

            mat_iner[0,1] = -np.sum(prop_clump['mass'][selec]*pos_rel[:,0][selec]*pos_rel[:,1][selec])
            mat_iner[0,2] = -np.sum(prop_clump['mass'][selec]*pos_rel[:,0][selec]*pos_rel[:,2][selec])
            mat_iner[1,2] = -np.sum(prop_clump['mass'][selec]*pos_rel[:,1][selec]*pos_rel[:,2][selec])

            mat_iner[1,0] = mat_iner[0,1]
            mat_iner[2,0] = mat_iner[0,2]
            mat_iner[2,1] = mat_iner[1,2]


            w,v=eigh(mat_iner)
            w = np.sum(w)/2. - w
            w = np.sort(w)            

#            AR_neig[i,ineig]  = np.sqrt(w[1]/w[2])
            AR_neig[i,ineig] = np.sqrt(np.max((w[1],0.))/w[2])
            AR2_neig[i,ineig] = np.sqrt(np.max((w[0],0.))/w[2])

#            AR2_neig[i,ineig] = np.sqrt(w[0]/w[2])

#            if(w[0] < 0.):
#                print 'negative eigen value ',w[0],w[0]/w[2],w[1]/w[2]
#                pdb.set_trace() # debug mode     

        #look for the nearest neighbour which is self-gravitating 
#        mask = prop_clump['mu'][ind_sort][1:n_neig[9]+1] > 1.
        for ineig in range(n_neig[9]):
            if(prop_clump['mu'][ind_sort][ineig+1] > 1.):
                ind_first_sgrav_neig[i] = ind_sort[ineig+1]
                break


    correl_clump.update({'ind_first_neig':np.copy(ind_first_neig)})

#    describ='gravitational virial parameter'
#    units='code units'
#    query=1
#    show =2
#    head_clump.update({'grav_vir':[describ,units,query,show]})


    correl_clump.update({'ind_first_sgrav_neig':np.copy(ind_first_sgrav_neig)})

    correl_clump.update({'n_neig':np.copy(n_neig)})

    correl_clump.update({'M_neig':np.copy(M_neig)})

    correl_clump.update({'rad_neig':np.copy(rad_neig)})

    correl_clump.update({'AR_neig':np.copy(AR_neig)})

    correl_clump.update({'AR2_neig':np.copy(AR2_neig)})


    #save the dictionary
    name_correl_clump=name+'_correl_struct_'+str(num)+'.save'
    f = open(directory_out+name_correl_clump,'w')
    pickle.dump(Nstneig,f)
    pickle.dump(n_neig,f)
    pickle.dump(correl_clump,f)
    f.close()


    return  name_correl_clump



########################################################################################
########################################################################################
## calculate the properties of clump correlation
########################################################################################
########################################################################################                 
def prop_correl_clump(name_prop_clump,name_correl_clump,name,path_in,num,path_out=None):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


#    if(subname != ''):
#        subname='_'+subname+'_'


    ##dictionnary with clump correlation properties
    prop_correl_clump=dict()


    f = open(directory_out+name_correl_clump,'r')
    Nstneig = pickle.load(f)
    n_neig = pickle.load(f)
    correl_clump = pickle.load(f)
    f.close()



    f = open(directory_out+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()

    nclump=len(prop_clump['mass'])

#    pdb.set_trace() # debug mode     

    #this mass contains the mass of the self-gravitating clumps to which the mass of the non-self-gravitating clumps 
    #for which the sg clumps is their first sg neighbours
    mass_clust_tot = np.copy(prop_clump['mass'])

    #this mass contains the mass of the self-gravitating clumps to which the mass of the non-self-gravitating clumps 
    #for which the sg clumps is their first sg neighbours and "close enough" has been added 
    mass_clust_abs1 = np.copy(prop_clump['mass'])
    mass_clust_abs2 = np.copy(prop_clump['mass'])

    mass_clust_rel1 = np.copy(prop_clump['mass'])
    mass_clust_rel2 = np.copy(prop_clump['mass'])

    for i in range(nclump):
        ind_fsg = correl_clump['ind_first_sgrav_neig'][i]

        rad = np.sqrt(np.sum((prop_clump['pos_n_max'][i,:]-prop_clump['pos_n_max'][ind_fsg,:])**2))

        if(prop_clump['mu'][i] < 1.):
            mass_clust_tot[ind_fsg] += prop_clump['mass'][i]

        if(prop_clump['mu'][i] < 1. and rad < 0.01):
            mass_clust_abs1[ind_fsg] += prop_clump['mass'][i]

        if(prop_clump['mu'][i] < 1. and rad < 0.05):
            mass_clust_abs2[ind_fsg] += prop_clump['mass'][i]

        if(prop_clump['mu'][i] < 1. and rad < 2.*np.max(prop_clump['size_iner'][i])):
            mass_clust_rel1[ind_fsg] += prop_clump['mass'][i]

        if(prop_clump['mu'][i] < 1. and rad < 5.*np.max(prop_clump['size_iner'][i])):
            mass_clust_rel2[ind_fsg] += prop_clump['mass'][i]



    prop_correl_clump.update({'M_cl_tot':np.copy(mass_clust_tot)})

    prop_correl_clump.update({'M_cl_abs1':np.copy(mass_clust_abs1)})

    prop_correl_clump.update({'M_cl_abs2':np.copy(mass_clust_abs2)})

    prop_correl_clump.update({'M_cl_rel1':np.copy(mass_clust_rel1)})

    prop_correl_clump.update({'M_cl_rel2':np.copy(mass_clust_rel2)})



    #save the dictionary
    name_prop_correl_clump=name+'_prop_correl_struct_'+str(num)+'.save'
    f = open(directory_out+name_prop_correl_clump,'w')
    pickle.dump(prop_correl_clump,f)
    f.close()



    ps=False
    subname=''
    mark=3
    mask = prop_clump['mass'] > 0.

#    P.clf()
#    P.plot(np.log10(prop_clump['mass'][mask]),np.log10(prop_correl_clump['M_cl_tot'][mask]/prop_clump['mass'][mask]),'.',markersize=mark)
#    P.xlabel(r'$log(M) \, (M_ \odot)$')
#    P.ylabel(r'$rap_{tot}$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_mass_rap_tot_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_rap_tot_'+str(num).zfill(5)+'.png')



#    P.clf()
#    P.plot(np.log10(prop_clump['mass'][mask]),np.log10(prop_correl_clump['M_cl_abs1'][mask]/prop_clump['mass'][mask]),'.',markersize=mark)
#    P.xlabel(r'$log(M) \, (M_ \odot)$')
#    P.ylabel(r'$rap_{abs1}$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_mass_rap_abs1_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_rap_abs1_'+str(num).zfill(5)+'.png')


#    P.clf()
#    P.plot(np.log10(prop_clump['mass'][mask]),np.log10(prop_correl_clump['M_cl_abs2'][mask]/prop_clump['mass'][mask]),'.',markersize=mark)
#    P.xlabel(r'$log(M) \, (M_ \odot)$')
#    P.ylabel(r'$rap_{abs2}$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_mass_rap_abs2_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_rap_abs2_'+str(num).zfill(5)+'.png')


#    P.clf()
#    P.plot(np.log10(prop_clump['mass'][mask]),np.log10(prop_correl_clump['M_cl_rel1'][mask]/prop_clump['mass'][mask]),'.',markersize=mark)
#    P.xlabel(r'$log(M) \, (M_ \odot)$')
#    P.ylabel(r'$rap_{rel1}$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_mass_rap_rel1_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_rap_rel1_'+str(num).zfill(5)+'.png')




#    P.clf()
#    P.plot(np.log10(prop_clump['mass'][mask]),np.log10(prop_correl_clump['M_cl_rel2'][mask]/prop_clump['mass'][mask]),'.',markersize=mark)
#    P.xlabel(r'$log(M) \, (M_ \odot)$')
#    P.ylabel(r'$rap_{rel2}$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_mass_rap_rel2_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_rap_rel2_'+str(num).zfill(5)+'.png')


    #cut on the mass to flux ratio, big mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    mass_vlog = np.log10(prop_clump['mass'][mask])
    nbin=50.
    log_min=-3.
    log_max=4.
    width=(log_max-log_min)/nbin
    P.clf()

    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width,color='none',linewidth=1)

    mass_vlog = np.log10(prop_correl_clump['M_cl_tot'][mask])
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='r',linewidth=1)

    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.plot([0.,2.3],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist_correc_mass_mucut_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist_correc_mass_mucut_'+str(num).zfill(5)+'.png')




    #cut on the mass to flux ratio, big mu and select only the non-collapsed objects
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    mass_vlog = np.log10(prop_clump['mass'][mask])
    nbin=50.
    log_min=-3.
    log_max=4.
    width=(log_max-log_min)/nbin
    P.clf()

    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    rho_mean = prop_clump['mass'][mask] / rad[mask]**3  

    mask_n = rho_mean < 1.e5
    if (np.sum(mask_n) > 0):
        hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',linewidth=1)

        mass_vlog = np.log10(prop_correl_clump['M_cl_tot'][mask])
        hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='r',linewidth=1)

    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.plot([0.,2.3],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist_correc_mass_hd_mucut_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist_correc_mass_hd_mucut_'+str(num).zfill(5)+'.png')






    #cut on the mass to flux ratio, big mu and select only the non-collapsed objects
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    mass_vlog = np.log10(prop_clump['mass'][mask])
    nbin=50.
    log_min=-3.
    log_max=4.
    width=(log_max-log_min)/nbin
    P.clf()


    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    rho_mean = prop_clump['mass'][mask] / rad[mask]**3  

    mask_n = rho_mean < 1.e5
    mask_iso = correl_clump['rad_neig'][:,5][mask] > 0.3 # mean distance of 10 closest neighbours must be > 0.1 pc
    mask_n = mask_n * mask_iso
    if (np.sum(mask_n) > 0):
        hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',linewidth=1)


    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_hist_iso_mass_hd_mucut_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_hist_iso_mass_hd_mucut_'+str(num).zfill(5)+'.png')

    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    rho_mean = prop_clump['mass'][mask] / rad[mask]**3  

    mask_n = rho_mean < 1.e5
    mask_iso = correl_clump['rad_neig'][:,5][mask] < 0.3 # mean distance of 10 closest neighbours must be > 0.1 pc
    mask_n = mask_n * mask_iso
    if (np.sum(mask_n) > 0):
        hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='r',linewidth=1)


    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.plot([0.,2.3],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist_iso_mass_hd_mucut_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist_iso_mass_hd_mucut_'+str(num).zfill(5)+'.png')





#    P.clf()
#    P.plot(np.log10(prop_clump['mass'][mask]),np.log10(correl_clump['rad_neig'][:,1][mask]),'.',markersize=mark)
#    P.xlabel(r'$log(M) \, (M_ \odot)$')
#    P.ylabel(r'$R_{2}$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_mass_rad2_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_rad2_'+str(num).zfill(5)+'.png')


    ##2D distogram of mass vs mean neighbour distance
    mass_rad_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['rad_neig'][:,1][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(mass_rad_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$R_{2} \, (pc)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_rad2_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_rad2_'+str(num).zfill(5)+'.png')



    mass_rad_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['rad_neig'][:,5][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(mass_rad_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$R_{2} \, (pc)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_rad10_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_rad10_'+str(num).zfill(5)+'.png')



    mass_rad_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['rad_neig'][:,7][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(mass_rad_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$R_{2} \, (pc)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_rad50_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_rad50_'+str(num).zfill(5)+'.png')



    ##2D histogram of the aspect ratio vs mass
    ar_mass_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['AR_neig'][:,4][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(ar_mass_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'aspect ratio')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_ar5_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_ar5_'+str(num).zfill(5)+'.png')



    ar_mass_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['AR_neig'][:,5][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(ar_mass_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'aspect ratio')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_ar10_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_ar10_'+str(num).zfill(5)+'.png')



    ar_mass_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['AR_neig'][:,7][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(ar_mass_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'aspect ratio')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_ar50_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_ar50_'+str(num).zfill(5)+'.png')



    ar_mass_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['AR2_neig'][:,4][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(ar_mass_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'aspect ratio')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_ar2_5_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_ar2_5_'+str(num).zfill(5)+'.png')



    ar_mass_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['AR2_neig'][:,5][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(ar_mass_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'aspect ratio')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_ar2_10_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_ar2_10_'+str(num).zfill(5)+'.png')



    ar_mass_neig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),correl_clump['AR2_neig'][:,7][mask],bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(ar_mass_neig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'aspect ratio')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_ar2_50_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+subname+'_hist2D_mass_ar2_50_'+str(num).zfill(5)+'.png')










    return name_prop_correl_clump

#######################################################################################
#######################################################################################
#loop over neighboors to find potential and density minimum 
#neighbours are extended to all cells within the neighbouring oct
#path is the path
#num output number
#thres_dens density threshold above which cells are considered
#thres_level level threshold above which cells are considered
#pos_zoom the center of the zoom coordinates
#size_zoom the 3 zoom extension (x, y and z)
#######################################################################################
#######################################################################################
def loop_far_neigh(path,num,thres_dens,thres_level,pos_zoom,size_zoom,gcomp=True):


    ro=pymses.RamsesOutput(path,num)

    amr = ro.amr_source(["rho","phi"], grav_compat=gcomp)

    list_max =[] #list of local maximum
    nmax = 0

    list_max_phi =[] #list of local maximum
    nmax_phi = 0

    for dset in amr.iter_dsets():

#    if(nmax >= 10):
#        break

        ngrids = dset.amr_struct["ngrids"]
        ndim = dset.amr_header["ndim"]
        dset.amr_struct["neighbors"] = -np.ones((ngrids, 2*ndim), dtype='i')
        octree_compute_neighbors(dset)


        dset.get_grid_levels()
        dset.amr_struct["cell_centers"] = dset.get_cell_centers()
        dset.amr_struct["idomain_grid"] = dset.get_idomain_grid()


        #    ndim = dset.amr_struct["ndim"]
        twotondim = 2**ndim

        ngridlevel = dset.amr_struct["ngridlevel"]
        ngridbound = dset.amr_struct["ngridbound"]
        ncpu = dset.amr_struct["ncpu"]
        nboundary = dset.amr_struct["nboundary"]

        son = dset.amr_struct["son_indices"]

        # Initialize array indexes
        iarr = 0

        # Loop over grids
        #    cdef int kcpu
        #    cdef int ind
        ilevel = 0
        igrid = 0
        ngridsremain = 0

        icpu = dset.icpu

    #    cdef numpy.ndarray[int, ndim=2] 
        __nbor_lut = np.array([[0, 2, 4, 6], [0, 1, 4, 5], [0, 1, 2, 3]], 'i')

        while(ilevel < dset.amr_struct["levelmax"]):

            # Skip all CPUs 0..icpu-1 for this level
            for kcpu in range(0, icpu):
                igrid += ngridlevel[kcpu, ilevel]

            # Loop over the CPU grids of this level
            ngridsremain = ngridlevel[icpu, ilevel]
            while(ngridsremain > 0):

                #            if(nmax >= 10):
                #               break


                # Loop over cells
                for ind in range(twotondim):

                    # Only stop at leaf cells
                    if (son[igrid, ind] >= 0):
                        continue

#                if var_array[igrid, ind] > thr:
#                    # Add the current cell to the tree
#                    avl.insert(hi, igrid, ind, iarr)
#                    iarr += 1

#                if (dset["rho"][igrid,ind] > thres_dens and  dset.amr_struct["cell_levels"][igrid] >= thres_level):
                    xg = dset.amr_struct["grid_centers"][igrid,0]
                    yg = dset.amr_struct["grid_centers"][igrid,1]
                    zg = dset.amr_struct["grid_centers"][igrid,2]
                #check that the grid is in the maximum resolution domain
                    if( (xg-pos_zoom[0]) < size_zoom[0]/2. and (yg-pos_zoom[1]) < size_zoom[1]/2. and (zg-pos_zoom[2]) < size_zoom[2]/2. and dset["rho"][igrid,ind] > thres_dens):

                        rho_cell = dset["rho"][igrid,ind]
                        rho_max = rho_cell

                        phi_cell = abs(dset["phi"][igrid,ind])
                        phi_max = phi_cell
                        phi_min = phi_cell
#                    print igrid, ind, dset["rho"][igrid,ind]
#                    print dset.amr_struct["grid_centers"][igrid,0:3]
#                    print    dset.amr_struct["cell_levels"][igrid]

#                    for ineib in range(6):
#                        igridn = dset.amr_struct["neighbors"][igrid,ineib]
#                        print 'ineib ', ineib, 'igrid ', igrid
#                        print dset.amr_struct["grid_centers"][igridn,0:3]
#                        for indn in range(twotondim):
#                             print 'igridn, indn, son ', igridn,indn,son[igridn,indn]



        # Loop over the neighbours
                        indc=ind
                        igridc=igrid
                        for idim in range(ndim):
                            for isgn in range(2):
                                idir = idim * 2 + isgn
                                if (indc >> idim) & 1 == isgn: # The neighbour is in another grid
                                    #                                igridn = nbor[igridc, idir]
                                    igridn = dset.amr_struct["neighbors"][igridc,idir]
                                else: # Same grid
                                    igridn = igridc
                                if igridn != -1: # Neighbour exists
#                                dl = cell_level[igridc] - cell_level[igridn] # dl is -1 or 0
                                    dl = dset.amr_struct["cell_levels"][igridc] - dset.amr_struct["cell_levels"][igridn] # dl is -1 or 0

#                                nsub = 1
                                    nsub = 8 
                                    if dl == 0: # The neighbour is on the same level
                                        indn = indc ^ (1 << idim)

                                        if son[igridn, indn] != -1: # There is one more level available in the neighbouring cell
                                            igridn = son[igridn, indn]
#                                        nsub = 1 << (ndim - 1)
                                            dl = 1
                                    else: # The neighbour is coarser
                                        indn = 0
                                        for idim2 in range(ndim):
                                            if idim2 != idim:
                                                if dset.amr_struct["grid_centers"][igridc, idim2] > dset.amr_struct["grid_centers"][igridn, idim2]:
                                                    indn |= 1 << idim2
                                            else:
                                                indn |= (1 - isgn) << idim2

#                                print 'neigh grid ',dset.amr_struct["grid_centers"][igridn,0:3]
#                                print 'nsub ',nsub
                                    for isub in range(nsub): # Iterate through neighbouring sub-cells
#                                    if dl > 0: # Finer level, more than one neighbour in this direction
#                                        indn = __nbor_lut[idim, isub] ^ ((1 - isgn) << idim)

#                                    print 'isub ',isub , ' ,indn ',  indn , ',igridn ',igridn
#                                    print 'neigh cell ',dset.amr_struct["cell_centers"][igridn,indn]-dset.amr_struct["cell_centers"][igrid,ind]
#                                    rho_neigh = dset["rho"][igridn,indn]
                                        rho_neigh = dset["rho"][igridn,isub]
                                        rho_max = max(rho_max,rho_neigh)

                                        phi_neigh = abs(dset["phi"][igridn,isub])
                                        phi_max = max(phi_max,phi_neigh)
                                        phi_min = min(phi_min,phi_neigh)

#                    pdb.set_trace() # debug mode     
                        if(rho_max == rho_cell): #this is a local maximum
                            list_max.append((igridn,indn,icpu,rho_cell,phi_cell,dset.amr_struct["cell_centers"][igrid,ind]))
                            nmax = nmax + 1
#                        print 'local dens max, nmax ',nmax, 'igridn ',igridn, 'indn ',indn , 'icpu ',icpu
#                        print rho_max, dset.amr_struct["cell_centers"][igrid,ind]
#                        print 'nmax_phi ',nmax_phi, 'phi_max ',phi_max, 'phi_cell ',phi_cell



                        if(phi_max == phi_cell): #this is a local maximum
                            list_max_phi.append((igridn,indn,icpu,rho_cell,phi_cell,dset.amr_struct["cell_centers"][igrid,ind],phi_min))
                            nmax_phi = nmax_phi + 1
                            print 'local phi max, nmax_phi ',nmax_phi, 'igridn ',igridn, 'indn ',indn , 'icpu ',icpu
                            print dset.amr_struct["cell_centers"][igrid,ind]

                            
            # Go to next grid at this CPU and level
                ngridsremain -= 1
                igrid += 1

        # Skip all CPUs icpu+1, ... ncpu for this level
            for kcpu in range(icpu+1, ncpu):
                igrid += ngridlevel[kcpu, ilevel]

        # Skip all boundary grids for this level
            for kcpu in range(nboundary):
                igrid += ngridbound[kcpu, ilevel]

            ilevel += 1


    name='local_maximum_'+str(num)+'.save'
    f = open(name,'w')
    pickle.dump(nmax,f)
    pickle.dump(list_max,f)
    pickle.dump(nmax_phi,f)
    pickle.dump(list_max_phi,f)
    f.close()



    return name

#######################################################################################
#######################################################################################
#loop over oct to compute the tidal forces
#path is the path
#num output number
#thres_dens density threshold above which cells are considered
#thres_level level threshold above which cells are considered
#pos_zoom the center of the zoom coordinates
#size_zoom the 3 zoom extension (x, y and z)
#######################################################################################
#######################################################################################
def loop_oct_tidal(path,num,path_out=None,gcomp=True):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in



    ro=pymses.RamsesOutput(path,num)

    amr = ro.amr_source(["rho","phi","g"], grav_compat=gcomp)
#    amr = ro.amr_source(["rho","phi","g"], grav_compat=False)

    list_max =[] #list of local maximum
    nmax = 0

    list_max_phi =[] #list of local maximum
    nmax_phi = 0


    mat_tidal=np.zeros((3,3))



    w_v=[]
    rho_v=[]
    x_v=[]
    y_v=[]
    z_v=[]
    dx_v=[]


    w2_v=[]
    rho2_v=[]


    ncell=0
    ncell_neg=0
    ncell_tid2_neg=0
    ncell_tid3_neg=0
    ncell2_neg=0
    ncell2_tid2_neg=0
    ncell2_tid3_neg=0



    for dset in amr.iter_dsets():

        ngrids = dset.amr_struct["ngrids"]
        ndim = dset.amr_header["ndim"]
        dset.amr_struct["neighbors"] = -np.ones((ngrids, 2*ndim), dtype='i')
        octree_compute_neighbors(dset)


        dset.get_grid_levels()
        dset.amr_struct["cell_centers"] = dset.get_cell_centers()
        dset.amr_struct["idomain_grid"] = dset.get_idomain_grid()


        #    ndim = dset.amr_struct["ndim"]
        twotondim = 2**ndim

        ngridlevel = dset.amr_struct["ngridlevel"]
        ngridbound = dset.amr_struct["ngridbound"]
        ncpu = dset.amr_struct["ncpu"]
        nboundary = dset.amr_struct["nboundary"]

        son = dset.amr_struct["son_indices"]

        # Initialize array indexes
        iarr = 0

        # Loop over grids
        #    cdef int kcpu
        #    cdef int ind
        ilevel = 0
        igrid = 0
        ngridsremain = 0

        icpu = dset.icpu

    #    cdef numpy.ndarray[int, ndim=2] 
#        __nbor_lut = np.array([[0, 2, 4, 6], [0, 1, 4, 5], [0, 1, 2, 3]], 'i')


        print 'ncell ',ncell , 'ncell_neg ',ncell_neg, 'ncell_tid2_neg ', ncell_tid2_neg , 'ncell_tid3_neg ', ncell_tid3_neg 
        print 'ncell ',ncell , 'ncell2_neg ',ncell2_neg, 'ncell2_tid2_neg ', ncell2_tid2_neg , 'ncell2_tid3_neg ', ncell2_tid3_neg 


        while(ilevel < dset.amr_struct["levelmax"]):


            print 'ilevel',ilevel

            dx = 0.5**ilevel

            # Skip all CPUs 0..icpu-1 for this level
            for kcpu in range(0, icpu):
                igrid += ngridlevel[kcpu, ilevel]

            # Loop over the CPU grids of this level
            ngridsremain = ngridlevel[icpu, ilevel]
            while(ngridsremain > 0):



                #select the oct which have only leaf cells
                find_leaf=True
                # Loop over cells to check whether all of them are leaf
                for ind in range(twotondim):
                    # Only stop at leaf cells
                    if (son[igrid, ind] >= 0):
                        find_leaf=False


                #select the oct which have at least one leaf cell
                find_leaf=False
                # Loop over cells to check whether all of them are leaf
                for ind in range(twotondim):
                    # Only stop at leaf cells
                    if (son[igrid, ind] == -1):
                       find_leaf=True



                if(find_leaf):
#                    print 'igrid',igrid

                    #position of the grid
                     xg = dset.amr_struct["grid_centers"][igrid,0]
                     yg = dset.amr_struct["grid_centers"][igrid,1]
                     zg = dset.amr_struct["grid_centers"][igrid,2]

                     
                    #average over the oct and calculate the eigen values
                     g_xm = 0.25*(dset["g"][igrid,0]+dset["g"][igrid,2]+dset["g"][igrid,4]+dset["g"][igrid,6])
                     g_xp = 0.25*(dset["g"][igrid,1]+dset["g"][igrid,3]+dset["g"][igrid,5]+dset["g"][igrid,7])

                     g_ym = 0.25*(dset["g"][igrid,0]+dset["g"][igrid,1]+dset["g"][igrid,4]+dset["g"][igrid,5])
                     g_yp = 0.25*(dset["g"][igrid,2]+dset["g"][igrid,3]+dset["g"][igrid,6]+dset["g"][igrid,7])

                     g_zm = 0.25*(dset["g"][igrid,0]+dset["g"][igrid,1]+dset["g"][igrid,2]+dset["g"][igrid,3])
                     g_zp = 0.25*(dset["g"][igrid,4]+dset["g"][igrid,5]+dset["g"][igrid,6]+dset["g"][igrid,7])

                     mat_tidal[:,0] = (g_xp - g_xm) / dx
                     mat_tidal[:,1] = (g_yp - g_ym) / dx
                     mat_tidal[:,2] = (g_zp - g_zm) / dx

                     w,v=eig(mat_tidal)

                     w = np.sort(w)

                     w_v.append(w)

                     rho = np.sum(dset["rho"][igrid,:])/8.
                     rho_v.append(rho)  

                     dx_v.append(dx)  

                     x_v.append(xg)  
                     y_v.append(yg)  
                     z_v.append(zg)  

                     ncell = ncell + 1. 

                     if(np.sum(w) <= 0.):
                         ncell_neg = ncell_neg + 1


                     if(w[2]<=0.):
                         ncell_tid3_neg = ncell_tid3_neg + 1


                     if(w[1]<=0.):
                         ncell_tid2_neg = ncell_tid2_neg + 1


                    #average over the oct and calculate the eigen values
                     g_xm = dset["g"][igrid,0]
                     g_xp = dset["g"][igrid,1]

                     g_ym = dset["g"][igrid,0]
                     g_yp = dset["g"][igrid,2]

                     g_zm = dset["g"][igrid,0]
                     g_zp = dset["g"][igrid,4]

                     mat_tidal[:,0] = (g_xp - g_xm) / dx
                     mat_tidal[:,1] = (g_yp - g_ym) / dx
                     mat_tidal[:,2] = (g_zp - g_zm) / dx

                     w,v=eig(mat_tidal)

                     w = np.sort(w)

                     w2_v.append(w)

                     rho = dset["rho"][igrid,0]
                     rho2_v.append(rho)  



#                     print 'trace ', mat_tidal[0,0]+mat_tidal[1,1]+ mat_tidal[2,2]
#                     print np.sum(w)

                     if(np.sum(w) <= 0.):
                         ncell2_neg = ncell2_neg + 1


                     if(w[2]<=0.):
                         ncell2_tid3_neg = ncell2_tid3_neg + 1


                     if(w[1]<=0.):
                         ncell2_tid2_neg = ncell2_tid2_neg + 1


#                     pdb.set_trace() # debug mode     



#                    for ind in range(twotondim):
#                        print dset["g"][igrid,ind]
#                mat_tidal[0,0] = dset["g"][igrid,ind][] - dset["g"][igrid,ind][]
                



#                if (dset["rho"][igrid,ind] > thres_dens and  dset.amr_struct["cell_levels"][igrid] >= thres_level):
#                    xg = dset.amr_struct["grid_centers"][igrid,0]
#                    yg = dset.amr_struct["grid_centers"][igrid,1]
#                    zg = dset.amr_struct["grid_centers"][igrid,2]


#                    print 'igrid: ',igrid,'ind: ',ind, dset["rho"][igrid,ind],dset["g"][igrid,ind]


#                     pos = dset.amr_struct["grid_centers"][igrid,:]

#                     pos_v.append(pos)

                            
            # Go to next grid at this CPU and level
                ngridsremain -= 1
                igrid += 1

        # Skip all CPUs icpu+1, ... ncpu for this level
            for kcpu in range(icpu+1, ncpu):
                igrid += ngridlevel[kcpu, ilevel]

        # Skip all boundary grids for this level
            for kcpu in range(nboundary):
                igrid += ngridbound[kcpu, ilevel]

            ilevel += 1



    rho_v  = np.asarray(rho_v)
    w_v = np.asarray(w_v)

    rho2_v  = np.asarray(rho_v)
    w2_v = np.asarray(w_v)

    dx_v = np.asarray(dx_v)

    x_v = np.asarray(x_v)
    y_v = np.asarray(y_v)
    z_v = np.asarray(z_v)

#    w_v    = np.asarray(w_v)


#    mask = [np.sum(wtot_v) < 0.]

#    print 'number of cells) ', len(wtot_v)
#    print 'number of cells with positive trace(tidal) ', np.sum(mask)



#    P.clf()
#    P.plot(np.log10(rho_v),np.log10(abs(wtot_v))*np.sign(wtot_v),'.')
#    P.savefig(directory_out+'/rho_tidal'+'_'+str(num)+'.ps')


    name=directory_out+'/tidal_'+str(num)+'.save'
#    f = open(name,'w')
#    pickle.dump(rho_v,f)
#    pickle.dump(wtot_v,f)
#    pickle.dump(w_v,f)
#    pickle.dump(nmax_phi,f)
#    pickle.dump(list_max_phi,f)
#    f.close()

    np.savez(name,rho_v=rho_v,w_v=w_v,x_v=x_v,y_v=y_v,z_v=z_v,w2_v=w2_v,rho2_v=rho2_v,dx_v=dx_v)


    return 

#######################################################################################
#######################################################################################
## display the tidal forces using the results of loop_oct_tidal
#######################################################################################
#######################################################################################
def display_tidal(path_in,num,path_out=None):


    if(path_out is None):
        directory_out=path_in
    else:
        directory_out=path_out

    name=path_in+'tidal_'+str(num)+'.save.npz'
    npzfile = np.load(name)

    rho_v = npzfile['rho_v']
    w_v   = npzfile['w_v']
    x_v = npzfile['x_v']
    y_v = npzfile['y_v']
    z_v = npzfile['z_v']

    dx_v = npzfile['dx_v']


    mask1 =  np.imag(w_v[:,1]) == 0. 
    mask2 = w_v[:,0] < 0.
    mask3 = w_v[:,1] != 0.
    mask4 = w_v[:,2] != 0.

    mask=mask1*mask2*mask3*mask4

    rho_v = rho_v[mask]
    w1_v   = np.real(w_v[:,0][mask])
    w2_v   = np.real(w_v[:,1][mask])
    w3_v   = np.real(w_v[:,2][mask])

    dx_v = dx_v[mask]

    mass_v = ((dx_v)**3)*rho_v

    wtot_v = w1_v+w2_v+w3_v



    logwt_v = np.log10(abs(wtot_v))*np.sign(wtot_v)

    logw_w1_v = np.log10(abs(wtot_v/w1_v))*np.sign(wtot_v)

    logw1_v = np.log10(abs(w1_v))*np.sign(w1_v)
    logw2_v = np.log10(abs(w2_v))*np.sign(w2_v)
    logw3_v = np.log10(abs(w3_v))*np.sign(w3_v)


    logw13_v = np.log10(abs(w3_v/w1_v))
    
    logw12_v = np.log10(abs(w2_v/w1_v))


##1D histogram
    log_rho_min = 3.
    log_rho_max = 13.

    nbin=50.

    width = (log_rho_max - log_rho_min) / nbin

    P.clf() 


    hist_mass , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v)

    P.bar(hist_edges[:-1],hist_mass,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$\log(M) \, (code units)$')


    P.savefig(directory_out+'hist_mass'+str(num)+'.ps')
    P.savefig(directory_out+'hist_mass'+str(num)+'.jpeg')


    P.clf() 

    hist_Tidal1 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*(w1_v))

    hist_Tidal1 = abs(hist_Tidal1) / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal1,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$\log(\lambda1) \, (code units)$')


    P.savefig(directory_out+'hist_T1'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T1'+'_'+str(num)+'.jpeg')


    P.clf() 

    hist_Tidal2 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*w2_v)

    hist_Tidal2 = abs(hist_Tidal2) / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal2,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$\log(\lambda1) \, (code units)$')


    P.savefig(directory_out+'hist_T2'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T2'+'_'+str(num)+'.jpeg')


    P.clf() 

    hist_Tidal3 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*w3_v)

    hist_Tidal3 = abs(hist_Tidal3) / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal3,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$\log(\lambda1) \, (code units)$')


    P.savefig(directory_out+'hist_T3'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T3'+'_'+str(num)+'.jpeg')



    P.clf() 

    P.bar(hist_edges[:-1],hist_Tidal2/hist_Tidal1,width)

#P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$\lambda2 / \lambda1) $')


    P.savefig(directory_out+'hist_T12'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T12'+'_'+str(num)+'.jpeg')




    P.clf() 

    mask_pos = w2_v > 0.

    hist_Tidal2_pos , hist_edges = P.histogram(np.log10(rho_v[mask_pos]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_pos]*(w2_v[mask_pos]))


    mask_neg = w2_v < 0.
    hist_Tidal2_neg , hist_edges = P.histogram(np.log10(rho_v[mask_neg]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_neg]*abs(w2_v[mask_neg]))

    hist_Tidal2_rap = hist_Tidal2_pos / hist_Tidal2_neg

    P.bar(hist_edges[:-1],hist_Tidal2_rap,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$r_2$')


    P.savefig(directory_out+'hist_T22'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T22'+'_'+str(num)+'.jpeg')




    P.clf() 

    P.bar(hist_edges[:-1],hist_Tidal3/hist_Tidal1,width)

#P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$\lambda3 / \lambda1) $')


    P.savefig(directory_out+'hist_T13'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T13'+'_'+str(num)+'.jpeg')



    P.clf() 

    w12_v = w1_v + w2_v
    hist_Tidal1p2 , hist_edges = P.histogram(np.log10(rho_v),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v*abs(w12_v))

    hist_Tidal1p2 = hist_Tidal1p2 / hist_mass

    P.bar(hist_edges[:-1],hist_Tidal3/hist_Tidal1p2,width)

#P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$\lambda3 / (\lambda1+\lambda2)) $')


    P.savefig(directory_out+'hist_T1p23'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T1p23'+'_'+str(num)+'.jpeg')




    P.clf() 

    mask_pos = w3_v > 0.

    hist_Tidal3_pos , hist_edges = P.histogram(np.log10(rho_v[mask_pos]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_pos]*(w3_v[mask_pos]))


    mask_neg = w3_v < 0.
    hist_Tidal3_neg , hist_edges = P.histogram(np.log10(rho_v[mask_neg]),bins=nbin ,range=(log_rho_min,log_rho_max),weights=mass_v[mask_neg]*abs(w3_v[mask_neg]))

    hist_Tidal3_rap = hist_Tidal3_pos / hist_Tidal3_neg

    P.bar(hist_edges[:-1],hist_Tidal3_rap,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\rho (cm^{-3})$')         
    P.ylabel(r'$r_3$')


    P.savefig(directory_out+'hist_T33'+'_'+str(num)+'.ps')
    P.savefig(directory_out+'hist_T33'+'_'+str(num)+'.jpeg')


    return
#######################################################################################
#######################################################################################
## do spatial map of the tidal forces using the results of loop_oct_tidal
#######################################################################################
#######################################################################################
def display_tidal_map(path_in,num,path_out=None,sinks=True,no_cvs=True):


    if(path_out is None):
        directory_out=path_in
    else:
        directory_out=path_out

    ro=pymses.RamsesOutput(path_in,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    print 'lbox ',lbox


    #path out is used because it is expected that the result of loop_oct_tidal has been written there
    name=path_out+'tidal_'+str(num)+'.save.npz'
    npzfile = np.load(name)

    rho_v = npzfile['rho_v']
    w_v   = npzfile['w_v']
    x_v = npzfile['x_v']
    y_v = npzfile['y_v']
    z_v = npzfile['z_v']

    dx_v = npzfile['dx_v']

    mask1 =  np.imag(w_v[:,1]) == 0. 
    mask2 = w_v[:,0] < 0.
    mask3 = w_v[:,1] != 0.
    mask4 = w_v[:,2] != 0.

    mask=mask1*mask2*mask3*mask4


    rho_v = rho_v[mask]
    w1_v   = np.real(w_v[:,0][mask])
    w2_v   = np.real(w_v[:,1][mask])
    w3_v   = np.real(w_v[:,2][mask])


    dx_v = dx_v[mask]
    x_v = x_v[mask]
    y_v = y_v[mask]
    z_v = z_v[mask]




#    w1_v   = np.real(w_v[:,0])
#    w2_v   = np.real(w_v[:,1])
#    w3_v   = np.real(w_v[:,2])
#    mask1 = w1_v == 0.
#    w1_v[mask1] = -abs(w3_v[mask1]) * 10.
#    mask1 = w1_v == 0.
#    w1_v[mask1]  = 1.

    rapw_v = w3_v / w1_v


    if(sinks):
    #read sink file

        pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

        sinks = read_sink_cvs(num,path_in,no_cvs=no_cvs)


        if(len(sinks) !=0):
            mass_v = sinks[:,1] * scale_mass / Ms / lbox**3
            xs_v = sinks[:,2]
            ys_v = sinks[:,3]
            zs_v = sinks[:,4]

#position of the most massive sink
            ind_max = np.argmax(sinks[:,1])
            center_max=[xs_v[ind_max]/lbox,ys_v[ind_max]/lbox,zs_v[ind_max]/lbox]


        else:
            center_max=[0.5,0.5,0.5]


    else:
            center_max=[0.5,0.5,0.5]
    

    #do a zoom
    nzoom=4
    eps=0.5**nzoom
    map_coldens , map_w13, xedges, yedges = make_hierarch_map(x_v,y_v,z_v,dx_v,rho_v,rapw_v,eps,center=center_max)

    map_coldens = map_coldens.T 
    map_w13 = map_w13.T

    xedges = xedges * lbox
    yedges = yedges * lbox

    mask = map_coldens != 0.
    map_w13[mask] = map_w13[mask] / map_coldens[mask]

    map_coldens = map_coldens* unit_col

    log_map_coldens = map_coldens
    log_map_coldens[mask] = np.log10(log_map_coldens[mask])


    P.clf()
    xmin = (center_max[0]-eps)*lbox
    xmax = (center_max[0]+eps)*lbox
    ymin = (center_max[1]-eps)*lbox
    ymax = (center_max[1]+eps)*lbox

    im = P.imshow(log_map_coldens,origin='lower',extent=[xmin,xmax,ymin,ymax] )

    cbar = P.colorbar(im)                                                                                                          

    cbar.set_label(r'$log(N) (cm^{-2})$')

    P.xlabel(r'$x$ (pc)')     
    P.ylabel(r'$y$ (pc)')

    P.savefig(directory_out+'/map_rho_tidal'+'_'+str(num)+'_'+str(nzoom)+'.ps')
    P.savefig(directory_out+'/map_rho_tidal'+'_'+str(num)+'_'+str(nzoom)+'.png')



    P.clf()
    im = P.imshow(map_w13,origin='lower',extent=[xmin,xmax,ymin,ymax] )

    cbar = P.colorbar(im)                                                                                                          

    cbar.set_label(r'$\lambda_3 / \lambda_1$')

    P.xlabel(r'$x$ (pc)')     
    P.ylabel(r'$y$ (pc)')

    P.savefig(directory_out+'/map_tidal13'+'_'+str(num)+'_'+str(nzoom)+'.ps')
    P.savefig(directory_out+'/map_tidal13'+'_'+str(num)+'_'+str(nzoom)+'.png')



    #do a deeper zoom
    nzoom=6
    eps=0.5**nzoom
    map_coldens , map_w13, xedges, yedges = make_hierarch_map(x_v,y_v,z_v,dx_v,rho_v,rapw_v,eps,center=center_max)

    map_coldens = map_coldens.T * unit_col
    map_w13 = map_w13.T

    xedges = xedges * lbox
    yedges = yedges * lbox

    mask = map_coldens != 0.
    map_w13[mask] = map_w13[mask] / map_coldens[mask]

    log_map_coldens = map_coldens
    log_map_coldens[mask] = np.log10(log_map_coldens[mask])


    P.clf()

    xmin = (center_max[0]-eps)*lbox
    xmax = (center_max[0]+eps)*lbox
    ymin = (center_max[1]-eps)*lbox
    ymax = (center_max[1]+eps)*lbox
    im = P.imshow(log_map_coldens,origin='lower',extent=[xmin,xmax,ymin,ymax] )

    cbar = P.colorbar(im)                                                                                                          

    cbar.set_label(r'$log(N) (cm^{-2})$')

    P.xlabel(r'$x$ (pc)')     
    P.ylabel(r'$y$ (pc)')

    P.savefig(directory_out+'/map_rho_tidal'+'_'+str(num)+'_'+str(nzoom)+'.ps')
    P.savefig(directory_out+'/map_rho_tidal'+'_'+str(num)+'_'+str(nzoom)+'.png')



    P.clf()
    im = P.imshow(map_w13,origin='lower', extent=[xmin,xmax,ymin,ymax] )

    cbar = P.colorbar(im)                                                                                                          

    cbar.set_label(r'$\lambda_3 / \lambda_1$')

    P.xlabel(r'$x$ (pc)')     
    P.ylabel(r'$y$ (pc)')

    P.savefig(directory_out+'/map_tidal13'+'_'+str(num)+'_'+str(nzoom)+'.ps')
    P.savefig(directory_out+'/map_tidal13'+'_'+str(num)+'_'+str(nzoom)+'.png')


    return

##########################################################################
##########################################################################
#build a column density map and column density weighted map of field
# by summing hierachically on cells
#x_v,y_v,rho_v,field_v are list of fiels from which the map is built
# x_v, y_v and z_v between 0 an 1
#nbin number of bins, eps is the spatial extansion 
# center is the center of the map (between 0 an 1)
##########################################################################
##########################################################################
def make_hierarch_map(x_v,y_v,z_v,dx_v,rho_v,field_v,eps,center=[0.5,0.5,0.5],nref = 9):

    mask1 = abs(x_v - center[0]) <= eps
    mask2 = abs(y_v - center[1]) <= eps
    mask3 = abs(z_v - center[2]) <= eps

    mask = mask1*mask2*mask3

    rho_v = rho_v[mask]
    field_v = field_v[mask]
    dx_v = dx_v[mask]

    x_v = x_v[mask]
    y_v = y_v[mask]


    xmin=center[0]-eps
    xmax=center[0]+eps
    ymin=center[1]-eps
    ymax=center[1]+eps

    #prepare maps by summing recursively starting from low resolution
    map_coldens_tot = np.zeros((1,1)) + 1.
    map_field_tot = np.zeros((1,1)) + 1.

    for i in range(1,nref,1):
        map_coldens_tot = np.repeat(map_coldens_tot,2,axis=0)
        map_coldens_tot = np.repeat(map_coldens_tot,2,axis=1)
        map_field_tot = np.repeat(map_field_tot,2,axis=0)
        map_field_tot = np.repeat(map_field_tot,2,axis=1)


        nbin=2**i
        dx_ref = 2.*eps / nbin
        mask = abs(dx_v - dx_ref) < dx_ref / 10.
        if(np.sum(mask) > 0): 
            map_coldens , xedges, yedges = np.histogram2d(x_v[mask],y_v[mask],bins=nbin,weights=rho_v[mask]*dx_v[mask],range=[[xmin,xmax],[ymin,ymax]])
            map_field , xedges, yedges = np.histogram2d(x_v[mask],y_v[mask],bins=nbin,weights=rho_v[mask]*dx_v[mask]*field_v[mask],range=[[xmin,xmax],[ymin,ymax]])

            map_coldens_tot = map_coldens_tot + map_coldens
            map_field_tot = map_field_tot + map_field


#            P.clf()
#            im = P.imshow(np.log10(map_coldens_tot),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]] )
#            cbar = P.colorbar(im)                                                  
#            cbar.set_label(r'$N$')
#            P.savefig('map_N'+'_'+str(i)+'.png')


    map_coldens_tot = np.repeat(map_coldens_tot,2,axis=0)
    map_coldens_tot = np.repeat(map_coldens_tot,2,axis=1)
    map_field_tot = np.repeat(map_field_tot,2,axis=0)
    map_field_tot = np.repeat(map_field_tot,2,axis=1)

    nbin=2**nref
    dx_ref = 2.*eps / nbin
    mask = dx_v  <= dx_ref * 1.1

    if(np.sum(mask) > 0):
        map_coldens , xedges, yedges = np.histogram2d(x_v[mask],y_v[mask],bins=nbin,weights=rho_v[mask]*dx_v[mask],range=[[xmin,xmax],[ymin,ymax]])
        map_field , xedges, yedges = np.histogram2d(x_v[mask],y_v[mask],bins=nbin,weights=rho_v[mask]*dx_v[mask]*field_v[mask],range=[[xmin,xmax],[ymin,ymax]])

        map_coldens_tot = map_coldens_tot + map_coldens
        map_field_tot = map_field_tot + map_field


#        P.clf()
#        im = P.imshow(np.log10(map_coldens_tot),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]] )
#        cbar = P.colorbar(im)                                                  
#        cbar.set_label(r'$N$')
#        P.savefig('map_N'+'_'+str(nref)+'.png')


#    pdb.set_trace() # debug mode     

    return map_coldens_tot, map_field_tot, xedges, yedges

##########################################################################
##########################################################################
#check whether a given cell is contained in a list
#the gravitational potential is used for this and the list
#has been sorted by potential
##########################################################################
##########################################################################
def test_neigh(ic,liste_cell,neigh_cell_am,cells):

    nc = len(liste_cell)
    imin=0
    imax=nc-1

    phi_test = neigh_cell_am["phi"][ic]


    if( cells["phi"][liste_cell[imax]] < phi_test):
        test=False
        return test

    #do an iterative procedure to check whether the potential matches a value in the list
    while (imax-imin > 1):
        imed = (imin+imax)/2

        if(cells["phi"][liste_cell[imed]] <= phi_test):
            imin = imed
        else: 
            imax = imed

    if(cells["phi"][liste_cell[imin]] == phi_test or cells["phi"][liste_cell[imax]] == phi_test):
        test=True
        return test

    test=False
    return test


########################################################################################
########################################################################################
#merge the local extremum in they are closer than merge_length
#path the path
#num the output number
#name the name of the local potential minimum position
#fait_image make image of the cores
#freq_dump frequency of dump onto the disk
########################################################################################
########################################################################################
#def merge_local_max(path,num,name,merge_length=1.e-5):
def merge_local_max(name,merge_length=1.e-5):

#    ro=pymses.RamsesOutput(path,num)
#    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)
#    lbox_pc = lbox*3.08e18
#    amr = ro.amr_source(["rho","vel","phi"], grav_compat=True)


    print 'WORK IN PROGRESS...'
    return

    f = open(name,'r')
    nmax = pickle.load(f)
    list_max = pickle.load(f)
    nmax_phi = pickle.load(f)
    list_max_phi = pickle.load(f)
    f.close()


    #start by sorting the local extremum by x
    def getkey(item):
        return item[5][0]

    list_max_phi = sorted(list_max_phi,reverse=True,key=getkey)


    nmax = len(list_max_phi)

    merge_lengthon2 = merge_length**2


    #loop over the extremum and define groups
    for imax in range(nmax):
        ip = imax +1
        testp=True
        while(ip < nmax and testp):
            if(radon2 < merge_lengthon2):
                ip = ip + 1 




    return name_merged_list

########################################################################################
########################################################################################
#visualise the local maximum found in loop_far_neigh in column densities map
#path the path
#num the output number
#name the name of the local potential minimum tuple which contains the positions
#pos_zoom the center of the zoom coordinates
#size_zoom the 3 zoom extension (x, y and z)
#n_image the number of images along each axis
#xx and yy if specified should be arrays specifying the coordinates of the points
# they are expected to be expressed in lbox (and not in 1)
##to be drawn (expected values between 0 and 1)
##rest_struct restore an existing structure file as produced in extract_struct
########################################################################################
########################################################################################
def visualise_loc_max(path,num,name,pos_zoom,size_zoom,n_image,xx=None,yy=None,path_out=None,rest_struct=False,gcomp=True,num_thr=100):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path


    ro=pymses.RamsesOutput(path,num)

    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)
    lbox_cm = lbox*3.08e18


    amr = ro.amr_source(["rho","vel","phi"], grav_compat=gcomp)


    if(rest_struct):

        name_prop_clump='prop_struct_'+str(num)+'.save'

        f = open(path+name_prop_clump,'r')
        prop_clump = pickle.load(f)
        thres_v = pickle.load(f)
        f.close()


        xpos_v = prop_clump['pos_n_max'][:,0] / lbox
        ypos_v = prop_clump['pos_n_max'][:,1] / lbox

        mu_v = prop_clump['mu']


    else:

        if(xx is None or yy is None):

            f = open(name,'r')
            nmax = pickle.load(f)
            list_max = pickle.load(f)
            nmax_phi = pickle.load(f)
            list_max_phi = pickle.load(f)
            f.close()

            def getkey(item):
                return item[3]

            #[i[3] for i in list_max]
            #list_max = sorted(list_max,reverse=True,key=getkey)
            list_max_phi = sorted(list_max_phi,reverse=True,key=getkey)

            #extract the x and y coordinates of the local extremum
            xpos_v =  np.asarray([i[5][0] for i in list_max_phi])
            ypos_v =  np.asarray([i[5][1] for i in list_max_phi])



        else:
            xpos_v = xx / lbox
            ypos_v = yy / lbox


    radius_x = size_zoom[0]/n_image/2.
    radius_y = size_zoom[1]/n_image/2.

    radius_z = size_zoom[2]/2. #we do not divide by n_image here because we integrate ovver the box

    center = np.zeros(3)

    for ii in range(n_image):
            for jj in range(n_image):


                center[0] = pos_zoom[0] + (ii - (n_image-1.)/2.)/(n_image/2.)*size_zoom[0]/2.
                center[1] = pos_zoom[1] + (jj - (n_image-1.)/2.)/(n_image/2.)*size_zoom[1]/2.
                center[2] = pos_zoom[2]
            

                apos = np.where( (xpos_v >= -radius_x+center[0]) & (xpos_v <= radius_x+center[0]) & (ypos_v >= -radius_y+center[1]) & (ypos_v <= radius_y+center[1]) == True)


                if(rest_struct):
                                    apos_muh = np.where( (xpos_v >= -radius_x+center[0]) & (xpos_v <= radius_x+center[0]) & (ypos_v >= -radius_y+center[1]) & (ypos_v <= radius_y+center[1]) & (mu_v >= 1.) == True)

                                    apos_mul = np.where( (xpos_v >= -radius_x+center[0]) & (xpos_v <= radius_x+center[0]) & (ypos_v >= -radius_y+center[1]) & (ypos_v <= radius_y+center[1]) & (mu_v <= 1.) == True)


                #draw only an image if it contains more tha 100 local extremmum
                if(len(apos[0]) <= num_thr):
                    continue



                cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius_x,2.*radius_y],distance=radius_z,far_cut_depth=radius_z,up_vector='y',map_max_size=512)

                rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
                rt = raytracing.RayTracer(amr,ro.info,rho_op)
                datamap = rt.process(cam_z,surf_qty=True)
                map_col = np.log10(datamap.map.T*lbox_cm)


                P.clf() 
                im = P.imshow(map_col,extent=[(-radius_x+center[0])*lbox,(radius_x+center[0])*lbox,(-radius_y+center[1])*lbox,(radius_y+center[1])*lbox],origin='lower')   

                P.xlabel('$x$ (pc)')     
                P.ylabel('$y$ (pc)')
                cbar = P.colorbar(im)                                                                                                          

                cbar.set_label(r'$log(N) \, cm^{-2}$')                                          


#                pdb.set_trace() # debug mode     



                if(rest_struct):
                    #distinguish between super and subcritical

                    P.scatter((xpos_v[apos_mul])*lbox,(ypos_v[apos_mul])*lbox,c='m',marker='d',s=10)

                    P.scatter((xpos_v[apos_muh])*lbox,(ypos_v[apos_muh])*lbox,c='r',marker='x',s=10)


                else:
                    P.scatter((xpos_v[apos])*lbox,(ypos_v[apos])*lbox,c='r',marker='x',s=10)

                P.savefig(directory_out+'visu_max_z'+'_'+str(n_image)+'_'+str(ii)+'_'+str(jj)+'_'+str(num)+'.jpeg')


#    pdb.set_trace() # debug mode     

    return
########################################################################################
########################################################################################
##now do the same but using sink particles
def visualise_sink(path,num,name,pos_zoom,size_zoom,n_image,path_out=None,conv_sink=1):


    ro=pymses.RamsesOutput(path,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)


    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

    sinks = read_sink_cvs(num,path,no_cvs=no_cvs)


    if(len(sinks) !=0):
        mass_v = sinks[:,1] * scale_mass / Ms / lbox**3
        x_v = sinks[:,2]
        y_v = sinks[:,3]
        z_v = sinks[:,4]

        if(conv_sink == 2):
            x_v = sinks[:,3]
            y_v = sinks[:,4]
            z_v = sinks[:,5]


    visualise_loc_max(path,num,name,pos_zoom,size_zoom,n_image,xx=x_v,yy=y_v,path_out=path_out,rest_struct=False,gcomp=True,num_thr=0)


    return
########################################################################################
########################################################################################
#extract core starting from the minimum potential found in loop_far_neigh
#path the path
#num the output number
#name the name of the local potential minimum position
#fait_image make image of the cores
#freq_dump frequency of dump onto the disk
#ntol number of times the potential may increase before core stops
########################################################################################
########################################################################################
def extract_core(path,num,name,fait_image=False,freq_dump=200,ntol=1,gcomp=True):


    ro=pymses.RamsesOutput(path,num)

    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)
    lbox_pc = lbox*3.08e18


    amr = ro.amr_source(["rho","vel","phi"], grav_compat=gcomp)


    f = open(name,'r')
    nmax = pickle.load(f)
    list_max = pickle.load(f)
    nmax_phi = pickle.load(f)
    list_max_phi = pickle.load(f)
    f.close()




    def getkey(item):
        return item[3]

#[i[3] for i in list_max]
#list_max = sorted(list_max,reverse=True,key=getkey)

    list_max_phi = sorted(list_max_phi,reverse=True,key=getkey)



##loop over local minimum to identify the core and do some images
#n_im = 500 #nmax_phi / 20
    n_im =  nmax_phi #/ 20
    list_core=[]
    list_core_int=[]

    o_p_tol = 1. #use to define core boundaries

    for i_im in range(n_im):
        print
        print  'core : ',i_im
        print
    #test whether the core file already exists
        num_loc = i_im / freq_dump
        name='phi-core-new_'+str(num_loc)+'_'+str(o_p_tol)+'_'+str(num)+'.save'
#    print name, glob.glob(name)
#    pdb.set_trace() # debug mode     
        if ( len(glob.glob(name)) == 1):
            print 'core already processed, skip'
            continue


        center =  list_max_phi[i_im][5]
#    radius = 1. / 5000.

        test=True
    #initial guess on the radius
        radius = 1. / 10000. / 2. 
    
        ##we check whether the core has reached the boundary of the queried area 
        while(test):

            #count the number of times the decreasing gradient is not satisfied
            #core propagation stops when ncount > =ntol
            ncount=0

            region = Sphere(center,radius)
            filt_amr = RegionFilter(region,amr)

            cell_source = CellsToPoints(filt_amr)
            cells = cell_source.flatten()
            ncells = cells.npoints

            if(ncells ==0):
                print 'problem with the filter 0 cells found, this should not happen...'
                test=False
                continue

    
            ccenters = cells.points

        #get center position with respect to local minimum
            ccenters[:,0] = ccenters[:,0] - center[0]
            ccenters[:,1] = ccenters[:,1] - center[1]
            ccenters[:,2] = ccenters[:,2] - center[2]

        #potential
            phi = cells["phi"]

        #ordered with respect to potential
            ord_phi = np.argsort(phi)

            phi = phi[ord_phi]
            cells = cells.reorder_points(ord_phi)

        ##I am not quite sure why in the pointer logic it is necessary to change order here but it seems to be needed
            ccenters = ccenters[ord_phi]

        #mini avoids dividing by zero
            mini=1.e-20
            rad = np.sqrt(np.sum(ccenters*ccenters,1))
            rad[rad==0]=mini

            dx = cells.get_sizes()

##we search for the largest contours of potential which encompasses the local minima
##for that purpose we use the two neighbours of the current cells located in the radial 
#direction. The furthest must have a lower (in absolute value) potential while the 
#closest must belong to the core already

    ##first step
    ##loop over the cells sorted by potential and verify whether the neighboor of this cell that is
    ##further away has a potential that is "significantly" higher (in absolute value).

            points_ncell_av= np.zeros([ncells,3])
            points_ncell_am= np.zeros([ncells,3])
            np.copyto(points_ncell_av,ccenters) 
            np.copyto(points_ncell_am,ccenters) 



            for idim in range(3):
        #coef 0.9 is a compromise between being <  dx (the current cell is coarse and the neighbour is refined) and > sqrt(3)/2 (the
        #direction) is along the longest cube diagonal
                points_ncell_av[:,idim] +=  0.9*dx * points_ncell_av[:,idim] / rad + center[idim]
                points_ncell_am[:,idim] += -0.9*dx * points_ncell_am[:,idim] / rad + center[idim]



            psamp = PointSamplingProcessor(amr)
            neigh_cell_av = psamp.process(points_ncell_av) 
            neigh_cell_am = psamp.process(points_ncell_am) 

    #start looping over the cells sorted by potential to find the local minimum 
            for ic in range(0,ncells):
                if(rad[ic] == mini):
                    list_cell = [ic] # the local minimum is the first cell to be stored
                    nstart=ic+1
                    break

            for ic in range(nstart,ncells):

        #test whether the cell is not a local maximum
        #factor o_p_tol is there to eliminate small fluctuations
                if (cells["phi"][ic] < neigh_cell_av["phi"][ic] / o_p_tol):
            #now test whether the neighbour that is close to the minimum of 
            #potential is in the list of registered cells within the potential well
                    if(test_neigh(ic,list_cell,neigh_cell_am,cells)):
                        list_cell.append(ic)
            
                else:
            #test whether the local maximum is connected to the rest of the potential well
                    if(test_neigh(ic,list_cell,neigh_cell_am,cells)):
                 #end of the potential well
                        print 'invert gradient ',cells["phi"][ic] , neigh_cell_av["phi"][ic],ncount
                        ncount = ncount + 1
                        if(ncount >= ntol):
                            break

        #check whether the core has reached the boundary of the domain which has been queried
            if( max(rad[list_cell]+dx[list_cell])  >= radius ):
            #increase the size of the domain
                radius = radius * 2.
                print 'increase the size of the domain ',radius 
            #pdb.set_trace() # debug mode     
            else:    
            #no need to expand the domain further
                test=False


        if(ncells > 0):
            list_cell = np.asarray(list_cell)

            nc = ic+1

            x_v = ccenters[0:nc,0][list_cell]
            y_v = ccenters[0:nc,1][list_cell]
            z_v = ccenters[0:nc,2][list_cell]

            list_core.append([x_v,y_v,z_v,list_max_phi[i_im]])
            list_core_int.append([x_v,y_v,z_v,list_max_phi[i_im]])



#    if(fait_image):
        if( (np.mod(i_im+1,freq_dump) == 0 or i_im == nmax_phi-1) ):

            ii_l = i_im/freq_dump
            fnew = open('phi-core-new_'+str(ii_l)+'_'+str(o_p_tol)+'_'+str(num)+'.save','w')
            pickle.dump(list_core_int,fnew)
            fnew.close()
            list_core_int=[]
        
        ncell_core = x_v.shape[0]
        if(fait_image and ncell_core >= 100):
            #do some images
            make_image_core(amr,ro,center,radius,x_v,y_v,num,i_im=i_im,directory=directory)


     
    ##now read all output sort them and damp them in a single file
    def getkey2(item):
        key=item[13:15]
        if(key[1] =='_'):
            key='0'+key[0]
        return key


    list_core_gath=[]
    for lphi in sorted(glob.glob("phi-core-new*"+"_"+str(num)+".save"),key=getkey2):
        print 'reading and merging ',lphi 
        fcore = open(lphi,'r')
        list_core = pickle.load(fcore)
        fcore.close()
        list_core_gath += list_core

    list_core = list_core_gath

    name_phi_core='phi-core-tot'+'_'+str(o_p_tol)+'_'+str(num)+'.save'
    fnew = open(name_phi_core,'w')
    pickle.dump(list_core,fnew)
    fnew.close()


    return name_phi_core

########################################################################################
########################################################################################
###examine the phi-core which have been obtained in extract_core
## and redefine them using virial parameters. Then calculate various quantities
#fait_image : to make images
#directory : where the images are written
#exten : to specify an extension of the name of the results file that is being written
# should be used to avoid overwritting a former one for example
########################################################################################
########################################################################################
def extract_vir_core(path,num,name,fait_image=False,directory='IMAGES',exten='',gcomp=True):



    ro=pymses.RamsesOutput(path,num)

    amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"], grav_compat=gcomp)


    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)


    print 'name ',name

    fcore = open(name,'r')
    list_core = pickle.load(fcore)
    fcore.close()


    list_max_phi = [li[3] for li in list_core]



#    pdb.set_trace() # debug mode 

    ncore = len(list_core)

    mass_v = np.zeros(ncore)
    mass1_v = np.zeros(ncore)
    mass2_v = np.zeros(ncore)
    mass3_v = np.zeros(ncore)
    mass4_v = np.zeros(ncore)

    npoints_v = np.zeros(ncore)
    npoints1_v = np.zeros(ncore)
    npoints2_v = np.zeros(ncore)
    npoints3_v = np.zeros(ncore)
    npoints4_v = np.zeros(ncore)

    rho_max_v  = np.zeros(ncore)
    rho_mean_v  = np.zeros(ncore)


    print 'ncore ',ncore

    for icore in range(ncore):


        center =  list_max_phi[icore][5]

#    center_v[:,icore] = center

        npoints =  list_core[icore][0].size


        print 'icore ', icore, ' npoints ',npoints
    
        list_core[icore][0] = list_core[icore][0] + center[0]
        list_core[icore][1] = list_core[icore][1] + center[1]
        list_core[icore][2] = list_core[icore][2] + center[2]

        points = np.concatenate((list_core[icore][0],list_core[icore][1],list_core[icore][2])).reshape(3,npoints).T

        psamp = PointSamplingProcessor(amr)
        cells = psamp.process(points,add_level=True) 

#    dx = cells.get_sizes()

#    cells = sample_points(amr,points,add_level=True) 


        list_core[icore][0] = list_core[icore][0] - center[0]
        list_core[icore][1] = list_core[icore][1] - center[1]
        list_core[icore][2] = list_core[icore][2] - center[2]


    #min and max values of the potential
        phimax = np.max(cells["phi"])
        phimin = np.min(cells["phi"])


        dx = 0.5**cells["level"]

        dV = dx**3
        mass_code = np.sum(cells["rho"]*dV) 

    #min core velocity
#    vx = np.sum(cells["rho"]*cells["vel"][:,0]*dV) / mass_code
#    vy = np.sum(cells["rho"]*cells["vel"][:,1]*dV) / mass_code
#    vz = np.sum(cells["rho"]*cells["vel"][:,2]*dV) / mass_code


    #velocity of the local minimum of potential 
        vx = cells["vel"][0,0]
        vy = cells["vel"][0,1]
        vz = cells["vel"][0,2]


    #compute the velocity in the radial direction. It must be removed from the "support" since it is 
    #high when the gas collapses
        rad = np.sqrt( list_core[icore][0]**2 + list_core[icore][1]**2 + list_core[icore][2]**2 ) + 1.e-10 

        vr = ( (cells["vel"][:,0]-vx)*list_core[icore][0]+(cells["vel"][:,1]-vy)*list_core[icore][1]+(cells["vel"][:,2]-vz)*list_core[icore][2] ) / rad 

        ar = np.where( vr > 0.)
        vr[ar] = 0.

    # 1e-10 is there to avoid / 0

        vr_x = vr * list_core[icore][0] / rad 
        vr_y = vr * list_core[icore][1] / rad
        vr_z = vr * list_core[icore][2] / rad

        support_cv = ( (cells["vel"][:,0]-vx)**2 + (cells["vel"][:,1]-vy)**2 + (cells["vel"][:,2]-vz)**2 ) * cells["rho"] / 2. 
        support_cv += 1.5*cells["P"] 
        support_cv += 0.125*((cells["Br"][:,0]+cells["Bl"][:,0])**2+(cells["Br"][:,1]+cells["Bl"][:,1])**2+(cells["Br"][:,2]+cells["Bl"][:,2])**2)


        support2_cv = ( (cells["vel"][:,0]-vx -vr_x)**2 + (cells["vel"][:,1]-vy-vr_y)**2 + (cells["vel"][:,2]-vz-vr_z)**2 ) * cells["rho"] / 2. 
        support2_cv += 1.5*cells["P"] 
        support2_cv += 0.125*((cells["Br"][:,0]+cells["Bl"][:,0])**2+(cells["Br"][:,1]+cells["Bl"][:,1])**2+(cells["Br"][:,2]+cells["Bl"][:,2])**2)


    #lbox is there because of its definition in ramses  
        enerpot_cv = (phimax-cells["phi"])*cells["rho"]*lbox

        gdotr_cv = cells["rho"]*(cells["g"][:,0]*list_core[icore][0]+cells["g"][:,1]*list_core[icore][1]+cells["g"][:,2]*list_core[icore][2])*lbox


#    a_vir = np.where( support_cv < potent_cv)
    

#    cells_vir=cells[a_vir]

        sum_supp = np.cumsum(support_cv*dV)
        sum_supp2 = np.cumsum(support2_cv*dV)
        sum_gdotr = np.cumsum(gdotr_cv*dV)

        cells.add_scalars("supp",support_cv)
        cells.add_scalars("supp2",support2_cv)
        cells.add_scalars("epot",enerpot_cv)
        cells.add_scalars("gdotr",gdotr_cv)

        cells.add_scalars("sum_supp",sum_supp)
        cells.add_scalars("sum_supp2",sum_supp2)
        cells.add_scalars("sum_gdotr",sum_gdotr)


#    vir_core_func = lambda dset: (dset["supp"]<dset["epot"])
#    vir_core = PointFunctionFilter(vir_core_func,cells)

#    cells_vir = vir_core.flatten()


        rho_max_v[icore] = np.max(cells["rho"])

        mass_v[icore] = mass_code
        npoints_v[icore] = cells.npoints

        rho_mean_v[icore] = mass_v[icore] / np.sum(dV)


        a1 = np.where( support_cv <= enerpot_cv )
        mass1_v[icore] = np.sum((cells["rho"]*dV)[a1])   
        npoints1_v[icore] = len(a1[0])

        a2 = np.where( support2_cv <= enerpot_cv )
        mass2_v[icore] = np.sum((cells["rho"]*dV)[a2])   
        npoints2_v[icore] = len(a2[0])
    

        a_cum = np.where( sum_supp <= -sum_gdotr )
        npt = (a_cum[0].shape)[0]
        if(npt > 0):
            ind_max = a_cum[0][npt-1]
            mass3_v[icore] = np.sum((cells["rho"]*dV)[0:ind_max])   
            npoints3_v[icore] = ind_max


        a_cum = np.where( sum_supp2 <= -sum_gdotr)
        npt2 = (a_cum[0].shape)[0]
        if(npt2 > 0):
            ind_max2 = a_cum[0][npt2-1]
            mass4_v[icore] = np.sum((cells["rho"]*dV)[0:ind_max2])   
            npoints4_v[icore] = ind_max2



        if(fait_image):
            #do some images
            rad = np.sqrt(max(list_core[icore][0]**2 + list_core[icore][1]**2 + list_core[icore][2]**2))
            radius = max(rad,5.e-5)
            x_v = list_core[icore][0]
            y_v = list_core[icore][1]
            make_image_core(amr,ro,center,radius,x_v,y_v,num,i_im=icore,directory=directory)


#    pdb.set_trace() # debug mode 

    name_vir_core = 'phi-core-vir_'+str(num)+exten+'.save'

    core_vir = open(name_vir_core,'w')

    pickle.dump(rho_max_v,core_vir)
    pickle.dump(rho_mean_v,core_vir)
    pickle.dump(mass_v,core_vir)
    pickle.dump(mass1_v,core_vir)
    pickle.dump(mass2_v,core_vir)
    pickle.dump(mass3_v,core_vir)
    pickle.dump(mass4_v,core_vir)
    pickle.dump(npoints_v,core_vir)
    pickle.dump(npoints1_v,core_vir)
    pickle.dump(npoints2_v,core_vir)
    pickle.dump(npoints3_v,core_vir)
    pickle.dump(npoints4_v,core_vir)

    core_vir.close()



    return name_vir_core

########################################################################################
########################################################################################
###show various quantities for core which have been obtained in extract_core and extract_ir_core
## and redefine them using virial parameters. 
########################################################################################
########################################################################################
def show_vir_core(path,num,name,n_cut=100,gcomp=True):
    ro=pymses.RamsesOutput(path,num)
    amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"],grav_compat=gcomp)
#amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"])



##units
#Ms=2.e33   #Ms in g
#pc=3.08e18 #pc in cm
#mp=1.4*1.66e-24 #mass per particles in g
#G=6.8e-8 #G in cgs


    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    mu  =  1.4
    mp =  mu * 1.660531e-24  #n gramme
    G = 6.7e-8
    kbol  =  1.38062e-16   # erg/degre
    pc=3.08e18 #cm in pc
    Ms = 2.e33 #solar mass in g
    Myrs =  3.15e13 # time in Myrs
#scale_n converts particle density from user unit into number of
#particle per cm^-3
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




    vir_core = open(name,'r')

    rho_max_v = pickle.load(vir_core)
    rho_mean_v = pickle.load(vir_core)

    mass_v = pickle.load(vir_core)*scale_mass/Ms
    mass1_v = pickle.load(vir_core)*scale_mass/Ms
    mass2_v = pickle.load(vir_core)*scale_mass/Ms
    mass3_v = pickle.load(vir_core)*scale_mass/Ms
    mass4_v = pickle.load(vir_core)*scale_mass/Ms

    npoints_v = pickle.load(vir_core)
    npoints1_v = pickle.load(vir_core)
    npoints2_v = pickle.load(vir_core)
    npoints3_v = pickle.load(vir_core)
    npoints4_v = pickle.load(vir_core)

    vir_core.close()



    mask0 = [npoints_v < 10]
    print 'point0 ',np.sum(mask0)

    mask1 = [npoints1_v < 10]
    print 'point1 ',np.sum(mask1)

    mask2 = [npoints2_v < 10]
    print 'point2 ',np.sum(mask2)
    
    mask3 = [npoints3_v < 10]
    print 'point3 ',np.sum(mask3)

    mask4 = [npoints4_v < 10]
    print 'point4 ',np.sum(mask4)


#    n_cut=10

    log_mmin =-2.
    log_mmax = 3.

    nbin=50.

    width = (log_mmax - log_mmin) / nbin

    P.clf() 

    a = np.where(npoints_v > n_cut)

    hist_mass , hist_edges = P.histogram(np.log10(mass_v[a]),bins=nbin ,range=(log_mmin,log_mmax))

    P.bar(hist_edges[:-1],hist_mass,width)

    P.yscale('log', nonposy='clip') 

    
    P.xlabel(r'$\log(M) \, (M_s)$')
    P.ylabel(r'$\log(dN/dlogM)$')

    P.savefig('hist_mass.ps')
    P.savefig('hist_mass.jpeg')



    P.clf()

    a = np.where(npoints_v != 0)

    P.plot(np.log10(npoints_v[a]),np.log10(mass_v[a]),'o')

    P.xlabel(r'$\log N$')
    P.ylabel(r'$\log(M) \, (M_s)$')

    P.savefig('point_mass.ps')
    P.savefig('point_mass.jpeg')



    P.clf() 

    a = np.where(npoints1_v > n_cut)

    hist_mass , hist_edges = P.histogram(np.log10(mass1_v[a]),bins=nbin ,range=(log_mmin,log_mmax))


    P.bar(hist_edges[:-1],hist_mass,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\log(M) \, (M_s)$')
    P.ylabel(r'$\log(dN/dlogM)$')

    P.savefig('hist_mass1.ps')
    P.savefig('hist_mass1.jpeg')


    P.clf()

    a = np.where(npoints1_v != 0)

    P.plot(np.log10(npoints1_v[a]),np.log10(mass1_v[a]),'o')

    P.xlabel(r'$\log N$')
    P.ylabel(r'$\log(M) \, (M_s)$')

    P.savefig('point_mass1.ps')
    P.savefig('point_mass1.jpeg')



    P.clf() 

    a = np.where(npoints2_v > n_cut)

    hist_mass , hist_edges = P.histogram(np.log10(mass2_v[a]),bins=nbin ,range=(log_mmin,log_mmax))


    P.bar(hist_edges[:-1],hist_mass,width)

    P.yscale('log', nonposy='clip') 


    P.xlabel(r'$\log(M) \, (M_s)$')
    P.ylabel(r'$\log(dN/dlogM)$')

    P.savefig('hist_mass2.ps')
    P.savefig('hist_mass2.jpeg')



    P.clf()

    a = np.where(npoints2_v != 0)

    P.plot(np.log10(npoints2_v[a]),np.log10(mass2_v[a]),'o')

    P.xlabel(r'$\log N$')
    P.ylabel(r'$\log(M) \, (M_s)$')

    P.savefig('point_mass2.ps')
    P.savefig('point_mass2.jpeg')


    P.clf() 

    a = np.where(npoints3_v > n_cut)

    hist_mass , hist_edges = P.histogram(np.log10(mass3_v[a]),bins=nbin ,range=(log_mmin,log_mmax))


    P.bar(hist_edges[:-1],hist_mass,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\log(M) \, (M_s)$')
    P.ylabel(r'$\log(dN/dlogM)$')

    P.savefig('hist_mass3.ps')
    P.savefig('hist_mass3.jpeg')

    P.clf()

    a = np.where(npoints3_v != 0)

    P.plot(np.log10(npoints3_v[a]),np.log10(mass3_v[a]),'o')

    P.savefig('point_mass3.ps')
    P.savefig('point_mass3.jpeg')



    P.clf() 

    a = np.where(npoints4_v > n_cut)

    hist_mass , hist_edges = P.histogram(np.log10(mass4_v[a]),bins=nbin ,range=(log_mmin,log_mmax))


    P.bar(hist_edges[:-1],hist_mass,width)

    P.yscale('log', nonposy='clip') 

    P.xlabel(r'$\log(M) \, (M_s)$')
    P.ylabel(r'$\log(dN/dlogM)$')

    P.savefig('hist_mass4.ps')
    P.savefig('hist_mass4.jpeg')


    P.clf()

    a = np.where(npoints4_v != 0)

    P.plot(np.log10(npoints4_v[a]),np.log10(mass4_v[a]),'o')

    P.xlabel(r'$\log N$')
    P.ylabel(r'$\log(M) \, (M_s)$')

    P.savefig('point_mass4.ps')
    P.savefig('point_mass4.jpeg')

    return
########################################################################################
########################################################################################
##make images for cores
##amr and ro pymses stuff
##center 3D array for coordinates center
##radius the size of the region to be plotted (between 0 and 1)
##x_v and y_v coordinates of the points to be overplotted
##num output number
##i_im image number
##directory : the directory in which the images are being written 
########################################################################################
########################################################################################
def  make_image_core(amr,ro,center,radius,x_v,y_v,num,directory,i_im=0,force=False,coldens=True):


    name = directory+'/coldens_z'+'_'+str(i_im)+'_'+str(num)+'.jpeg'
    if (len(glob.glob(name))==1 and not force):
        return 

    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)
    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

    lbox_cm = lbox * pc #lbox in cm

    cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=512)

    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_col = np.log10(datamap.map.T*lbox_cm)


    P.clf() 
    im = P.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    P.xlabel('$x$ (pc)')     
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(N) \, cm^{-2}$')                                          

    P.savefig(directory+'/coldens2_z'+'_'+str(i_im)+'_'+str(num)+'.jpeg')



    P.clf() 
    im = P.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    P.xlabel('$x$ (pc)')     
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(N) \, cm^{-2}$')                                          
#    P.scatter((x_v+center[0])*lbox,(y_v+center[1])*lbox,c='r',marker='x',s=10)

    P.scatter((x_v+center[0])*lbox,(y_v+center[1])*lbox,c='r',marker='x',s=10)

    P.savefig(directory+'/coldens2_z'+'_'+str(i_im)+'_nscat_'+str(num)+'.jpeg')






    #do an integration and mass weight the results
    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_rho = datamap.map.T


    vx_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,vx_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_vx = datamap.map.T / map_rho


    vy_op = ScalarOperator(lambda dset: dset["vel"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,vy_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_vy = datamap.map.T / map_rho

    map_vx_red=map_vx[::20,::20]
    map_vy_red=map_vy[::20,::20]

    P.clf() 
    im = P.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_vx_red.shape[0]
    ny = map_vy_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_vx_red,map_vy_red)

#    P.title(titre)
    P.xlabel('$x$ (pc)')          
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(N) \, cm^{-2}$')                                          

    P.savefig(directory+'/coldens2_vel_z'+'_'+str(i_im)+'_'+str(num)+'.jpeg')
    P.savefig(directory+'/coldens2_vel_z'+'_'+str(i_im)+'_'+str(num)+'.ps')
    P.savefig(directory+'/coldens2_vel_z'+'_'+str(i_im)+'_'+str(num)+'.pdf')



    if( coldens):
        return 


    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    dmap_rho = slicing.SliceMap(amr, cam_z, rho_op, z=0.) #,verbose=False)
    map_rho = np.log10(dmap_rho.map)

    map_rho=map_rho.T

    vx_op = ScalarOperator(lambda dset: dset["vel"][:,0] ,  ro.info["unit_velocity"])
    dmap_vx = slicing.SliceMap(amr, cam_z, vx_op, z=0.) #,verbose=False)
    map_vx_red=dmap_vx.map[::20,::20]

    map_vx_red=map_vx_red.T

    vy_op = ScalarOperator(lambda dset: dset["vel"][:,1],  ro.info["unit_velocity"])
    dmap_vy = slicing.SliceMap(amr, cam_z, vy_op, z=0.) #,verbose=False)
    map_vy_red=dmap_vy.map[::20,::20]

    map_vy_red=map_vy_red.T




    P.clf() 
    im = P.imshow(map_rho,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_vx_red.shape[0]
    ny = map_vx_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_vx_red,map_vy_red)

    P.xlabel('$x$ (pc)')          
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(\rho)$')                                                                                                             
    P.scatter((x_v+center[0])*lbox,(y_v+center[1])*lbox,c='r',marker='x',s=10)

    P.savefig(directory+'/rho_z'+'_'+str(i_im)+'_'+str(num)+'.jpeg')



#    phi_op = ScalarOperator(lambda dset: np.log10(abs(dset["phi"])) ,ro.info["unit_density"])
#    dmap_phi = slicing.SliceMap(amr, cam_z, phi_op, z=0.) #,verbose=False)
#    map_phi = dmap_phi.map


#    map_phi=map_phi.T

#    P.clf() 
#    im = P.imshow(map_phi,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

#    P.xlabel('$x$ (pc)')                                                     
#    P.ylabel('$y$ (pc)')
#    cbar = P.colorbar(im)                                                                                                                     
#    cbar.set_label(r'$log(\phi)$')                              
                                                                              
#    P.scatter((x_v+center[0])*lbox,(y_v+center[1])*lbox,c='r',marker='x',s=10) 

#    P.savefig(directory+'/phi_scat'+'_'+str(i_im)+'_'+str(num)+'.jpeg')

    return 

########################################################################################
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
########################################################################################
def make_image(amr,ro,center,radius,num,path,x_v=None,y_v=None,z_v=None,mass_v=None,i_im=0,force=False,path_out=None,col_dens_only=False,map_size=512,vel_red=20,im_pos=0.,tag='',ps=False,all_sink=True):

##    #must be reorganised...
#    Myr=3.15e7*1.e6
#    mu  =  1.4
#    mp =  mu * 1.660531e-24  #n gramme
#    G = 6.7e-8
#    scale_n = 1.
#    scale_d = scale_n*mp  #code density unit is in particle/cc
#    scale_t = 1.0/np.sqrt(G*scale_d)

    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

    lbox_cm = lbox * pc #lbox in cm

    time=ro.info['time'] #time in codeunits 
    time = time * scale_t / Myr


    titre = 't='+ str(time)[0:5] +' (Myr)'

    if( path_out is not None):
        directory = path_out
    else:
        directory = path

    name = directory+'/coldens_z'+'_'+str(i_im)+'_'+format(num,'05')+'.jpeg'
    if (len(glob.glob(name))==1 and not force):
        return 
    
    ## figures along the z-axis

    cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=map_size)

    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_col = np.log10(datamap.map.T*lbox_cm)

    P.clf() 
    im = P.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$x$ (pc)')     
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(N) \, cm^{-2}$')                                          
    P.savefig(directory+'/coldens_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
        P.savefig(directory+'/coldens_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')

    if(x_v is not None and y_v is not None):
        P.clf() 
        im = P.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   
#        P.title=(titre,horizontalalignment='right')
        P.title=(titre)
        P.xlabel('$x$ (pc)')     
        P.ylabel('$y$ (pc)')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(N) \, cm^{-2}$')                                          

#        P.scatter((x_v+center[0])*lbox,(y_v+center[1])*lbox,c='r',marker='x',s=10)

        if(all_sink):
            mask = np.where( (x_v/lbox >=-radius+center[0]) & (x_v/lbox <= radius+center[0]) &  (y_v/lbox >= -radius+center[1]) & (y_v/lbox <= radius+center[1]) == True)

        else:
            mask = np.where( (x_v/lbox >=-radius+center[0]) & (x_v/lbox <= radius+center[0]) &  (y_v/lbox >= -radius+center[1]) & (y_v/lbox <= radius+center[1]) &  (z_v/lbox >=-radius+center[2]) & (z_v/lbox <= radius+center[2]) == True)


        siz_sym=np.sqrt(mass_v[mask])*20.
        P.scatter(x_v[mask],y_v[mask],marker='o',s=siz_sym,color='red') #,alpha=0.1)

        P.savefig(directory+'/coldens_z'+'_'+tag+str(i_im)+'_nscat_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/coldens_z'+'_'+tag+str(i_im)+'_nscat_'+format(num,'05')+'.ps')


    if(not col_dens_only):


        rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
        dmap_rho = slicing.SliceMap(amr, cam_z, rho_op, z=0.) #,verbose=False)
        map_rho = np.log10(dmap_rho.map)

        map_rho=map_rho.T

        vx_op = ScalarOperator(lambda dset: dset["vel"][:,0] ,  ro.info["unit_velocity"])
        dmap_vx = slicing.SliceMap(amr, cam_z, vx_op, z=0.) #,verbose=False)
        map_vx_red=dmap_vx.map[::vel_red,::vel_red]

        map_vx_red=map_vx_red.T

        vy_op = ScalarOperator(lambda dset: dset["vel"][:,1],  ro.info["unit_velocity"])
        dmap_vy = slicing.SliceMap(amr, cam_z, vy_op, z=0.) #,verbose=False)
        map_vy_red=dmap_vy.map[::vel_red,::vel_red]

        map_vy_red=map_vy_red.T

        P.clf() 
        im = P.imshow(map_rho,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

        nx = map_vx_red.shape[0]
        ny = map_vx_red.shape[1]
        vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
        vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
        xx,yy = np.meshgrid(vec_x,vec_y)

        P.quiver(xx,yy,map_vx_red,map_vy_red)

#        P.title=(titre,horizontalalignment='right')
        P.title=(titre)
        P.xlabel('$x$ (pc)')          
        P.ylabel('$y$ (pc)')
        cbar = P.colorbar(im)                                                                                                                     
        cbar.set_label(r'$log(n) \, (cm^{-3})$')                                                                                                   

        P.savefig(directory+'/rho_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/rho_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')




        P_op = ScalarOperator(lambda dset: dset["P"] ,  ro.info["unit_pressure"])
        dmap_P = slicing.SliceMap(amr, cam_z, P_op, z=0.) #,verbose=False)

        P.clf() 
        map_T = np.log10(dmap_P.map.T/dmap_rho.map.T*scale_T2)

        im = P.imshow(map_T,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$x$ (pc)')        
        P.ylabel('$y$ (pc)')
        cbar = P.colorbar(im)   
        cbar.set_label(r'$log(T) \, (K)$')                                                                                                   

        P.savefig(directory+'/T_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/T_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')



        level_op = MaxLevelOperator() #ScalarOperator(lambda dset: dset["level"] ,  ro.info["unit_pressure"])
        amr.set_read_levelmax(20)
        rt = raytracing.RayTracer(amr,ro.info,level_op)
        datamap = rt.process(cam_z, surf_qty=True)
        map_level = (datamap.map.T)

#        dmap_level = slicing.SliceMap(amr, cam_z, level_op, z=0.) #,verbose=False)
#        map_level = dmap_level.map.T

        P.clf() 

        im = P.imshow(map_level,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$x$ (pc)')        
        P.ylabel('$y$ (pc)')
        cbar = P.colorbar(im)   
        cbar.set_label(r'level')

        P.savefig(directory+'/level_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/level_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')




        cpu_op = ScalarOperator(lambda dset: dset.icpu*(np.ones(dset["P"].shape)) ,  ro.info["unit_pressure"])
        rt = raytracing.RayTracer(amr,ro.info,cpu_op)
        datamap = rt.process(cam_z, surf_qty=True)
        map_cpu = datamap.map.T

#        dmap_cpu = slicing.SliceMap(amr, cam_z, cpu_op, z=0.) #,verbose=False)
#        map_cpu = dmap_cpu.map.T


        P.clf() 

        im = P.imshow(map_cpu,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$x$ (pc)')        
        P.ylabel('$y$ (pc)')
        cbar = P.colorbar(im)   
        cbar.set_label(r'level')

        P.savefig(directory+'/cpu_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/cpu_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')




    ## figures along the y-axis

    cam_y= Camera(center=center,line_of_sight_axis='y',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='x',map_max_size=map_size)

    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_y, surf_qty=True)
#    map_col = np.log10(datamap.T*lbox_cm)
    map_col = np.log10(datamap.map*lbox_cm)

    P.clf() 
    im = P.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$x$ (pc)')     
    P.ylabel('$z$ (pc)')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(N) \, cm^{-2}$')                                          
    P.savefig(directory+'/coldens_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/coldens_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')

    if(x_v is not None and z_v is not None):
        P.clf() 
        im = P.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   
        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$x$ (pc)')     
        P.ylabel('$z$ (pc)')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(N) \, cm^{-2}$')                                          




        if(all_sink):
            mask = np.where( (x_v/lbox >= -radius+center[0]) & (x_v/lbox <= radius+center[0]) &  (z_v/lbox >= -radius+center[2]) & (z_v/lbox <= radius+center[2]) == True)

        else:
            mask = np.where( (x_v/lbox >=-radius+center[0]) & (x_v/lbox <= radius+center[0]) &  (y_v/lbox >= -radius+center[1]) & (y_v/lbox <= radius+center[1]) &  (z_v/lbox >=-radius+center[2]) & (z_v/lbox <= radius+center[2]) == True)



        siz_sym=np.sqrt(mass_v[mask])*20.
        P.scatter(x_v[mask],z_v[mask],marker='o',s=siz_sym,color='red')#,alpha=0.1)
#        P.scatter(x_v,z_v,c='r',marker='x',s=10)

        P.savefig(directory+'/coldens_y'+'_'+tag+str(i_im)+'_nscat_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/coldens_y'+'_'+tag+str(i_im)+'_nscat_'+format(num,'05')+'.ps')




    if(not col_dens_only):

        rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
        dmap_rho = slicing.SliceMap(amr, cam_y, rho_op, z=0.) #,verbose=False)
        map_rho = np.log10(dmap_rho.map)

#        map_rho = map_rho.T

        vx_op = ScalarOperator(lambda dset: dset["vel"][:,0] ,  ro.info["unit_velocity"])
        dmap_vx = slicing.SliceMap(amr, cam_y, vx_op, z=0.) #,verbose=False)
        map_vx_red=dmap_vx.map[::vel_red,::vel_red]

#        map_vx_red = map_vx_red.T

        vz_op = ScalarOperator(lambda dset: dset["vel"][:,2],  ro.info["unit_velocity"])
        dmap_vz = slicing.SliceMap(amr, cam_y, vz_op, z=0.) #,verbose=False)
        map_vz_red=dmap_vz.map[::vel_red,::vel_red]

#        map_vz_red = map_vz_red.T


        P.clf() 
        im = P.imshow(map_rho,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   

        nx = map_vx_red.shape[0]
        nz = map_vx_red.shape[1]
        vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
        vec_z = (np.arange(nz)*2./nz*radius - radius + center[2] + radius/nz)*lbox
        xx,zz = np.meshgrid(vec_x,vec_z)

        P.quiver(xx,zz,map_vx_red,map_vz_red)

        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$x$ (pc)')          
        P.ylabel('$z$ (pc)')
        cbar = P.colorbar(im)                                                                                                                     
        cbar.set_label(r'$log(n) \, (cm^{-3})$')                                                                                                   

        P.savefig(directory+'/rho_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/rho_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')



        P_op = ScalarOperator(lambda dset: dset["P"] ,  ro.info["unit_pressure"])
        dmap_P = slicing.SliceMap(amr, cam_y, P_op, z=0.) #,verbose=False)

        P.clf() 
#        map_T = np.log10(dmap_P.map.T/dmap_rho.map.T*scale_T2)
        map_T = np.log10(dmap_P.map/dmap_rho.map*scale_T2)

        im = P.imshow(map_T,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   

        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$x$ (pc)')        
        P.ylabel('$z$ (pc)')
        cbar = P.colorbar(im)   
        cbar.set_label(r'$log(T) \, (K)$')                                                                                                   

        P.savefig(directory+'/T_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/T_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')




    ## figures along the x-axis

    cam_x= Camera(center=center,line_of_sight_axis='x',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='z',map_max_size=map_size)

    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_x, surf_qty=True)
    map_col = np.log10(datamap.map.T*lbox_cm)

    P.clf() 
    im = P.imshow(map_col,extent=[(-radius+center[1])*lbox,(radius+center[1])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$y$ (pc)')     
    P.ylabel('$z$ (pc)')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(N) \, cm^{-2}$')                                          
    P.savefig(directory+'/coldens_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/coldens_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')


    if(y_v is not None and z_v is not None):
        P.clf() 
        im = P.imshow(map_col,extent=[(-radius+center[1])*lbox,(radius+center[1])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   
        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$y$ (pc)')     
        P.ylabel('$z$ (pc)')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(N) \, cm^{-2}$')                                          



        if(all_sink):
            mask = np.where( (y_v/lbox >= -radius+center[1]) & (y_v/lbox <= radius+center[1]) &  (z_v/lbox >= -radius+center[2]) & (z_v/lbox <= radius+center[2]) == True)
        else:
            mask = np.where( (x_v/lbox >=-radius+center[0]) & (x_v/lbox <= radius+center[0]) &  (y_v/lbox >= -radius+center[1]) & (y_v/lbox <= radius+center[1]) &  (z_v/lbox >=-radius+center[2]) & (z_v/lbox <= radius+center[2]) == True)


        siz_sym=np.sqrt(mass_v[mask])*20.
        P.scatter(y_v[mask],z_v[mask],marker='o',s=siz_sym,color='red')#,alpha=0.1)
#        P.scatter(y_v,z_v,c='r',marker='x',s=10)

        P.savefig(directory+'/coldens_x'+'_'+tag+str(i_im)+'_nscat_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/coldens_x'+'_'+tag+str(i_im)+'_nscat_'+format(num,'05')+'.ps')




    if(not col_dens_only):

        rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
        dmap_rho = slicing.SliceMap(amr, cam_x, rho_op, z=0.) #,verbose=False)
        map_rho = np.log10(dmap_rho.map)

        map_rho = map_rho.T


        vy_op = ScalarOperator(lambda dset: dset["vel"][:,1] ,  ro.info["unit_velocity"])
        dmap_vy = slicing.SliceMap(amr, cam_x, vy_op, z=0.) #,verbose=False)
        map_vy_red=dmap_vy.map[::vel_red,::vel_red]

        map_vy_red = map_vy_red.T


        vz_op = ScalarOperator(lambda dset: dset["vel"][:,2],  ro.info["unit_velocity"])
        dmap_vz = slicing.SliceMap(amr, cam_x, vz_op, z=0.) #,verbose=False)
        map_vz_red=dmap_vz.map[::vel_red,::vel_red]

        map_vz_red = map_vz_red.T
        P.clf() 
        im = P.imshow(map_rho,extent=[(-radius+center[1])*lbox,(radius+center[1])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   

        ny = map_vx_red.shape[0]
        nz = map_vx_red.shape[1]
        vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/ny)*lbox
        vec_z = (np.arange(nz)*2./ny*radius - radius + center[2] + radius/nz)*lbox
        yy,zz = np.meshgrid(vec_y,vec_z)

        P.quiver(yy,zz,map_vy_red,map_vz_red)

        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$y$ (pc)')          
        P.ylabel('$z$ (pc)')
        cbar = P.colorbar(im)                                                                                                                     
        cbar.set_label(r'$log(n) \, (cm^{-3})$')                                                                                                   

        P.savefig(directory+'/rho_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/rho_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')





        P_op = ScalarOperator(lambda dset: dset["P"] ,  ro.info["unit_pressure"])
        dmap_P = slicing.SliceMap(amr, cam_x, P_op, z=0.) #,verbose=False)

        P.clf() 
        map_T = np.log10(dmap_P.map/dmap_rho.map*scale_T2)

        map_T = map_T.T

        im = P.imshow(map_T,extent=[(-radius+center[1])*lbox,(radius+center[1])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   

        P.title=(titre)
#        P.title=(titre,horizontalalignment='right')
        P.xlabel('$y$ (pc)')        
        P.ylabel('$z$ (pc)')
        cbar = P.colorbar(im)   
        cbar.set_label(r'$log(T) \, (K)$')                                                                                                   

        P.savefig(directory+'/T_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
        if(ps):
            P.savefig(directory+'/T_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')


    return 


########################################################################################
########################################################################################
##make images of yje magnetic field
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
########################################################################################
def make_image_B(amr,ro,center,radius,num,path,x_v=None,y_v=None,z_v=None,mass_v=None,i_im=0,force=False,path_out=None,col_dens_only=False,map_size=512,vel_red=20,im_pos=0.,tag='',ps=False):


    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

    lbox_cm = lbox * pc #lbox in cm

    time=ro.info['time'] #time in codeunits 
    time = time * scale_t / Myr

#    titre = 't='+ str(time)[0:5] +' (Myr)'

    titre = 't='+ str(time)[0:5] +' (Myr)'

    if( path_out is not None):
        directory = path_out
    else:
        directory = path

    name = directory+'/B_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg'
    if (len(glob.glob(name))==1 and not force):
        return 
    
    ## figures along the z-axis
    cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=map_size)


    #do an integration and mass weight the results
    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_rho = datamap.map.T


    Bx_op = ScalarOperator(lambda dset: dset["Br"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,Bx_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_Bx = datamap.map.T / map_rho


    By_op = ScalarOperator(lambda dset: dset["Br"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,By_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_By = datamap.map.T / map_rho


    Bz_op = ScalarOperator(lambda dset: dset["Br"][...,2]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,Bz_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_Bz = datamap.map.T / map_rho


    map_B = np.log10(np.sqrt(map_Bx**2+map_By**2+map_Bz**2)*scale_mag*microG )

    map_norm = np.sqrt(map_Bx**2+map_By**2)
    map_Bx = map_Bx / map_norm
    map_By = map_By / map_norm

    map_Bx_red=map_Bx[::vel_red,::vel_red]
    map_By_red=map_By[::vel_red,::vel_red]

    P.clf() 
    im = P.imshow(map_B,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_Bx_red.shape[0]
    ny = map_By_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_Bx_red,map_By_red)

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$x$ (pc)')          
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(B) \, (\mu G)$')                                                                                                   

    P.savefig(directory+'/B_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/B_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')



#now do a cut
    Bx_op = ScalarOperator(lambda dset: dset["Br"][:,0] ,  ro.info["unit_velocity"])
    dmap_Bx = slicing.SliceMap(amr, cam_z, Bx_op, z=0.) #,verbose=False)
    map_Bx = dmap_Bx.map.T


    By_op = ScalarOperator(lambda dset: dset["Br"][:,1] ,  ro.info["unit_velocity"])
    dmap_By = slicing.SliceMap(amr, cam_z, By_op, z=0.) #,verbose=False)
    map_By = dmap_By.map.T


    Bz_op = ScalarOperator(lambda dset: dset["Br"][:,2] ,  ro.info["unit_velocity"])
    dmap_Bz = slicing.SliceMap(amr, cam_z, Bz_op, z=0.) #,verbose=False)
    map_Bz = dmap_Bz.map.T


    map_B = np.log10(np.sqrt(map_Bx**2+map_By**2+map_Bz**2)*scale_mag*microG )

    map_norm = np.sqrt(map_Bx**2+map_By**2)
    map_Bx = map_Bx / map_norm
    map_By = map_By / map_norm

    map_Bx_red=map_Bx[::vel_red,::vel_red]
    map_By_red=map_By[::vel_red,::vel_red]

    P.clf() 
    im = P.imshow(map_B,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_Bx_red.shape[0]
    ny = map_By_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_Bx_red,map_By_red)

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$x$ (pc)')          
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(B) \, (\mu G)$')                                                                                                   

    P.savefig(directory+'/B_z_cut'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/B_z_cut'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')

    ## figures along the y-axis
    cam_y= Camera(center=center,line_of_sight_axis='y',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='x',map_max_size=map_size)


    #do an integration and mass weight the results
    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_y, surf_qty=True)
#    map_rho = datamap.T
    map_rho = datamap.map


    Bx_op = ScalarOperator(lambda dset: dset["Br"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,Bx_op)
    datamap = rt.process(cam_y, surf_qty=True)
#    map_Bx = datamap.T / map_rho
    map_Bx = datamap.map / map_rho


    By_op = ScalarOperator(lambda dset: dset["Br"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,By_op)
    datamap = rt.process(cam_y, surf_qty=True)
#    map_By = datamap.T / map_rho
    map_By = datamap.map / map_rho


    Bz_op = ScalarOperator(lambda dset: dset["Br"][...,2]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,Bz_op)
    datamap = rt.process(cam_y, surf_qty=True)
    map_Bz = datamap.map / map_rho


    map_B = np.log10(np.sqrt(map_Bx**2+map_By**2+map_Bz**2)*scale_mag*microG )

    map_norm = np.sqrt(map_Bx**2+map_Bz**2)
    map_Bx = map_Bx / map_norm
    map_Bz = map_Bz / map_norm

    map_Bx_red=map_Bx[::vel_red,::vel_red]
    map_Bz_red=map_Bz[::vel_red,::vel_red]

    P.clf() 
    im = P.imshow(map_B,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_Bx_red.shape[0]
    ny = map_Bz_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_Bx_red,map_Bz_red)

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$x$ (pc)')          
    P.ylabel('$z$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(B) \, (\mu G)$')                                                                                                   

    P.savefig(directory+'/B_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/B_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')



#now do a cut
    Bx_op = ScalarOperator(lambda dset: dset["Br"][:,0] ,  ro.info["unit_velocity"])
    dmap_Bx = slicing.SliceMap(amr, cam_y, Bx_op, z=0.) #,verbose=False)
    map_Bx = dmap_Bx.map#.T


    By_op = ScalarOperator(lambda dset: dset["Br"][:,1] ,  ro.info["unit_velocity"])
    dmap_By = slicing.SliceMap(amr, cam_y, By_op, z=0.) #,verbose=False)
    map_By = dmap_By.map#.T


    Bz_op = ScalarOperator(lambda dset: dset["Br"][:,2] ,  ro.info["unit_velocity"])
    dmap_Bz = slicing.SliceMap(amr, cam_y, Bz_op, z=0.) #,verbose=False)
    map_Bz = dmap_Bz.map#.T


    map_B = np.log10(np.sqrt(map_Bx**2+map_By**2+map_Bz**2)*scale_mag*microG )

    map_norm = np.sqrt(map_Bx**2+map_Bz**2)
    map_Bx = map_Bx / map_norm
    map_Bz = map_Bz / map_norm

    map_Bx_red=map_Bx[::vel_red,::vel_red]
    map_Bz_red=map_Bz[::vel_red,::vel_red]

    P.clf() 
    im = P.imshow(map_B,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_Bx_red.shape[0]
    ny = map_Bz_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_Bx_red,map_Bz_red)

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$x$ (pc)')          
    P.ylabel('$z$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(B) \, (\mu G)$')                                                                                                   

    P.savefig(directory+'/B_y_cut'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/B_y_cut'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')




    ## figures along the x-axis
    cam_x= Camera(center=center,line_of_sight_axis='x',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='z',map_max_size=map_size)


    #do an integration and mass weight the results
    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_x, surf_qty=True)
    map_rho = datamap.map.T


    Bx_op = ScalarOperator(lambda dset: dset["Br"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,Bx_op)
    datamap = rt.process(cam_x, surf_qty=True)
    map_Bx = datamap.map.T / map_rho


    By_op = ScalarOperator(lambda dset: dset["Br"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,By_op)
    datamap = rt.process(cam_x, surf_qty=True)
    map_By = datamap.map.T / map_rho


    Bz_op = ScalarOperator(lambda dset: dset["Br"][...,2]*dset["rho"] ,  ro.info["unit_velocity"])
    rt = raytracing.RayTracer(amr,ro.info,Bz_op)
    datamap = rt.process(cam_x, surf_qty=True)
    map_Bz = datamap.map.T / map_rho


    map_B = np.log10(np.sqrt(map_Bx**2+map_By**2+map_Bz**2)*scale_mag*microG )

    map_norm = np.sqrt(map_By**2+map_Bz**2)
    map_By = map_By / map_norm
    map_Bz = map_Bz / map_norm

    map_By_red=map_By[::vel_red,::vel_red]
    map_Bz_red=map_Bz[::vel_red,::vel_red]

    P.clf() 
    im = P.imshow(map_B,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_Bx_red.shape[0]
    ny = map_By_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_By_red,map_Bz_red)

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$y$ (pc)')          
    P.ylabel('$z$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(B) \, (\mu G)$')                                                                                                   

    P.savefig(directory+'/B_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/B_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')



#now do a cut
    Bx_op = ScalarOperator(lambda dset: dset["Br"][:,0] ,  ro.info["unit_velocity"])
    dmap_Bx = slicing.SliceMap(amr, cam_x, Bx_op, z=0.) #,verbose=False)
    map_Bx = dmap_Bx.map.T


    By_op = ScalarOperator(lambda dset: dset["Br"][:,1] ,  ro.info["unit_velocity"])
    dmap_By = slicing.SliceMap(amr, cam_x, By_op, z=0.) #,verbose=False)
    map_By = dmap_By.map.T


    Bz_op = ScalarOperator(lambda dset: dset["Br"][:,2] ,  ro.info["unit_velocity"])
    dmap_Bz = slicing.SliceMap(amr, cam_x, Bz_op, z=0.) #,verbose=False)
    map_Bz = dmap_Bz.map.T


    map_B = np.log10(np.sqrt(map_Bx**2+map_By**2+map_Bz**2)*scale_mag*microG )

    map_norm = np.sqrt(map_By**2+map_Bz**2)
    map_By = map_By / map_norm
    map_Bz = map_Bz / map_norm

    map_By_red=map_By[::vel_red,::vel_red]
    map_Bz_red=map_Bz[::vel_red,::vel_red]

    P.clf() 
    im = P.imshow(map_B,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    nx = map_By_red.shape[0]
    ny = map_Bz_red.shape[1]
    vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox
    xx,yy = np.meshgrid(vec_x,vec_y)

    P.quiver(xx,yy,map_By_red,map_Bz_red)

    P.title=(titre)
#    P.title=(titre,horizontalalignment='right')
    P.xlabel('$y$ (pc)')          
    P.ylabel('$z$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(B) \, (\mu G)$')                                                                                                   

    P.savefig(directory+'/B_x_cut'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.jpeg')
    if(ps):
            P.savefig(directory+'/B_x_cut'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.ps')




    return 

########################################################################################
#estimate the 3 components of B in the equatorial plane (I&H2017 paper)
def estimate_Bx_By_Bz_z0(path,num):


    ro=pymses.RamsesOutput(path,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)


    amr = ro.amr_source(["Br"])


    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

    lbox_cm = lbox * pc #lbox in cm

    time=ro.info['time'] #time in codeunits 
    time = time * scale_t / Myr

    titre = 't='+ str(time)[0:5] +' (Myr)'


    center=[0.5,0.5,0.5]
    radius=0.5
    cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=map_size)


    Bx_op = ScalarOperator(lambda dset: dset["Br"][:,0] ,  ro.info["unit_velocity"])
    dmap_Bx = slicing.SliceMap(amr, cam_z, Bx_op, z=0.) #,verbose=False)
    map_Bx = np.abs(dmap_Bx.map.T)*scale_mag*microG 
    Bx=np.mean(map_Bx)


    By_op = ScalarOperator(lambda dset: dset["Br"][:,1] ,  ro.info["unit_velocity"])
    dmap_By = slicing.SliceMap(amr, cam_z, By_op, z=0.) #,verbose=False)
    map_By = np.abs(dmap_By.map.T)*scale_mag*microG 
    By=np.mean(map_By)


    Bz_op = ScalarOperator(lambda dset: dset["Br"][:,2] ,  ro.info["unit_velocity"])
    dmap_Bz = slicing.SliceMap(amr, cam_z, Bz_op, z=0.) #,verbose=False)
    map_Bz = np.abs(dmap_Bz.map.T)*scale_mag*microG
    Bz=np.mean(map_Bz)


    return Bx,By,Bz,time 


########################################################################################
########################################################################################
##################################################################
##make images for gravity
##amr and ro pymses stuff
##center 3D array for coordinates center
##radius the size of the region to be plotted (between 0 and 1)
##num output number
##directory : the directory in which the images are being written 
##x_v, y_v and z_v coordinates of the points to be overplotted
##i_im image number
##force test whether the first of the image series exists already. If yes the return
##################################################################
########################################################################################
########################################################################################
def  make_image_phi(amr,ro,center,radius,num,path,x_v=None,y_v=None,z_v=None,i_im=0,force=False,path_out=None):

##    #must be reorganised...
#    Myr=3.15e7*1.e6
#    mu  =  1.4
#    mp =  mu * 1.660531e-24  #n gramme
#    G = 6.7e-8
#    scale_n = 1.
#    scale_d = scale_n*mp  #code density unit is in particle/cc
#    scale_t = 1.0/np.sqrt(G*scale_d)

    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

    lbox_cm = lbox * pc #lbox in cm

    time=ro.info['time'] #time in codeunits 
    time = time * scale_t / Myr

    titre = 't='+ str(time)[0:5] +' (Myr)'


    if( path_out is not None):
        directory = path_out
    else:
        directory = path

    name = directory+'/coldens_z'+'_'+str(i_im)+'_'+format(num,'05')+'.jpeg'
    if (len(glob.glob(name))==1 and not force):
        return 
    
    ## figures along the z-axis
    cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=512)


    ## figures along the y-axis
    cam_y= Camera(center=center,line_of_sight_axis='y',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='z',map_max_size=512)

    ## figures along the x-axis
    cam_x= Camera(center=center,line_of_sight_axis='x',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='z',map_max_size=512)





    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])

    phi_op = ScalarOperator(lambda dset: np.log10(abs(dset["phi"])) ,ro.info["unit_density"])
    dmap_phi = slicing.SliceMap(amr, cam_z, phi_op, z=0.) #,verbose=False)
    map_phi = dmap_phi.map

    P.clf() 
    im = P.imshow(map_phi,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   

    P.xlabel('$x$ (pc)')                                                     
    P.ylabel('$y$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
    cbar.set_label(r'$log(\phi)$')                              

    if(x_v != ''):                                                                              
        P.scatter(x_v,y_v,c='r',marker='x',s=10) 

    P.savefig(directory+'/phi_scat'+'_'+str(i_im)+'_'+format(num,'05')+'.jpeg')



    phi_op = ScalarOperator(lambda dset: np.log10(abs(dset["phi"])) ,ro.info["unit_density"])
    dmap_phi = slicing.SliceMap(amr, cam_x, phi_op, z=0.) #,verbose=False)
    map_phi = dmap_phi.map

    gy_op = ScalarOperator(lambda dset: dset["g"][:,1] ,  ro.info["unit_velocity"])
    dmap_gy = slicing.SliceMap(amr, cam_x, gy_op, z=0.) #,verbose=False)
#    map_gy_red=dmap_gy.map[::20,::20]
#    map_gy_red = map_gy_red.T

    map_gy=dmap_gy.map.T
    map_gy_red = map_gy[::20,::20]


    gz_op = ScalarOperator(lambda dset: dset["g"][:,2],  ro.info["unit_velocity"])
    dmap_gz = slicing.SliceMap(amr, cam_x, gz_op, z=0.) #,verbose=False)
#    map_gz_red=dmap_gz.map[::20,::20]
#    map_gz_red = map_gz_red.T

    map_gz=dmap_gz.map.T
    map_gz_red = map_gz[::20,::20]

    map_g_norm = np.sqrt(map_gy_red**2+map_gz_red**2)

    map_gy_red = map_gy_red / map_g_norm
    map_gz_red = map_gz_red / map_g_norm


    P.clf() 
    im = P.imshow(map_phi,extent=[(-radius+center[1])*lbox,(radius+center[1])*lbox,(-radius+center[2])*lbox,(radius+center[2])*lbox],origin='lower')   


    ny = map_gy_red.shape[0]
    nz = map_gy_red.shape[1]
    vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/ny)*lbox
    vec_z = (np.arange(nz)*2./ny*radius - radius + center[2] + radius/nz)*lbox
    yy,zz = np.meshgrid(vec_y,vec_z)

    P.quiver(yy,zz,map_gy_red,map_gz_red)

    P.title=(titre)
    P.xlabel('$y$ (pc)')          
    P.ylabel('$z$ (pc)')
    cbar = P.colorbar(im)                                                                                                                     
#    cbar.set_label(r'$log(n) \, (cm^{-3})$')                                                                                                   
    cbar.set_label(r'$log(\phi) \, (code units)$')

    P.savefig(directory+'/phi_x'+'_'+str(i_im)+'_'+format(num,'05')+'.jpeg')
    P.savefig(directory+'/phi_x'+'_'+str(i_im)+'_'+format(num,'05')+'.ps')


    dmap_rho = slicing.SliceMap(amr, cam_x, rho_op, z=0.) #,verbose=False)
    map_rho = np.log10(dmap_rho.map)
    map_rho = map_rho.T

#    map_rho_red=map_rho[::20,::20]

    #calculate the tidal forces in 2D
#    tidal_yy = np.roll(map_gy_red,1,axis=0) - np.roll(map_gy_red,-1,axis=0) 
#    tidal_zz = np.roll(map_gz_red,1,axis=1) - np.roll(map_gz_red,-1,axis=1) 

#    tidal_yz = np.roll(map_gy_red,1,axis=1) - np.roll(map_gy_red,-1,axis=1) 
#    tidal_zy = np.roll(map_gz_red,1,axis=0) - np.roll(map_gz_red,-1,axis=0) 


    tidal_yy = np.roll(map_gy,1,axis=0) - np.roll(map_gy,-1,axis=0) 
    tidal_zz = np.roll(map_gz,1,axis=1) - np.roll(map_gz,-1,axis=1) 

    tidal_yz = np.roll(map_gy,1,axis=1) - np.roll(map_gy,-1,axis=1) 
    tidal_zy = np.roll(map_gz,1,axis=0) - np.roll(map_gz,-1,axis=0) 

    ii=0
    nsize=map_rho.size
    wtot = np.zeros(nsize)
    rho  = np.zeros(nsize)
    for j in range(len(map_rho)-2):    
            for k in range(len(map_rho)-2):    

                if(tidal_yy[j+1,k+1]*tidal_yz[j+1,k+1]*tidal_zy[j+1,k+1]*tidal_zz[j+1,k+1] !=0):
                    mat=np.array( [ [ tidal_yy[j+1,k+1] , tidal_yz[j+1,k+1] ] , [ tidal_zy[j+1,k+1] , tidal_zz[j+1,k+1] ] ])
                    w,v=eig(mat)


                    wtot[ii] = np.sum(w)
#                    rho[ii]  = map_rho_red[j,k]
                    rho[ii]  = map_rho[j+1,k+1]
                    ii=ii+1


    mask = [ np.sign(wtot[0:ii]) < 0. ]
    print 'total number of point',ii
    print 'number of point with negative tidal',sum(mask[0])

    P.clf()
    P.plot(rho[0:ii],np.log10(abs(wtot[0:ii]))*np.sign(wtot[0:ii]),'.')
    P.savefig(directory+'/rho_tidal'+'_'+str(i_im)+'_'+format(num,'05')+'.ps')

    return 
##############################################################################
##############################################################################
### get the properties of the particules with just the particules' positions
##############################################################################
##############################################################################
def read_sink_cvs(num,directory,no_cvs=None):

    name = directory+'output_'+str(num).zfill(5)+'/sink_'+str(num).zfill(5)+'.csv'
    print 'read sinks from ', name


    if(no_cvs is None):
        sinks = np.loadtxt(name,delimiter=',',ndmin=2,usecols=(0,1,2,3,4,5,6,7,8)) #,9,10,11,12))
    else:
        sinks = np.loadtxt(name,ndmin=2,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12))

    return sinks 
##############################################################################
##############################################################################
### get normalisation for code quantities
##############################################################################
##############################################################################
def normalisation(lbox):

    mu  =  1.4
    mp =  mu * 1.660531e-24  #n gramme
    G = 6.7e-8
    kbol  =  1.38062e-16   # erg/degre
    pc=3.08e18 #cm in pc
    Ms = 2.e33 #solar mass in g
    Myr =  3.15e13 # time in Myrs
#scale_n converts particle density from user unit into number of
#particle per cm^-3
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
##############################################################################
##############################################################################
### get the mass within sinks as a function of time
##############################################################################
##############################################################################
def draw_mass_sink_cvs(directory,file_exclude=None,path_out=None,no_cvs=None):

    name_search = directory+'output*'

    n_out=len(glob.glob(name_search))

    mass_v = np.zeros(n_out)
    time_v = np.zeros(n_out)

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = directory

    isink=0
    for name  in glob.glob(name_search):

        siz=len(name)
        num=int(name[siz-5:siz])

        skip=False
        if(file_exclude is not None):
            for num_file in file_exclude:
                if(num_file == num):
                    skip=True
            if(skip):
                continue

        print 'read sinks and time for file ', name

#        pdb.set_trace() # debug mode     

        ro=pymses.RamsesOutput(directory,num)
        lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)


        pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

        sinks = read_sink_cvs(num,directory,no_cvs=no_cvs)

        nsink =  sinks.shape [0]

        if (nsink >= 1):

            if(nsink == 1):
                mass = sinks[0,1] * scale_mass / Ms / lbox**3
            else:
                mass = np.sum(sinks[0:nsink,1]) * scale_mass / Ms / lbox**3
#            print mass
            mass_v[isink] = mass
        else:
            mass=0.
            mass_v[isink] = 0.

 
        time=ro.info['time'] #time in codeunits 
        time = time * scale_t / Myr
        time_v[isink] = time

        isink=isink+1

        print 'time ', time , 'mass ', mass


    P.clf()


    P.plot(time_v,mass_v,'o')

    P.xlabel(r'$t \, (Myr)$')
    P.ylabel(r'$M \, (M_s)$')

    P.savefig(directory_out+'sink_mass_time.ps')
    P.savefig(directory_out+'sink_mass_time.jpeg')
    

    return mass_v,time_v

##############################################################################
##############################################################################
### do a series of zoom
### path : the path
### num the output number
### zoom_v an arrays which contains the zoom
### sinks particle can be overplotted and used to place the zoom
##############################################################################
##############################################################################
def make_image_zoom(path,num,zoom_v,path_out=None,sinks=False,force=True,center_zoom=None,make_phi=False,no_cvs=None,col_dens_only=False,tag='',gcomp=True,mag_im=True,vel_red=20,map_size=512,im_pos=0.,conv_sink=1,center_dmax=False,order='<',ps=False,ind_sink=0,all_sink=True):


    ro=pymses.RamsesOutput(path,num,order=order)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    if(not make_phi):
        amr = ro.amr_source(["rho","vel","P","Bl","Br"])
    else:
        amr = ro.amr_source(["rho","vel","P","phi","g"], grav_compat=gcomp)
#        amr = ro.amr_source(["rho","vel","P","phi","g","Bl","Br"], grav_compat=True)
#        amr = ro.amr_source(["rho","vel","P"])

    if(center_dmax):
        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        pos = cells.points
        
        arg_max = np.argmax(cells["rho"])
        center_max=[pos[:,0][arg_max],pos[:,1][arg_max],pos[:,2][arg_max]]

        print 'center_max'
        print center_max


    else:
        center_max=[0.5,0.5,0.5]

    if(sinks):
    #read sink file

        pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(lbox)

        sinks = read_sink_cvs(num,path,no_cvs=no_cvs)


        if(len(sinks) !=0):
            mass_v = sinks[:,1] * scale_mass / Ms / lbox**3
            x_v = sinks[:,2]
            y_v = sinks[:,3]
            z_v = sinks[:,4]

            if(conv_sink == 2):
                x_v = sinks[:,3]
                y_v = sinks[:,4]
                z_v = sinks[:,5]


#position of the most massive sink
#            ind_max = np.argmax(sinks[:,1])

            args = np.argsort(sinks[:,1])[::-1]
            ind_sink = np.min((ind_sink,len(sinks)))
            ind_max = args[ind_sink]
            if(not center_dmax):
                center_max=[x_v[ind_max]/lbox,y_v[ind_max]/lbox,z_v[ind_max]/lbox]
                tag= 'ns'+str(ind_sink)+'_'
            
        else:
            mass_v=None
            x_v=None
            y_v=None
            z_v=None




    else:
        mass_v=None
        x_v=None
        y_v=None
        z_v=None





    i_im=0
    for rad in zoom_v:
        i_im = i_im + 1
        if(rad == 0.5):
           center=[0.5,0.5,0.5]
        else:
            if(center_zoom is None):
                center = center_max
            else:
                center = center_zoom

        print 'center ', center
        print 'rad ', rad 

        if(not make_phi):
            make_image(amr,ro,center,rad,num,path,x_v=x_v,y_v=y_v,z_v=z_v,mass_v=mass_v,i_im=i_im,force=force,path_out=path_out,col_dens_only=col_dens_only,tag=tag,vel_red=vel_red,map_size=map_size,im_pos=im_pos,ps=ps,all_sink=all_sink)
            if (not col_dens_only and mag_im):
                make_image_B(amr,ro,center,rad,num,path,x_v=None,y_v=None,z_v=None,mass_v=None,i_im=i_im,force=force,path_out=path_out,tag=tag,ps=ps)
        else:
            make_image_phi(amr,ro,center,rad,num,path,x_v=x_v,y_v=y_v,z_v=z_v,i_im=i_im,force=force,path_out=path_out)


    return
##############################################################################
##############################################################################
### do images for all output in path
### path_in : the path
### num the output number
### zoom_v an arrays which contains the zoom
### sinks particle can be overplotted and used to place the zoom
##############################################################################
##############################################################################
def make_image_zoom_all(path_in,zoom_v,path_out=None,sinks=False,force=True,center_zoom=None,make_phi=False,no_cvs=None,col_dens_only=False,tag=''):

    name_search= path_in+'output*'
    for name  in glob.glob(name_search):

        siz=len(name)
        num=int(name[siz-5:siz])

        make_image_zoom(path_in,num,zoom_v,path_out=path_out,sinks=sinks,force=force,center_zoom=center_zoom,make_phi=make_phi,no_cvs=no_cvs,col_dens_only=col_dens_only,tag=tag)

    return
##############################################################################
##############################################################################
def look2(path,num):
    ro=pymses.RamsesOutput(path,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)


    amr = ro.amr_source(["rho","passive"])

    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    dx = cells.get_sizes()

    pos = cells.points
    for ip in range(4):

        p = cells["passive"][:,ip]#/cells["rho"]
        imax = np.argmax(p)
        print 'P%i'%ip, p[imax]
        print pos[imax,:]

        mask = p != 0.
        print 'number of cells non equl to 0',np.sum(mask)

    pdb.set_trace() # debug mode     

#    maskx = np.abs(pos[:,0]-0.5-dx/4) < dx / 2.
#    masky = np.abs(pos[:,1]-0.5-dx/4) < dx / 2.
#    maskz = np.abs(pos[:,2]-0.5-dx/4) < dx / 2.
#
#    maskxy=maskx*masky
#    maskxz=maskx*maskz
#    maskyz=masky*maskz
#
#
#    ind_sort = np.argsort(pos[:,2][maskxy])
#
#    P.clf()
#    P.plot(pos[:,2][maskxy][ind_sort],cells["phi"][maskxy][ind_sort])
#    P.savefig(path+'/cut_xy_phi_'+ str(num).zfill(5) +'.ps')
#    P.clf()
#    P.plot(pos[:,2][maskxy][ind_sort],cells["rho"][maskxy][ind_sort])
#    P.savefig(path+'/cut_xy_rho_'+ str(num).zfill(5) +'.ps')
#    P.clf()
#    P.plot(pos[:,2][maskxy][ind_sort],cells["g"][:,2][maskxy][ind_sort])
#    P.savefig(path+'/cut_xy_gz_'+ str(num).zfill(5) +'.ps')
#    P.clf()
#    P.plot(pos[:,2][maskxy][ind_sort],cells["g"][:,1][maskxy][ind_sort])
#    P.savefig(path+'/cut_xy_gy_'+ str(num).zfill(5) +'.ps')
#    P.clf()
#    P.plot(pos[:,2][maskxy][ind_sort],cells["g"][:,0][maskxy][ind_sort])
#    P.savefig(path+'/cut_xy_gx_'+ str(num).zfill(5) +'.ps')
##############################################################################
##############################################################################
def look(path,num):

    ro=pymses.RamsesOutput(path,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)


    amr = ro.amr_source(["rho","vel","P","phi","g"])

    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    dx = cells.get_sizes()

    pos = cells.points


    maskx = np.abs(pos[:,0]-0.5-dx/4) < dx / 2.
    masky = np.abs(pos[:,1]-0.5-dx/4) < dx / 2.
    maskz = np.abs(pos[:,2]-0.5-dx/4) < dx / 2.

    maskxy=maskx*masky
    maskxz=maskx*maskz
    maskyz=masky*maskz


    ind_sort = np.argsort(pos[:,2][maskxy])

    P.clf()
    P.plot(pos[:,2][maskxy][ind_sort],cells["phi"][maskxy][ind_sort])
    P.savefig(path+'/cut_xy_phi_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,2][maskxy][ind_sort],cells["rho"][maskxy][ind_sort])
    P.savefig(path+'/cut_xy_rho_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,2][maskxy][ind_sort],cells["g"][:,2][maskxy][ind_sort])
    P.savefig(path+'/cut_xy_gz_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,2][maskxy][ind_sort],cells["g"][:,1][maskxy][ind_sort])
    P.savefig(path+'/cut_xy_gy_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,2][maskxy][ind_sort],cells["g"][:,0][maskxy][ind_sort])
    P.savefig(path+'/cut_xy_gx_'+ str(num).zfill(5) +'.ps')

    ind_sort = np.argsort(pos[:,1][maskxz])

    P.clf()
    P.plot(pos[:,1][maskxz][ind_sort],cells["phi"][maskxz][ind_sort])
    P.savefig(path+'/cut_xz_phi_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,1][maskxz][ind_sort],cells["rho"][maskxz][ind_sort])
    P.savefig(path+'/cut_xz_rho_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,1][maskxz][ind_sort],cells["g"][:,2][maskxz][ind_sort])
    P.savefig(path+'/cut_xz_gz_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,1][maskxz][ind_sort],cells["g"][:,1][maskxz][ind_sort])
    P.savefig(path+'/cut_xz_gy_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,1][maskxz][ind_sort],cells["g"][:,0][maskxz][ind_sort])
    P.savefig(path+'/cut_xz_gx_'+ str(num).zfill(5) +'.ps')

    ind_sort = np.argsort(pos[:,0][maskyz])

    P.clf()
    P.plot(pos[:,0][maskyz][ind_sort],cells["phi"][maskyz][ind_sort])
    P.savefig(path+'/cut_yz_phi_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,0][maskyz][ind_sort],cells["rho"][maskyz][ind_sort])
    P.savefig(path+'/cut_yz_rho_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,0][maskyz][ind_sort],cells["g"][:,2][maskyz][ind_sort])
    P.savefig(path+'/cut_yz_gz_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,0][maskyz][ind_sort],cells["g"][:,1][maskyz][ind_sort])
    P.savefig(path+'/cut_yz_gy_'+ str(num).zfill(5) +'.ps')
    P.clf()
    P.plot(pos[:,0][maskyz][ind_sort],cells["g"][:,0][maskyz][ind_sort])
    P.savefig(path+'/cut_yz_gx_'+ str(num).zfill(5) +'.ps')

#    pdb.set_trace() # debug mode     

    return
