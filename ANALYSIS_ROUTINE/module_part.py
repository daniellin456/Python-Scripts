##extraction pipeline
import numpy as np
import module_extract as me
import pymses
from pymses.filters import CellsToPoints

import pylab as P
import glob as glob
import pdb
import pickle as pickle




tick_size = 20
fontlabel_size = 15
axislabel_size = 20


params = {'backend': 'wxAgg', 'lines.markersize' : 3, 'axes.labelsize': axislabel_size, 'font.size': fontlabel_size, 'legend.fontsize': fontlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'axes.linewidth' : 3}
P.rcParams.update(params)



#fig_width=8                                                                                                \                                                                                              
#fig_height=4                                                                                               \                                                                                              
#fig_size = [fig_width,fig_height]                                                                                                                                                                         
P.figure(figsize=(8,6))


P.gcf().subplots_adjust(bottom=0.2)
P.gcf().subplots_adjust(left=0.2)



############################################################################################
############################################################################################
def make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=False):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    def getkey2(item):
        siz=len(item)
        key=int(item[siz-5:siz])

        return key

    name_search= path_in+'output*0'

    for name_field  in sorted(glob.glob(name_search),key=getkey2):
        siz=len(name_field)
        num=int(name_field[siz-5:siz])

        ro  = pymses.RamsesOutput(path_in,num)

        make_image_part(ro,num,directory_out,name,subname=subname,nbin=20,ps=True,pdf=False,center_dmax=center_dmax)


    return



############################################################################################
############################################################################################
def make_image_part_num(path_in,name,num_v,subname='',path_out=None,ps=True,pdf=False,center_dmax=False):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    for num  in num_v:
        ro  = pymses.RamsesOutput(path_in,num)

        make_image_part(ro,num,directory_out,name,subname=subname,nbin=20,ps=True,pdf=False,center_dmax=center_dmax)


    return




############################################################################################
############################################################################################
def make_image_part(ro,num,directory_out,name,subname='',nbin=20,ps=True,pdf=False,center_dmax=False):


    if(center_dmax):
        amr = ro.amr_source(["rho"])

        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        pos = cells.points
        
        arg_max = np.argmax(cells["rho"])
        center_max=[pos[:,0][arg_max],pos[:,1][arg_max],pos[:,2][arg_max]]

        print 'center_max'
        print center_max


    else:
        center_max=[0.5,0.5,0.5]


    cx = center_max[0]
    cy = center_max[1]
    cz = center_max[2]


    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    part_src  = ro.particle_source(["id","vel","level","rhop", "tpgp"])

    part = part_src.flatten()

    #ind_sort = np.argsort(part["id"])

    x_v = part.points[:,0]#[ind_sort]
    y_v = part.points[:,1]#[ind_sort]
    z_v = part.points[:,2]#[ind_sort]

#    pdb.set_trace() # debug mode     

    rad_v = np.sqrt((x_v-cx)**2+(y_v-cy)**2+(z_v-cz)**2)

    print 'min rad', np.min(rad_v) , 'max_rad ',np.max(rad_v)

    rad_v = rad_v*lbox

    log_rad_v = np.log10(rad_v)
    log_cmax = np.log10(lbox)
    log_cmin = np.log10(lbox/1e5)
    width=(log_cmax-log_cmin)/nbin

    P.clf()
    hist_part , hist_edges = P.histogram(log_rad_v,bins=nbin,range=(log_cmin,log_cmax))
    P.bar(hist_edges[:-1],np.log10(hist_part),width)

    print 'hist part', hist_part


    if(center_dmax):
        cd = 'cd_'
    else:
        cd=''


#    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(r) \, (pc)$')
    P.ylabel(r'$log(N_{part})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist_npart_'+cd+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist_npart_'+cd+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist_npart_'+cd+str(num).zfill(5)+'.png')


#    pdb.set_trace() # debug mode     

    return



############################################################################################
#this routine looks at the distribution of particles at a reference time and compute some 
#statistics of their spatial distribution through earlier times
############################################################################################
def make_trace_path(path_in,num_ref,path_out=None,ps=True,pdf=False,center_dmax=True,radthr_v=[0.,1.e-4,1.e-3,1.e-2,1.e-1],force=False):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in

    name_stat = 'stat_'+format(num_ref,'05')+'.save'

    if( len(glob.glob(directory_out+name_stat))==1 and not force):
        return


    ro_ref  = pymses.RamsesOutput(path_in,num_ref)

    if(center_dmax):
        amr = ro_ref.amr_source(["rho"])

        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        pos = cells.points
        
        arg_max = np.argmax(cells["rho"])
        center_max=[pos[:,0][arg_max],pos[:,1][arg_max],pos[:,2][arg_max]]

        print 'center_max'
        print center_max


    else:
        center_max=[0.5,0.5,0.5]


    cx = center_max[0]
    cy = center_max[1]
    cz = center_max[2]


    lbox=ro_ref.info['boxlen'] #boxlen in codeunits (=>pc)

    part_src  = ro_ref.particle_source(["id","vel","level","rhop", "tpgp"])

    part_ref = part_src.flatten()

    ind_sort_ref = np.argsort(part_ref["id"])

    x_ref_v = part_ref.points[:,0][ind_sort_ref]
    y_ref_v = part_ref.points[:,1][ind_sort_ref]
    z_ref_v = part_ref.points[:,2][ind_sort_ref]

#    pdb.set_trace() # debug mode     

    rad_ref_v = np.sqrt((x_ref_v-cx)**2+(y_ref_v-cy)**2+(z_ref_v-cz)**2)

    rad_ref_v = rad_ref_v*lbox




    def getkey2(item):
        siz=len(item)
        key=int(item[siz-5:siz])

        return key

    name_search= path_in+'output*'

    #loop over output 
    liste_field=[]
    for name_field  in sorted(glob.glob(name_search),key=getkey2):

        print 'field', name_field

        siz=len(name_field)
        num=int(name_field[siz-5:siz])

        ro  = pymses.RamsesOutput(path_in,num)


        part_src  = ro.particle_source(["id","vel","level","rhop", "tpgp"])
        part = part_src.flatten()

        ind_sort = np.argsort(part["id"])

        x_v = part.points[:,0][ind_sort]
        y_v = part.points[:,1][ind_sort]
        z_v = part.points[:,2][ind_sort]

        #loop over radius 
        #particle located within this radius at the reference time are selected
        dd=[]
        nthr = len(radthr_v)-1
        for ithr in range(nthr):

            mask = rad_ref_v < radthr_v[ithr+1]*lbox 
            maskm = rad_ref_v > radthr_v[ithr]*lbox 
            mask = mask* maskm

            npart = np.sum(mask)

#            pdb.set_trace() # debug mode     

#            id_ref = part_ref["id"][ind_sort][mask]

            id_ref = mask

            min_x=0. ; min_y=0. ; min_z=0. ; min_rad=0.   
            max_x=0. ; max_y=0. ; max_z=0. ; max_rad=0.   
            mean_x=0. ; mean_y=0. ; mean_z=0. ; mean_rad=0.   
            if(npart > 0):

                min_x = np.min(np.abs(x_v[id_ref]-0.5))*lbox
                max_x = np.max(np.abs(x_v[id_ref]-0.5))*lbox
                mean_x = np.mean(np.abs(x_v[id_ref]-0.5))*lbox

                min_y = np.min(np.abs(y_v[id_ref]-0.5))*lbox
                max_y = np.max(np.abs(y_v[id_ref]-0.5))*lbox
                mean_y = np.mean(np.abs(y_v[id_ref]-0.5))*lbox

                min_z = np.min(np.abs(z_v[id_ref]-0.5))*lbox
                max_z = np.max(np.abs(z_v[id_ref]-0.5))*lbox
                mean_z = np.mean(np.abs(z_v[id_ref]-0.5))*lbox


                rad_v = np.sqrt((x_v[id_ref]-0.5)**2+(y_v[id_ref]-0.5)**2+(z_v[id_ref]-0.5)**2)*lbox

                min_rad = np.min(rad_v)
                max_rad = np.max(rad_v)
                mean_rad = np.mean(rad_v)

            data={'min_x':min_x,'min_y':min_y,'min_z':min_z,'min_rad':min_rad}
            data.update({'max_x':max_x,'max_y':max_y,'max_z':max_z,'max_rad':max_rad})
            data.update({'mean_x':mean_x,'mean_y':mean_y,'mean_z':mean_z,'mean_rad':mean_rad})
            data.update({'npart':npart})
        
            dd.append({'rad':radthr_v[ithr+1],'data':data})


        time = ro.info['time']

        pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = me.normalisation(lbox)

        time = time*scale_t / Myr

        liste_field.append({'num':num,'time':time,'stat':dd})


    name_stat = 'stat_'+format(num_ref,'05')+'.save'
    f = open(directory_out+name_stat,'w')
    pickle.dump(liste_field,f)
    f.close()




#    pdb.set_trace() # debug mode     

    return




############################################################################################
#this routine looks at the distribution of particles at a reference time and  
#and plot statistics of their spatial distribution through earlier times
#restoration of calculations done in make_trace_path are done
############################################################################################
def plot_trace_path(path_in,num_ref,path_out=None,ps=True,pdf=False):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    name_stat = 'stat_'+format(num_ref,'05')+'.save'
    f = open(directory_out+name_stat,'r')
    liste_field=pickle.load(f)
    f.close()



    time_v = [lf['time']  for lf in liste_field]

    print 'time' 
    print time_v

    P.clf()
    rad_num = len(liste_field[0]['stat'])
    col_v=['k','r','b','c','g']
    for rn in range(rad_num):

        print

        print 'npart' , liste_field[0]['stat'][rn]['data']['npart'] 

        rad1_mean_v = [lf['stat'][rn]['data']['mean_rad']  for lf in liste_field]
        rad1_min_v = [lf['stat'][rn]['data']['min_rad']  for lf in liste_field]
        rad1_max_v = [lf['stat'][rn]['data']['max_rad']  for lf in liste_field]

        print 'rad_mean' , liste_field[0]['stat'][rn]['rad'] 
        print rad1_mean_v


        P.plot(time_v,np.log10(rad1_mean_v),'-',color=col_v[rn],markersize=2)
        P.plot(time_v,np.log10(rad1_max_v),color=col_v[rn],markersize=1,linestyle='dotted')
        P.plot(time_v,np.log10(rad1_min_v),color=col_v[rn],markersize=1,linestyle='dotted')

    P.xlabel(r'$t \, (Myr)$')
    P.ylabel(r'$log(r) \, (pc)$')

    P.ylim(ymin=np.log10(liste_field[0]['stat'][0]['rad']/10.))


    xw = time_v[0]

    P.text(xw,-4.9,'r='+str(liste_field[0]['stat'][0]['rad']),color=col_v[0])
    P.text(xw,-4.6,'r='+str(liste_field[0]['stat'][1]['rad']),color=col_v[1])
    P.text(xw,-4.3,'r='+str(liste_field[0]['stat'][2]['rad']),color=col_v[2])
    P.text(xw,-4,'r='+str(liste_field[0]['stat'][3]['rad']),color=col_v[3])


    if(ps):
        P.savefig(directory_out+'mean_radius_'+str(num_ref).zfill(5)+'.ps')
    P.savefig(directory_out+'mean_radius_'+str(num_ref).zfill(5)+'.png')





    P.clf()
    for rn in range(rad_num):

        z1_mean_v = [lf['stat'][rn]['data']['mean_z']  for lf in liste_field]
        z1_min_v = [lf['stat'][rn]['data']['min_z']  for lf in liste_field]
        z1_max_v = [lf['stat'][rn]['data']['max_z']  for lf in liste_field]


        P.plot(time_v,np.log10(z1_mean_v),'-',color=col_v[rn],markersize=2)
        P.plot(time_v,np.log10(z1_max_v),color=col_v[rn],markersize=1,linestyle='dotted')
        P.plot(time_v,np.log10(z1_min_v),color=col_v[rn],markersize=1,linestyle='dotted')

    P.xlabel(r'$t \, (Myr)$')
    P.ylabel(r'$log(z) \, (pc)$')

    P.ylim(ymin=np.log10(liste_field[0]['stat'][0]['rad']/10.))



    xw = time_v[0]

    P.text(xw,-4.9,'r='+str(liste_field[0]['stat'][0]['rad']),color=col_v[0])
    P.text(xw,-4.6,'r='+str(liste_field[0]['stat'][1]['rad']),color=col_v[1])
    P.text(xw,-4.3,'r='+str(liste_field[0]['stat'][2]['rad']),color=col_v[2])
    P.text(xw,-4,'r='+str(liste_field[0]['stat'][3]['rad']),color=col_v[3])


    if(ps):
        P.savefig(directory_out+'mean_z_'+str(num_ref).zfill(5)+'.ps')
    P.savefig(directory_out+'mean_z_'+str(num_ref).zfill(5)+'.png')



    P.clf()
    for rn in range(rad_num):

        x1_mean_v = [lf['stat'][rn]['data']['mean_x']  for lf in liste_field]
        x1_min_v = [lf['stat'][rn]['data']['min_x']  for lf in liste_field]
        x1_max_v = [lf['stat'][rn]['data']['max_x']  for lf in liste_field]


        P.plot(time_v,np.log10(x1_mean_v),'-',color=col_v[rn],markersize=2)
        P.plot(time_v,np.log10(x1_max_v),color=col_v[rn],markersize=1,linestyle='dotted')
        P.plot(time_v,np.log10(x1_min_v),color=col_v[rn],markersize=1,linestyle='dotted')

    P.xlabel(r'$t \, (Myr)$')
    P.ylabel(r'$log(x) \, (pc)$')

    P.ylim(ymin=np.log10(liste_field[0]['stat'][0]['rad']/10.))



    xw = time_v[0]

    P.text(xw,-4.9,'r='+str(liste_field[0]['stat'][0]['rad']),color=col_v[0])
    P.text(xw,-4.6,'r='+str(liste_field[0]['stat'][1]['rad']),color=col_v[1])
    P.text(xw,-4.3,'r='+str(liste_field[0]['stat'][2]['rad']),color=col_v[2])
    P.text(xw,-4,'r='+str(liste_field[0]['stat'][3]['rad']),color=col_v[3])


    if(ps):
        P.savefig(directory_out+'mean_x_'+str(num_ref).zfill(5)+'.ps')
    P.savefig(directory_out+'mean_x_'+str(num_ref).zfill(5)+'.png')




#    pdb.set_trace() # debug mode     

    return
