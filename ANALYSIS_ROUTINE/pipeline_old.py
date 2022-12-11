##extraction pipeline
import numpy as np
import module_extract as me
import pymses
import sys
import glob as glob
import os 

#path='/gpfs/data1/phennebe/RB_DENSE_CORE2/'
#num_v = [400,350,300,250,200,150,100,50]



path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [10,20,30,40,50,60,70,80,100,117,128,140,150,160,166,176,188,210,248]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowturb_part2/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,20,24,30,34,40,42,46,50,54,56]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowrot_turb_part2/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,16,18,20,24,28,32,36,40,42,50,60,70,90,100,110]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowrot_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,16,18,20,24,28,30,34,40,50,60,70,78,82,86,90,92]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,16,18,20,24,28,30,34,40,42,50,60,70,80,90,110,110]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,16,20,30,40,50,60,70,80,90,100,110,118,140,160,180,200]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)









path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [1600,1500,1400,1300,1200,1120,1060,1000,960,900,860,808]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_sink/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [294,295,296,297,298,299,300,301,310,320,330,340,286,287,288,289,290,291,292,293]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=True,force=False,center_dmax=True,ps=True)






path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [60,80,100,120,140,150,180,210,230,250,270,290,400,500,600,700]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint2_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [60,80,100,150,170,190,210,230,250,300,400,500,600,700]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_part/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [10,20,30,40,60,80,100,120,140,160]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True,conv_sink=2)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [10,20,30,40,50,100,150,200,250,300]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True,conv_sink=2)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_sink/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [40,50,60,70,80,115,150,200,250,300,350,400]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=True,force=False,center_dmax=True,ps=True,conv_sink=2)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [40,50,60,70,80,115,150,200,250,300,350,400]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=True,force=False,center_dmax=True,ps=True,conv_sink=2)







#me.make_image_zoom_all(path_in,zoom_v,path_out=path_out)




path='/gpfs/data1/phennebe/B335/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [820,800,750,700,650,620,600,580,530,500,450,400,350,300]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hex/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [386]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)









path='/gpfs/data1/phennebe/RB_DENSE_CORE4_qua/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [984]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [300,400,500,600]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [2,4,8,10,20,30]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_ter/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#zoom_v=[0.5/16.,0.5/64.]

num_v = [1000]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True,ps=True)






path='/gpfs/data1/bcommerc/patrick/'

path_out='/gpfs/data1/phennebe/COEUR_BENOIT/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

#num_v = [180,150,100,50]
num_v = [4001,3001,2001]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True,path_out=path_out,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_bis/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

#num_v = [180,150,100,50]
num_v = [50,40,30]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True)



path='/gpfs/data1/phennebe/RB_DENSE_CORE4/'
#path='/gpfs/data1/phennebe/RB_DENSE_CORE3/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE/'


#mass_v, time_v = me.draw_mass_sink_cvs(path)

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

#num_v = [180,150,100,50]
num_v = [70,60,50,40,30]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True)


STOP



path='/gpfs/data1/phennebe/RB_DENSE_CORE_hydro/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE/'


#mass_v, time_v = me.draw_mass_sink_cvs(path)

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

num_v = [20,40,60,100,180]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,mag_im=False)


