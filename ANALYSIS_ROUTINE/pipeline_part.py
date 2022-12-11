##extraction pipeline
import numpy as np
import module_extract as me

import module_part as mpart
import pymses


import pylab as P
import glob as glob
import pdb
import pickle as pickle





path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part2/'
num_ref=249
mpart.make_trace_path(path_in,num_ref,ps=True,pdf=False,force=False) #True)
mpart.plot_trace_path(path_in,num_ref,ps=True,pdf=False)



path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part2/'

name='RB4'

mpart.make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)



path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowturb_part2/'
num_ref=57
mpart.make_trace_path(path_in,num_ref,ps=True,pdf=False,force=False) #True)
mpart.plot_trace_path(path_in,num_ref,ps=True,pdf=False)


path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowturb_part2/'

name='RB4'

mpart.make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)






path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part2/'
num_ref=188
mpart.make_trace_path(path_in,num_ref,ps=True,pdf=False,force=False) #True)
mpart.plot_trace_path(path_in,num_ref,ps=True,pdf=False)


path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part2/'

name='RB4'

mpart.make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)












#
path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part/'
num_ref=800
mpart.make_trace_path(path_in,num_ref,ps=True,pdf=False,force=False) #True)
mpart.plot_trace_path(path_in,num_ref,ps=True,pdf=False)




path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part/'
num_ref=300
mpart.make_trace_path(path_in,num_ref,ps=True,pdf=False,force=False) #True)
mpart.plot_trace_path(path_in,num_ref,ps=True,pdf=False)





path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part/'

name='RB4'

#num_v=[700]
#mpart.make_image_part_num(path_in,name,num_v,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)


mpart.make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)



path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint2_part/'

name='RB4'

mpart.make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)



path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part/'


name='RB4'

mpart.make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)



path_in='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_part/'


name='RB4'

mpart.make_image_part_series(path_in,name,subname='',path_out=None,ps=True,pdf=False,center_dmax=True)

