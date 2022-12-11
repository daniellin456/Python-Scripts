from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
import os
from module_local_dir import local_dir

def write_pymsesrc(dir):
    iout = search_ro(dir)[-1]
    keys = {'density': 'rho', 'velocit': 'vel', 'B_left_': 'Bl', 'B_right': 'Br', 'thermal': 'P', \
            'radiati': 'e_rad', 'passive': 'passive', 'tempera': 'T','interna': 'e_int','current': 'j'}
    filename = dir+'/output_'+str(iout).zfill(5)+'/hydro_file_descriptor.txt' 
    if not os.path.exists(filename): return
    lines = tuple(open(filename, 'r'))
    vars = [line[14:21] for line in lines[1:]]
    vars_set = []
    [vars_set.append(x) for x in vars if x not in vars_set]
    nvar = len(vars_set)
    inds = [i for i,val in enumerate(vars) if val=='velocit']
    ndim = len(inds)

    f = open('pymsesrc', 'w')
    f.write('{\n')
    f.write('    "Version": 1,\n')
    f.write('    "Multiprocessing max. nproc": 8,\n')
    f.write('    "RAMSES": {\n')
    f.write('        "ndimensions": '+str(ndim)+',\n')
    f.write('        "amr_field_descr": [\n')
    for ivar in range(0,nvar):
        var = vars_set[ivar]
        inds = [i for i,val in enumerate(vars) if val==var]
        inds_var = [int(lines[ind+1][10:12])-1 for ind in inds]
        key = keys[var] 
        if ivar == nvar-1: tail = '},\n'
        else: tail = '},\n'
        if len(inds)==1:
            line = '            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "'+key+'", "ivar": '+str(inds_var[0])+tail
            f.write(line)
        else:
            line = '            {"__type__": "vector_field", "__file_type__": "hydro", "name": "'+key+'", "ivars": '+str(inds_var)+tail
            f.write(line)
    f.write('            {"__type__": "scalar_field", "__file_type__": "grav", "name": "phi", "ivar": 0},\n')
    f.write('            {"__type__": "vector_field", "__file_type__": "grav", "name": "g", "ivars": [1, 2, 3]}\n')
    f.write('        ]\n')
    f.write('    }\n')
    f.write('}\n')
    f.close()





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Find ellipsoids for a series of radii")
    parser.add_argument("read_dir", type=str, help="RAMSES output repository")
    args = parser.parse_args()
    #readdir = '/dsm/anais/storageA/ylee/'+args.read_dir
    #if not os.path.exists(readdir): readdir = '/gpfs/data1/ylee/'+args.read_dir
    readdir = local_dir(args.read_dir)[0]
    write_pymsesrc(readdir) 
