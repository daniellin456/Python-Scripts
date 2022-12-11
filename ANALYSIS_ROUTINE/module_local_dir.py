### This module defines the local directory for different machines ###

import os

def local_dir(dir):
    #ldirs = ['/gpfs/users/', '/drf/projets/alfven-data/','/dsm/anais/storageA/','/drf/projets/capucine/']
    ldirs = ['/data/']
    for ldir in ldirs:
        readdir = ldir+'ynlee/'+dir
        writedir = ldir+'ynlee/'+dir
        #if ldir == ldirs[3]:  writedir = ldir[2]+'ylee/'+dir
        if os.path.exists(readdir): 
            if not os.path.exists(writedir): os.makedirs(writedir) 
            return readdir, writedir
        else:
            if os.path.exists(writedir):
                readdir = writedir
                return readdir, writedir
    return dir, dir



