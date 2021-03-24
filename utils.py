# This script checks the mask coordinates obtained with get_mask_coords.py.
#
# This chips are plotted as:
# 1 2 3 4
# 6 5 8 7
#
# Written by Nestor Espinoza, Benjamin Rackham

from matplotlib.pyplot import plot,subplot,imshow,colorbar,show,title,scatter,text
import numpy as np
from astropy.io import fits as pyfits

###################################################
#Classes and Functions

def biassec(x):
    x = x.split(',')
    fnum = ((x[0])[1:]).split(':')
    snum = ((x[1])[:len(x[1])-1]).split(':')
    fnum[0] = int(fnum[0])
    fnum[1] = int(fnum[1])
    snum[0] = int(snum[0])
    snum[1] = int(snum[1])
    return fnum,snum

def zero_oscan(d):
    zrows = len(np.where(d[:,0]==0)[0])
    zcols = len(np.where(d[0,:]==0)[0])
    if(zrows>zcols):
       mode = 'r' 
       length = d.shape[1]
    else:
       mode = 'c' 
       length = d.shape[0]
    cmed = []
    for i in range(length):
        if(mode == 'r'):
           data = d[:,i]
        else:
           data = d[i,:]
        I = np.where(data!=0)[0]
        cmed.append(np.median(data[I]))
    return np.median(np.array(cmed))

def BiasTrim(d,c,h,otype,datasec=None):
    """
    Overscan/Bias correct and Trim an IMACS chip
    """
    # bias has no significant structure, so a single median suffices, I think
    # overscan = [0:49] [2097:2145]
    oxsec,oysec = biassec(h['biassec'])
    if(datasec == None):
       dxsec,dysec = biassec(h['datasec'])
    else:
       dxsec,dysec = biassec(datasec)
    if(otype=='ift'):
       oscan_data = d[(oysec[0]-1):oysec[1],(oxsec[0]-1):oxsec[1]]
       overscan = np.median(oscan_data)
       if(overscan == 0):
          overscan = zero_oscan(oscan_data)
       newdata = d[(dysec[0]-1):dysec[1],(dxsec[0]-1):dxsec[1]] - overscan
    else:
       d = d.transpose()
       oscan_data = d[oxsec[0]-1:oxsec[1],oysec[0]-1:oysec[1]]
       overscan = np.median(oscan_data)
       if(overscan == 0):
          overscan = zero_oscan(oscan_data)
       newdata = d[dxsec[0]-1:dxsec[1],dysec[0]-1:dysec[1]] - overscan
    #overscan = np.median(d[:,2048:2112])
    #newdata = d[:4095,0:2048] - overscan
    if ((c == 'c5') or (c == 'c6') or (c == 'c7') or (c == 'c8')):
       if(otype == 'iff'):
          newdata = newdata[::-1,:]
       else:
          newdata = newdata[::-1,::-1]

    return newdata

class InputCoords:
    """
    A simple class to hold coordinate data.
    """
    def __init__(self, filename, skiplines=0):
        self.fname = filename
        self.obj = np.array([])
        self.chip = np.array([])
        self.x = np.array([])
        self.y = np.array([])

        with open(filename) as f:
            for _ in range(skiplines):
                next(f)
            for line in f:
                splitted = line.split()
                if(len(splitted)==0):
                    break
                self.obj = np.append(self.obj, splitted[0])
                self.chip = np.append(self.chip, splitted[1][-1])
                self.x = np.append(self.x, splitted[2])
                self.y = np.append(self.y, splitted[3])
        self.chip = self.chip.astype(np.int)
        self.x = self.x.astype(np.float)
        self.y = self.y.astype(np.float)


