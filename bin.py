# apply bin to thermally corrected data saved in .npy format
# takes two directoty address for 'obs.hdr' and 'rdn.hdr' location respectively and band no. as inputs
# apply bin for a single band



import os
import numpy as np
from spectral import *
from tqdm import tqdm
import os

def main(obs_dir, rad_dir, band):
    binned = np.zeros(shape=(120,30,90,2))
    #obs_dir = '../../Data/Data for hapke phase von karman/'
    #rad_dir = '../Initial_thermally_corrected/'
    file_list  = []
    for i in os.listdir(obs_dir):
        filename, extension = os.path.splitext(i)

        if 'rdn' in filename and '.hdr' in extension:
            file_list.append(filename+extension)

    for filename in tqdm(file_list):
        filename = filename.replace('_rdn.hdr', '_rdn.npy')
        filename = filename.replace('_rdn.npy', '_obs.hdr')
        img = envi.open(obs_dir+filename)
        img_obs = img.open_memmap(writable=False)
        img_open = np.load(rad_dir+'corr_i_'+filename)

        for row in range(len(img_open)):

            for col in range(len(img_open[0])):
                phAngle = int(round(img_obs[row,col,4],0))
                emAngle = int(round(img_obs[row,col,3],0))
                inAngle = int(round(img_obs[row,col,1],0))
    
                binned[phAngle, emAngle, inAngle][0] += img_open[row,col,band]
                binned[phAngle, emAngle, inAngle][1] += 1
    return binned
    #np.save(os.path.join('Binned','binned_%d.npy' %band), binned)

