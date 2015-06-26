import os
import numpy as np
import matplotlib.pyplot as plt
import ismrmrd
import h5py

# Pull the images out of the output files

rfiles = ['synth', 'bruker', 'ge', 'philips', 'siemens']

for fname in rfiles:

    dset = ismrmrd.Dataset('%s.h5'%fname)
    im = dset.read_image('cpp',0)
    fig = plt.figure(1, figsize=(5,5))    
    plt.axis('off')
    plt.imshow(np.squeeze(im.data[0,:,-1:0:-1]), cmap='gray')
    plt.savefig('recon_cpp_%s.png'%fname)
    plt.close(fig)

    fid = h5py.File('%s.h5'%fname,'r')
    data = fid.get('/dataset/matlab')
    data = np.array(data)
    fig = plt.figure(1, figsize=(5,5))    
    plt.axis('off')
    plt.imshow(data, cmap='gray')
    plt.savefig('recon_matlab_%s.png'%fname)
    plt.close(fig)
    
# Move the images to the figures directory
os.system("cp *.png ../figures")

