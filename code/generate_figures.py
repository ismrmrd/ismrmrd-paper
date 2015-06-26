import os
import numpy as np
import matplotlib.pyplot as plt
import ismrmrd

# Pull the images out of the output files
fig = plt.figure(1, figsize=(5,5))
plt.axis('off')

dset = ismrmrd.Dataset('synth.h5')
im = dset.read_image('myimage',0)
plt.imshow(np.squeeze(im.data[0,:,-1:0:-1]), cmap='gray')
plt.savefig("recon_cpp_synth.png")

dset = ismrmrd.Dataset('synth.h5')
im = dset.read_image('myimage',0)
plt.figure(1, figsize=(5,5))
plt.imshow(np.squeeze(im.data[0,:,-1:0:-1]), cmap='gray')
plt.axis('off')
plt.savefig("recon_cpp_synth.png")

dset = ismrmrd.Dataset('bruker.h5')
im = dset.read_image('myimage',0)
plt.figure(1, figsize=(5,5))
plt.imshow(np.squeeze(im.data[0,:,-1:0:-1]), cmap='gray')
plt.axis('off')
plt.savefig("recon_cpp_bruker.png")

dset = ismrmrd.Dataset('ge.h5')
im = dset.read_image('myimage',0)
plt.figure(1, figsize=(5,5))
plt.imshow(np.squeeze(im.data[0,:,-1:0:-1]), cmap='gray')
plt.axis('off')
plt.savefig("recon_cpp_ge.png")

dset = ismrmrd.Dataset('philips.h5')
im = dset.read_image('myimage',0)
plt.figure(1, figsize=(5,5))
plt.imshow(np.squeeze(im.data[0,:,-1:0:-1]), cmap='gray')
plt.axis('off')
plt.savefig("recon_cpp_philips.png")

dset = ismrmrd.Dataset('siemens.h5')
im = dset.read_image('myimage',0)
plt.figure(1, figsize=(5,5))
plt.imshow(np.squeeze(im.data[0,:,-1:0:-1]), cmap='gray')
plt.axis('off')
plt.savefig("recon_cpp_siemens.png")

plt.close(fig)

# Move the images to the figures directory
os.system("cp *.png ../figures")

