#
# In order to run this code, you must:
# 1. Install the ismrmrd python api from:
#     https://github.com/ismrmrd/ismrmrd-python
# 2. Download the test data (see do_get_data.sh)
#

import pybits
import h5py


# Get the noise data
noise = pybits.reconstruct_noise_scan('epi.h5','Noise')

# Reconstruct the fully sampled GRE data
gre = pybits.reconstruct_gre('epi.h5','GRE',noise)

# Reconstruct the accelerated EPI data
images = pybits.reconstruct_epi('epi.h5','EPI',noise,gre)
  
# Stick the array into the hdf5 file
fid = h5py.File('epi.h5', 'r+')
fid.create_dataset(name='/dataset/python',data=images)
fid.close()
