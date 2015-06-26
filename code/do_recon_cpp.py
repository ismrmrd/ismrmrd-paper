import os

# Reconstructions:
# 1. C++ using the ismrmrd simple 2d recon
os.system("ismrmrd_recon_cartesian_2d synth.h5")
os.system("ismrmrd_recon_cartesian_2d bruker.h5")
os.system("ismrmrd_recon_cartesian_2d ge.h5")
os.system("ismrmrd_recon_cartesian_2d philips.h5")
os.system("ismrmrd_recon_cartesian_2d siemens.h5")





