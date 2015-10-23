#!/bin/sh

# Get/synthesize the data
./do_get_data.sh

# CPP recon
./do_recon_cpp.sh

# Matlab recon
matlab -r do_recon_matlab
matlab -r do_spiral_recon_matlab

# Python recon
python do_recon_python.py
python do_epi_recon_python.py

# Extract the PNGs
python do_makefigs.py

# Move the images to the figures directory
mv *.eps ../figures

# Clean up
./do_clean.sh
