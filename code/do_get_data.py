import os

# Clean up from any previous runs
os.system("rm -f *.png *.h5")

# Generate synthetic data using C++ program
os.system("ismrmrd_generate_cartesian_shepp_logan -o synth.h5")

# Copy the scanner data here
os.system("cp ../data/*.h5 .")






