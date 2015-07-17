#!/bin/sh

# Generate synthetic data using C++ program
ismrmrd_generate_cartesian_shepp_logan -o synth.h5

# Copy the scanner data here
cp ../data/*.h5 .






