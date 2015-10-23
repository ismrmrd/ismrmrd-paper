# oct-setup.sh
# extra flags needed for compiling mex files for octave

# version = 3.6.4
# version = 3.8.1
version = 3.8.2

octlibdir = /opt/local/lib/octave/${version}/
octincdir = /opt/local/include/octave-${version}/octave

# http://stackoverflow.com/questions/7806418/using-setenv-in-makefile
export XTRA_CFLAGS=-std=c99 -UCountAlloc -DMmex -DUse_simd -DUse_thread -O3
