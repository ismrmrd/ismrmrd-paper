#!/bin/sh

# Copy the scanner data here
curl -L -o ismrmrd_data.zip https://zenodo.org/record/33166/files/ismrmrd_data.zip
unzip ismrmrd_data.zip
find ismrmrd_data -name "*.h5" -exec cp {} . \; -print






