ISMRMRD Manuscript
-------------------

This repository contains the manuscript and associated example code for the ISMRM Raw Data format.

If you have latex installed, you can generate a PDF of the manuscript with

    $ make

To run the example reconstructions:

    $ cd code
    $ ./doit.sh

This will download the test data and run reconstructions (C++, Matlab, and Python).

For all reconstructions to work, you must have:

1. Compiled and installed ISMRMRD (https://github.com/ismrmrd/ismrmrd)
2. Have added the path (``<ISMRMRD_SOURCE>/matlab``) to the Matlab API to your Matlab installation.
3. Have installed the Python ISMRMRD API (https://github.com/ismrmrd/ismrmrd-python)

On Windows computers or for individual steps (download and reconstructions), please see the ``code/doit.sh`` script for details.
