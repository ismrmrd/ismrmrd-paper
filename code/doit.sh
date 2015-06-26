#!/bin/bash

python do_get_data.py
python do_recon_cpp.py
matlab -r do_recon_matlab
python do_makefigs.py