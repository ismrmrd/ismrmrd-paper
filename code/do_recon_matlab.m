
im = recon_dataset('synth.h5');
s = 1.0 / max(max(im{1}));
imwrite(s*rot90(im{1}),'recon_matlab_synth.png')

im = recon_dataset('../data/bruker.h5');
s = 1.0 / max(max(im{1}));
imwrite(s*rot90(im{1}),'recon_matlab_bruker.png')

im = recon_dataset('../data/ge.h5');
s = 1.0 / max(max(im{1}));
imwrite(s*rot90(im{1}),'recon_matlab_ge.png')

im = recon_dataset('../data/philips.h5');
s = 1.0 / max(max(im{1}));
imwrite(s*rot90(im{1}),'recon_matlab_philps.png')

im = recon_dataset('../data/siemens.h5');
s = 1.0 / max(max(im{1}));
imwrite(s*rot90(im{1}),'recon_matlab_siemens.png')
