import numpy as np
import matplotlib.pyplot as plt
import ismrmrd
import h5py

def window_image(im,win_low=1,win_high=95):
    imhist,bins = np.histogram(im.flatten(),100,normed=True)
    im2 = im;
    im2[im2 < bins[win_low]] = bins[win_low]
    im2[im2 > bins[win_high]] = bins[win_high]
    im2 = 255.0*(im2 - bins[win_low])/bins[win_high]
    return im2

fid = h5py.File('bruker.h5','r')
bruk_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
bruk_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
bruk_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
bruk_cpp = window_image(bruk_cpp)
bruk_mat = window_image(bruk_mat)
bruk_py = window_image(bruk_py)

fid = h5py.File('ge.h5','r')
ge_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
ge_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
ge_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
# GE-data has chop so artifact at the edge of the image
ge_cpp[0:1,:]=ge_cpp[2,:]
ge_cpp = window_image(ge_cpp,win_low=20,win_high=50)
ge_mat[0:1,:]=ge_mat[2,:]
ge_mat = window_image(ge_mat,win_low=20,win_high=50)
ge_py[0:1,:]=ge_py[2,:]
ge_py = window_image(ge_py,win_low=20,win_high=50)

fid = h5py.File('philips.h5','r')
phil_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
phil_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
phil_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
phil_cpp = window_image(phil_cpp,win_low=12,win_high=75)
phil_mat = window_image(phil_mat,win_low=12,win_high=75)
phil_py = window_image(phil_py,win_low=12,win_high=75)

fid = h5py.File('siemens.h5','r')
siem_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
siem_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
siem_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
siem_cpp = window_image(siem_cpp,win_low=12,win_high=65)
siem_mat = window_image(siem_mat,win_low=12,win_high=65)
siem_py = window_image(siem_py,win_low=12,win_high=65)

fid = h5py.File('synth.h5','r')
syn_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
syn_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
syn_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
syn_cpp = window_image(syn_cpp)
syn_mat = window_image(syn_mat)
syn_py  = window_image(syn_py)

w, h = 5.,3.
dw, dh = 0.75, 0.4
g=.05
wt = w+dw
ht = h+dh

fig = plt.figure(1,(wt,ht),dpi=600,frameon=False)
ax = fig.add_axes([dw/wt,0,1.-dw/wt,h/ht])
ax.set_axis_off()

ax.imshow(np.hstack(
        (np.vstack((bruk_cpp, bruk_mat, bruk_py)),
         np.vstack((ge_cpp,   ge_mat,   ge_py)),
         np.vstack((phil_cpp, phil_mat, phil_py)),
         np.vstack((siem_cpp, siem_mat, siem_py)),
         np.vstack((syn_cpp,  syn_mat,  syn_py)))),
         cmap='gray')

fig.text(dw/wt-2*g/wt, 1./6.*h/ht,'C++',ha='right',va='center',size=12)
fig.text(dw/wt-2*g/wt, 3./6.*h/ht,'MATLAB',ha='right',va='center',size=12)
fig.text(dw/wt-2*g/wt, 5./6.*h/ht,'Python',ha='right',va='center',size=12)

fig.text(dw/wt+0.1*w/wt,h/ht+g/ht,'Bruker',ha='center',va='bottom',size=12)
fig.text(dw/wt+0.3*w/wt,h/ht+g/ht,'GE',ha='center',va='bottom',size=12)
fig.text(dw/wt+0.5*w/wt,h/ht+g/ht,'Philips',ha='center',va='bottom',size=12)
fig.text(dw/wt+0.7*w/wt,h/ht+g/ht,'Siemens',ha='center',va='bottom',size=12)
fig.text(dw/wt+0.9*w/wt,h/ht+g/ht,'Synthetic',ha='center',va='bottom',size=12)

fig.savefig('figure4_recon_demo.eps',format='eps',dpi=600)
plt.close()


fid = h5py.File('spiral.h5','r')
spiral = np.squeeze(np.array(fid.get('/dataset/matlab')))
spiral = np.reshape(spiral,[spiral.shape[0]*spiral.shape[1], spiral.shape[2]])
fid.close()
spiral = window_image(spiral,win_low=30,win_high=100)

w, h = 3.6,1.0
dw, dh = 0.0, 0.3
g = 0.03
wt = w+dw
ht = h+dh

fig = plt.figure(1,(wt,ht),dpi=600,frameon=False)
ax = fig.add_axes([dw/wt,0,1.-dw/wt,h/ht])
ax.set_axis_off()
ax.imshow(spiral.transpose(),cmap='gray',)
fig.text(dw/wt+0.15*w/wt,h/ht+g/ht,'Nominal',ha='center',va='bottom',size=12)
fig.text(dw/wt+0.50*w/wt,h/ht+g/ht,'Corrected',ha='center',va='bottom',size=12)
fig.text(dw/wt+0.83*w/wt,h/ht+g/ht,'Difference',ha='center',va='bottom',size=12)
fig.savefig('figure5_spiral_demo.eps',format='eps',dpi=600)
plt.close()

fid = h5py.File('epi.h5','r')
epi = np.squeeze(np.array(fid.get('/dataset/python'))[2,:,:])
fid.close()
epi = window_image(epi,win_low=12,win_high=75)
w, h = 3,3
fig = plt.figure(1,(3,3),dpi=600,frameon=False)
ax = fig.add_axes([0,0,1,1])
ax.set_axis_off()
ax.imshow(epi[::-1,:],cmap='gray',)
fig.savefig('figure6_epi_demo.eps',format='eps',dpi=600)

