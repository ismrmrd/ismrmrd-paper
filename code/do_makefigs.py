import numpy as np
import matplotlib.pyplot as plt
import ismrmrd
import h5py

fid = h5py.File('bruker.h5','r')
bruk_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
bruk_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
bruk_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
bruk_cpp *= 255/(0.9*np.max(bruk_cpp))
bruk_mat *= 255/(0.9*np.max(bruk_mat))
bruk_py  *= 255/(0.9*np.max(bruk_py))

fid = h5py.File('ge.h5','r')
ge_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
ge_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
ge_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
# GE-data has chop so artifact at the edge of the image
ge_cpp[0:1,:]=ge_cpp[2,:]
ge_cpp *= 255/(0.9*np.max(ge_cpp))
ge_mat[0:1,:]=ge_mat[2,:]
ge_mat *= 255/(0.9*np.max(ge_mat))
ge_py[0:1,:]=ge_py[2,:]
ge_py  *= 255/(0.9*np.max(ge_py))

fid = h5py.File('philips.h5','r')
phil_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
phil_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
phil_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
phil_cpp *= 255/(0.9*np.max(phil_cpp))
phil_mat *= 255/(0.9*np.max(phil_mat))
phil_py  *= 255/(0.9*np.max(phil_py))

fid = h5py.File('siemens.h5','r')
siem_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
siem_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
siem_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
siem_cpp *= 255/(0.9*np.max(siem_cpp))
siem_mat *= 255/(0.9*np.max(siem_mat))
siem_py  *= 255/(0.9*np.max(siem_py))

fid = h5py.File('synth.h5','r')
syn_cpp = np.squeeze(np.array(fid.get('/dataset/cpp/data')))
syn_mat = np.squeeze(np.array(fid.get('/dataset/matlab')))
syn_py  = np.squeeze(np.array(fid.get('/dataset/python')))
fid.close()
syn_cpp = 255/(0.9*np.max(syn_cpp)) * syn_cpp[::-1,:]
syn_mat = 255/(0.9*np.max(syn_mat)) * syn_mat[::-1,:]
syn_py  = 255/(0.9*np.max(syn_py))  * syn_py[::-1,:]

w, h = 5,3
dw, dh = 0.5, 0.5
wt = w+dw
ht = h+dh

fig = plt.figure(1,(wt,ht),dpi=600,frameon=False)
ax = fig.add_axes([dw/wt,0,1,h/ht],frameon=False)
ax.axis('off')
ax.imshow(np.hstack(
        (np.vstack((bruk_cpp, bruk_mat, bruk_py)),
         np.vstack((ge_cpp,   ge_mat,   ge_py)),
         np.vstack((phil_cpp, phil_mat, phil_py)),
         np.vstack((siem_cpp, siem_mat, siem_py)),
         np.vstack((syn_cpp,  syn_mat,  syn_py)))),
         cmap='gray')

fig.text(0.5*dw/wt,.777*h/ht,'C++',ha='center',va='center',size=12)
fig.text(0.5*dw/wt,0.5*h/ht,'MATLAB',ha='center',va='center',size=12)
fig.text(0.5*dw/wt,.133*h/ht,'Python',ha='center',va='center',size=12)

fig.text(dw/wt+0.1*w/wt,(h+dh/2)/ht,'Bruker',ha='center',va='center',size=12)
fig.text(dw/wt+0.3*w/wt,(h+dh/2)/ht,'GE',ha='center',va='center',size=12)
fig.text(dw/wt+0.5*w/wt,(h+dh/2)/ht,'Philips',ha='center',va='center',size=12)
fig.text(dw/wt+0.7*w/wt,(h+dh/2)/ht,'Siemens',ha='center',va='center',size=12)
fig.text(dw/wt+0.9*w/wt,(h+dh/2)/ht,'Synthetic',ha='center',va='center',size=12)

fig.savefig('recon_demo.eps',format='eps',dpi=1200)
