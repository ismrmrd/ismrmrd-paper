from __future__ import division

import ismrmrd
import ismrmrd.xsd
import math
import numpy as np
import scipy as sp
import scipy.interpolate as interp

import sys
sys.path.insert(1,'./extern')
from ismrmrdtools import coils, transform, show, sense

class noise_struct:
    preWMtx=[]
    bw=0
    def __init__(self,Mtx,bw):
        self.preWMtx = Mtx
        self.bw = bw
        
def reconstruct_noise_scan(filename, noisesetname):
    dset = ismrmrd.Dataset(filename,noisesetname)
    header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
    bw = header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth
    
    # Read the first acquisition to check the number of samples and coils
    nacq_noise = dset.number_of_acquisitions()
    acq = dset.read_acquisition(0)
    nsamp_noise = acq.number_of_samples
    ncoils_noise = acq.active_channels
    sampletime_noise = acq.sample_time_us

    # Accumulate the noise data into a big array
    noisedata = np.zeros([ncoils_noise,nsamp_noise,nacq_noise],dtype='complex')
    for n in range(nacq_noise):
        acq = dset.read_acquisition(n)
        noisedata[:,:,n] = acq.data 
    noisedata = np.reshape(noisedata, [ncoils_noise, nsamp_noise*nacq_noise])

    # Calculate the pre-whitening matrix        
    Mtx = coils.calculate_prewhitening(noisedata, scale_factor=1.0)

    # Close the noise data set
    dset.close()

    return noise_struct(Mtx, bw)

    
def reconstruct_calibration(filename, datasetname, noise=None):
    
    # Handle the imaging data
    dset = ismrmrd.Dataset(filename, datasetname)
    header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
    enc = header.encoding[0]
    # Matrix size
    eNx = enc.encodedSpace.matrixSize.x
    eNy = enc.encodedSpace.matrixSize.y
    rNx = enc.reconSpace.matrixSize.x
    rNy = enc.reconSpace.matrixSize.y
    # Number of Slices
    if enc.encodingLimits.slice != None:
        nslices = enc.encodingLimits.slice.maximum + 1
    else:
        nslices = 1

    # Loop through the acquisitions ignoring the noise scans
    firstscan = 0
    while True:
        acq = dset.read_acquisition(firstscan)
        if acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            firstscan += 1
        else:
            break

    acq = dset.read_acquisition(firstscan)
    ncoils = acq.active_channels
    gre_bw = header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth
    # The calibration data may be have fewer points than the full k
    refNx = acq.number_of_samples
    x0 = (eNx - refNx)// 2
    x1 = eNx - (eNx - refNx)//2
    
    # Reconsparallel imaging calibration scans which are GRE based
    # Initialiaze a storage array for the reference data
    ref_data = np.zeros((nslices, ncoils, eNy, eNx), dtype=np.complex64)
    # Loop             
    # Prewhiten and stuff into the buffer
    scan = firstscan
    while True:
        acq = dset.read_acquisition(scan)
        if acq.isFlagSet(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION):
            slice = acq.idx.slice
            ky = acq.idx.kspace_encode_step_1
            if noise:
                ref_data[slice, :, ky, x0:x1] = coils.apply_prewhitening(acq.data,noise.preWMtx)
            else:
                ref_data[slice, :, ky, x0:x1] = acq.data
            scan += 1
        else:
            break

    # Reconstruct calibration images
    # 2D FFT
    im_ref = transform.transform_kspace_to_image(ref_data, [2,3]) # [slice,coil,x,y]

    # Remove oversampling if needed
    if (eNx != rNx):
        x0 = (eNx - rNx)//2
        x1 = eNx - (eNx - rNx)//2
        im_ref = im_ref[:,:,:,x0:x1]

    # close the data set
    dset.close()

    return im_ref

def reconstruct_gre(filename, datasetname, noise):
    
    # Handle the imaging data
    dset = ismrmrd.Dataset(filename, datasetname)
    header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
    enc = header.encoding[0]
    # Matrix size
    eNx = enc.encodedSpace.matrixSize.x
    eNy = enc.encodedSpace.matrixSize.y
    rNx = enc.reconSpace.matrixSize.x
    rNy = enc.reconSpace.matrixSize.y
    # Number of Slices
    if enc.encodingLimits.slice != None:
        nslices = enc.encodingLimits.slice.maximum + 1
    else:
        nslices = 1

    # Loop through the acquisitions ignoring the noise scans
    firstscan = 0
    while True:
        acq = dset.read_acquisition(firstscan)
        if acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            firstscan += 1
        else:
            break

    acq = dset.read_acquisition(firstscan)
    ncoils = acq.active_channels
    gre_bw = header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth
    
    # Initialiaze a storage array for the data
    data = np.zeros((nslices, ncoils, eNy, eNx), dtype=np.complex64)
    # Loop             
    # Prewhiten and stuff into the buffer
    nacq = dset.number_of_acquisitions()
    for scan in range(firstscan,nacq):
        acq = dset.read_acquisition(scan)
        slice = acq.idx.slice
        ky = acq.idx.kspace_encode_step_1
        data[slice, :, ky, :] = coils.apply_prewhitening(acq.data,noise.preWMtx)

    # Reconstruct calibration images
    # 2D FFT
    im = transform.transform_kspace_to_image(data, [2,3]) # [slice,coil,x,y]

    # Remove oversampling if needed
    if (eNx != rNx):
        x0 = (eNx - rNx)//2
        x1 = eNx - (eNx - rNx)//2
        im = im[:,:,:,x0:x1]

    # close the data set
    dset.close()

    return im

def reconstruct_epi(filename, datasetname, noise, gre):
    
    # Read the epi data
    dset = ismrmrd.Dataset(filename,datasetname)

    ##############################
    # Scan Parameters and Layout #
    ##############################
    header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
    enc = header.encoding[0]
    nkx = enc.encodedSpace.matrixSize.x
    nky = enc.encodedSpace.matrixSize.y
    ncoils = header.acquisitionSystemInformation.receiverChannels
    epi_noise_bw = header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth
    acc_factor = enc.parallelImaging.accelerationFactor.kspace_encoding_step_1
    
    # Number of Slices
    if enc.encodingLimits.slice != None:
        nslices = enc.encodingLimits.slice.maximum + 1
    else:
        nslices = 1

    # Loop through the acquisitions ignoring the noise scans and the
    # parallel imaging calibration scans which are EPI based
    firstscan = 0
    while True:
        acq = dset.read_acquisition(firstscan)
        if acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT) or acq.isFlagSet(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION):
            firstscan += 1
        else:
            break

    #print('First imaging scan at:', firstscan)
    nsamp = acq.number_of_samples
    ncoils = acq.active_channels
    sampletime = acq.sample_time_us

    # The lines are labeled with flags as follows:
    # - Noise or Imaging using ACQ_IS_NOISE_MEASUREMENT
    # - Parallel calibration using ACQ_IS_PARALLEL_CALIBRATION
    # - Forward or Reverse using the ACQ_IS_REVERSE flag
    # - EPI navigator using ACQ_IS_PHASECORR_DATA
    # - First or last in a slice using ACQ_FIRST_IN_SLICE and ACQ_LAST_IN_SLICE
    # - The first navigator in a shot is labeled as first in slice
    # - The first imaging line in a shot is labeled as firt in slice
    # - The last imaging line in a show is labeled as last in slice
    # for n in range(firstscan-1,firstscan+60):
    #   acq = dset.read_acquisition(n)
    #   print(acq.idx.kspace_encode_step_1)
    #   if acq.isFlagSet(ismrmrd.ACQ_FIRST_IN_SLICE):
    #       print('First')
    #   elif acq.isFlagSet(ismrmrd.ACQ_LAST_IN_SLICE):
    #       print('Last')
    #   else:
    #       print('Middle')
    #   if acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE):
    #       print('Reverse')
    #   else:
    #       print('Forward')
    #   if acq.isFlagSet(ismrmrd.ACQ_IS_PHASECORR_DATA):
    #       print('Navigator')

    # The EPI trajectory is described in the XML header
    # for o in enc.trajectoryDescription.userParameterLong:
    #    print(o.name, o.value_)
    #    
    # for o in enc.trajectoryDescription.userParameterDouble:
    #     print(o.name, o.value_)
    tup = tdown = tflat = tdelay = nsamp = nnav = etl = 0
    for o in enc.trajectoryDescription.userParameterLong:
        if o.name == 'rampUpTime':
            tup = o.value_
        if o.name == 'rampDownTime':
            tdown = o.value_
        if o.name == 'flatTopTime':
            tflat = o.value_
        if o.name == 'acqDelayTime':
            tdelay = o.value_
        if o.name == 'numSamples':
            nsamp = o.value_
        if o.name == 'numberOfNavigators':
            nnav = o.value_
        if o.name == 'etl':
            etl = o.value_

    #print(tup, tdown, tflat, tdelay, nsamp, nnav, etl)

    ####################################
    # Calculate the gridding operators #
    ####################################
    nkx = enc.encodedSpace.matrixSize.x
    nx = enc.reconSpace.matrixSize.x
    t = tdelay + sampletime*np.arange(nsamp)
    x = np.arange(nx)/nx-0.5
    up = t<=tup
    flat = (t>tup)*(t<(tup+tflat))
    down = t>=(tup+tflat)

    #Integral of trajectory (Gmax=1.0)
    k = np.zeros(nsamp)
    k[up] = 0.5/tup*t[up]**2
    k[flat] = 0.5*tup + (t[flat] - tup)
    k[down] = 0.5*tup + tflat + 0.5*tdown-0.5/tdown*(tup+tflat+tdown-t[down])**2
    #Scale to match resolution
    k *= nkx/(k[-1]-k[0])
    #Center
    k -= k[nsamp//2]
    kpos = k
    kneg = -1.0*k
    #Corresponding even range
    keven = np.arange(nkx)
    keven -= keven[nkx//2]
    #Forward model
    Qpos = np.zeros([nsamp,nkx])
    Qneg = np.zeros([nsamp,nkx])
    for p in range(nsamp):
        Qpos[p,:] = np.sinc(kpos[p]-keven)
        Qneg[p,:] = np.sinc(kneg[p]-keven)
    #Inverse
    Rpos = np.linalg.pinv(Qpos)
    Rneg = np.linalg.pinv(Qneg)
    #Take transpose because we apply from the right
    Rpos = Rpos.transpose()
    Rneg = Rneg.transpose()

    #################################
    # Calculate the kspace filter   #
    # Hanning filter after gridding #
    #################################
    import scipy.signal
    kfiltx = scipy.signal.hann(nkx)
    kfilty = scipy.signal.hann(nky)
    Rpos = np.dot(Rpos, np.diag(kfiltx))
    Rneg = np.dot(Rneg, np.diag(kfiltx))

    ####################################
    # Calculate SENSE unmixing weights #
    ####################################
    # Some basic checks
    if gre.shape[0] != nslices:
        raise ValueError('Calibration and EPI data have different number of slices')
    if gre.shape[1] != ncoils:
        raise ValueError('Calibration and EPI data have different number of coils')

    # Estimate coil sensitivites from the GRE data
    csm_orig = np.zeros(gre.shape,dtype=np.complex)
    for z in range(nslices):
        (csmtmp, actmp, rhotmp) = coils.calculate_csm_inati_iter(gre[z,:,:,:])
        weight = rhotmp**2 / (rhotmp**2 + .01*np.median(rhotmp.ravel())**2)
        csm_orig[z,:,:,:] = csmtmp*weight
 
    # Deal with difference in resolution
    # Up/down sample the coil sensitivities to the resolution of the EPI
    xcsm = np.arange(gre.shape[3])/gre.shape[3]
    ycsm = np.arange(gre.shape[2])/gre.shape[2]
    xepi = np.arange(nx)/nx
    yepi = np.arange(nky)/nky
    csm = np.zeros([nslices,ncoils,nky,nx],dtype=np.complex)
    for z in range(nslices):
        for c in range(ncoils):
            # interpolate the real part and imaginary part separately
            i_real = interp.RectBivariateSpline(ycsm,xcsm,np.real(csm_orig[z,c,:,:]))
            i_imag = interp.RectBivariateSpline(ycsm,xcsm,np.imag(csm_orig[z,c,:,:]))
            csm[z,c,:,:] = i_real(yepi,xepi) + 1j*i_imag(yepi,xepi)

    # SENSE weights
    unmix = np.zeros(csm.shape,dtype=np.complex)
    for z in range(nslices):
        unmix[z,:,:,:] = sense.calculate_sense_unmixing(acc_factor, csm[z,:,:,:])[0]
    
    ###############
    # Reconstruct #
    ###############
    # Initialize the array for a volume's worth of data
    H = np.zeros([nslices, ncoils, nky, nx],dtype=np.complex)
    # Loop over the slices
    scan = firstscan
    for z in range(nslices):
        #print('Slice %d starts at scan %d.'%(z,scan))
        # Navigator 1
        acq = dset.read_acquisition(scan)
        #print(scan,acq.idx.slice,acq.idx.kspace_encode_step_1,acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE))
        currslice = acq.idx.slice # keep track of the slice number
        data = coils.apply_prewhitening(acq.data,noise.preWMtx)
        if acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE):
            rnav1 = transform.transform_kspace_to_image(np.dot(data, Rneg),dim=[1])
            sgn = -1.0
        else:
            rnav1 = transform.transform_kspace_to_image(np.dot(data, Rpos),dim=[1])
            sgn = 1.0
        scan += 1

        # Navigator 2
        acq = dset.read_acquisition(scan)
        #print(scan,acq.idx.slice,acq.idx.kspace_encode_step_1,acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE))
        data = coils.apply_prewhitening(acq.data,noise.preWMtx)
        if acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE):
            rnav2 = transform.transform_kspace_to_image(np.dot(data, Rneg),dim=[1])
        else:
            rnav2 = transform.transform_kspace_to_image(np.dot(data, Rpos),dim=[1])
        scan += 1

        # Navigator 3
        acq = dset.read_acquisition(scan)
        #print(scan,acq.idx.slice,acq.idx.kspace_encode_step_1,acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE))
        data = coils.apply_prewhitening(acq.data,noise.preWMtx)
        if acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE):
            rnav3 = transform.transform_kspace_to_image(np.dot(data, Rneg),dim=[1])
        else:
            rnav3 = transform.transform_kspace_to_image(np.dot(data, Rpos),dim=[1])
        scan += 1

        # Phase correction
        delta = np.conj(rnav1+rnav3) * rnav2
        fdelta = np.tile(np.mean(delta,axis=0),[ncoils,1])
        corr = np.exp(sgn*1j*np.angle(np.sqrt(fdelta)))

        for j in range(nky):
            acq = dset.read_acquisition(scan)
            slice = acq.idx.slice              
            if slice != currslice:
                # end of this slice
                break

            ky = acq.idx.kspace_encode_step_1
            data = coils.apply_prewhitening(acq.data,noise.preWMtx)
            if acq.isFlagSet(ismrmrd.ACQ_IS_REVERSE):
                rho = transform.transform_kspace_to_image(np.dot(data, Rneg),dim=[1])
                H[slice,:,ky,:] = kfilty[ky]*np.conj(corr)*rho
            else:
                rho = transform.transform_kspace_to_image(np.dot(data, Rpos),dim=[1])
                H[slice,:,ky,:] = kfilty[ky]*corr*rho        
            scan += 1

    # Close the data set
    dset.close()
    
    # Recon in along y
    H = transform.transform_kspace_to_image(H,dim=[2])
    
    # Combine with SENSE weights
    epi_im = np.abs(np.squeeze(np.sum(H*unmix,axis=1)))
    
    return epi_im
