import random
import os
import cStringIO as StringIO
import gzip
import numpy as np
import astropy.io.fits as aif
import astropy.table as at
import pandas as pd
import multiprocessing as mp

from tempfile import NamedTemporaryFile
import tspot

import sewpy

import unwise_psf
import crowdsource
import psf as cs_psf


#import logging
#logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)
sew = sewpy.SEW(params=["XPSF_IMAGE","YPSF_IMAGE",
                        "XPSF_WORLD","YPSF_WORLD",
                        "MAG_PSF","MAGERR_PSF",
                        "CHI2_PSF",
                        "XWIN_IMAGE","YWIN_IMAGE",
                        "XWIN_WORLD","YWIN_WORLD",
                        "AWIN_IMAGE","ERRAWIN_IMAGE",
                        "BWIN_IMAGE","ERRBWIN_IMAGE",
                        "THETAWIN_IMAGE","ERRTHETAWIN_IMAGE",
                        "X2WIN_IMAGE","ERRX2WIN_IMAGE",
                        "Y2WIN_IMAGE","ERRY2WIN_IMAGE",
                        "XYWIN_IMAGE","ERRXYWIN_IMAGE",
                        "MAG_AUTO","MAGERR_AUTO",
                        "MAG_BEST","MAGERR_BEST",
                        "FLUX_BEST","FLUXERR_BEST",
                        "FLUX_APER(1)","FLUXERR_APER(1)",
                        "BACKGROUND","MU_MAX","FWHM_IMAGE","FLAGS"],
                config={"DETECT_MINAREA":3,
                        "PIXEL_SCALE":2.75,
                        "DETECT_THRESH":0.8,
                        "MAG_ZEROPOINT":22.5,
                        #"FILTER_NAME":"/home/ubuntu/backyard-worlds-scripts-dev/byw/post_hunt/sub.conv"
                },
                sexpath="sex")
sew2 = sewpy.SEW(params=["MAG_PSF","MAGERR_PSF",
                         "CHI2_PSF",
                         "MAG_AUTO","MAGERR_AUTO",
                         "MAG_BEST","MAGERR_BEST",
                         "FLUX_BEST","FLUXERR_BEST",
                         "FLUX_APER(1)","FLUXERR_APER(1)",
                        "BACKGROUND","MU_MAX","FWHM_IMAGE","FLAGS"],
                 config={"DETECT_MINAREA":3,
                         "PIXEL_SCALE":2.75,
                         "DETECT_THRESH":0.8,
                         "MAG_ZEROPOINT":22.5,
                         #"FILTER_NAME":"/home/ubuntu/backyard-worlds-scripts-dev/byw/post_hunt/sub.conv"
                 },
                 sexpath="sex")


def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)
    

def extract_sources(cutout,cutout2=None,psf=None):
    """
    Trying sewpy for getting windowed ellipses
    """
    tf = NamedTemporaryFile()
    tf.write(cutout)
    tf.seek(0)
    if cutout2 is None:
        if psf is not None:
            open(sew._get_psf_filepath(),"wb").write(psf)
        res = sew(tf.name)
    else:
        tf2 = NamedTemporaryFile()
        tf2.write(cutout2)
        tf2.seek(0)
        if psf is not None:
            open(sew2._get_psf_filepath(),"wb").write(psf)
        res = sew2(tf2.name,imgfilepath2=tf.name)
    #os.unlink(res["catfilepath"])
    #os.unlink(res["logfilepath"])
    return res["table"]


def get_coadd(coadd_id,epoch,band):
    """Download coadd from touchspot"""
    # Build the URL
    path_ = "/".join(
        ("unwise/data/timeresolved",
         "e%03d"%int(epoch), # epoch in e### form
         coadd_id[:3], # first 3 digits of RA
         coadd_id, # the tile name itself
         "unwise-%s-w%d-img-m.fits"%(coadd_id,band)))

    print path_
    # Get content from S3
    sio = StringIO.StringIO()
    tspot.bucket.download_fileobj(path_,sio)
    sio.seek(0)
    
    return sio


def get_covmap(coadd_id,epoch,band):
    """Download coverag emap from touchspot"""
    # Build coverage map path
    path_ = "/".join(
        ("unwise/data/timeresolved",
         "e%03d"%int(epoch), # epoch in e### form
         coadd_id[:3], # first 3 digits of RA
         coadd_id, # the tile name itself
         "unwise-%s-w%d-n-m.fits.gz"%(coadd_id,band)))
    print path_
    # Get content from S3
    sio = StringIO.StringIO()
    tspot.bucket.download_fileobj(path_,sio)
    sio.seek(0)
    
    gz = gzip.GzipFile(fileobj=sio,mode="rb").read()
    sio = StringIO.StringIO(gz)
    sio.seek(0)
    
    return sio


def make_fits_image(hdr,coadd):
    """
    Convert header and image data to fits format
    """
    fim = aif.PrimaryHDU(coadd,header=hdr)
    sio = StringIO.StringIO()
    fim.writeto(sio)
    sio.seek(0)
    return sio.getvalue()


def make_fits_psf(psf):
    """
    Convert PSF data array into fits format for source extractor
    """
    print psf.shape
    fim = aif.BinTableHDU.from_columns([aif.Column(name="PSF_MASK",
                                                   format="%dE"%psf.size,
                                                   #dim="(%d, %d, 1)"%psf.shape,
                                                   array=psf)])
    hdr=fim.header
    hdr.set('EXTNAME' , 'PSF_DATA', 'TABLE NAME')
    """
    # TODO: This is from sewpy source. what of this do we need?
    hdr=tbhdu.header
    hdr.set('EXTNAME' , 'PSF_DATA', 'TABLE NAME')
    hdr.set('LOADED' , 36 , 'Number of loaded sources')
    hdr.set('ACCEPTED' , 32 , 'Number of accepted sources')
    hdr.set('CHI2' , 1.4190 , 'Final Chi2')
    hdr.set('POLNAXIS' , 0 , 'Number of context parameters')
    hdr.set('POLNGRP' , 0 , 'Number of context groups')
    hdr.set('PSF_FWHM' , 2.5813 , 'PSF FWHM')
    hdr.set('PSF_SAMP' , 0.5000 , 'Sampling step of the PSF data')
    hdr.set('PSFNAXIS' , 3 , 'Dimensionality of the PSF data')
    hdr.set('PSFAXIS1' , 31 , 'Number of element along this axis')
    hdr.set('PSFAXIS2' , 31 , 'Number of element along this axis')
    hdr.set('PSFAXIS3' , 1 , 'Number of element along this axis')
    """
    hdr.set('PSFNAXIS' , 3 , 'Dimensionality of the PSF data')
    hdr.set('PSFAXIS1' , psf.shape[0] , 'Number of element along this axis')
    hdr.set('PSFAXIS2' , psf.shape[1] , 'Number of element along this axis')
    hdr.set('PSFAXIS3' , 1 , 'Number of element along this axis')
    sio = StringIO.StringIO()
    fim.writeto(sio)
    sio.seek(0)
    return sio.getvalue()


def make_coadd_multi(coadds):
    """
    TODO: CURRENTLY NOT USED

    Old code for coadding coadds and adding synths
    """
    
    # Download TR coadds
    ims = [aif.open(get_coadd(x[0],x[1],x[2]))
           for x in coadds.index]
        
    hdr = ims[0][0].header
    
    ims = [x[0].data for x in ims]
    
    # Add in synthesized planetx
    for i in xrange(len(ims)):
        im = ims[i]
        meta = coadds.iloc[i]
        coadd_id,epoch,band = meta.name
        if "SYNTH_X" in meta and hasattr(meta["SYNTH_X"],"__len__"):
            for j in xrange(len(meta["SYNTH_X"])):
                x,y,flx = meta["SYNTH_X"][j],meta["SYNTH_Y"][j],meta["SYNTH_FLUX"][j]
                gauss = makeGaussian(2048,fwhm=(6.1/2.75),center=(x,y))
                gauss *= flx
                im += gauss

            #open("/tmp/planetx_synthesis_%s_%s_b%d_e%d.fits"
            #     %(coadd_id,
            #       "FORWARD" if meta["FORWARD"] == 1 else "BACKWARD",
            #       band,
            #       epoch),
            #     "wb").write(make_fits_image(hdr,im))
        
        ims[i] = im
        
    
    # Download TR coverage maps
    covmaps = [aif.open(get_covmap(x[0],x[1],x[2]))[0].data
               for x in coadds.index]

    for covmap in covmaps:
        # Set to 0 if <= 2
        covmap[covmap <= 2] = 0

    # Coadd the coadds together w/ weighted average
    # For tiles w/ lots of coadds, this exhausts memory
    # Instead, slice them into rows then re-join
    # TODO: split and array
    return hdr,np.ma.average(ims,weights=covmaps,axis=0)


def add_synths(meta):
    """
    Download coadd and add synths
    """
    coadd_id,epoch,band = meta.name
    
    # Download TR coadds
    im = aif.open(get_coadd(coadd_id,epoch,band))

    # Separate header and image data
    hdr = im[0].header
    im = im[0].data
    
    # Add in synthesized planetx
    if "SYNTH_X" in meta and hasattr(meta["SYNTH_X"],"__len__"):
        for j in xrange(len(meta["SYNTH_X"])):
            x,y,flx = meta["SYNTH_X"][j],meta["SYNTH_Y"][j],meta["SYNTH_FLUX"][j]
            gauss = makeGaussian(2048,fwhm=(6.1/2.75),center=(x,y))
            gauss *= flx
            im += gauss
    
    return hdr,im


def run_extract_sources(w1_coadd,w2_coadd,w1_psf=None,w2_psf=None):
    """
    Run source extractor using sewpy.

    Uses W2 for identifying sources, and W1 for photometry only
    """
    w2_sources = extract_sources(w2_coadd,psf=w2_psf)
    w1_sources = extract_sources(w1_coadd,cutout2=w2_coadd,psf=w1_psf) # TODO: w1 or w2 psf???
    w2_sources["MAG_AUTO_W1"] = w1_sources["MAG_AUTO"]
    w2_sources["MAGERR_AUTO_W1"] = w1_sources["MAGERR_AUTO"]
    w2_sources["MAG_BEST_W1"] = w1_sources["MAG_BEST"]
    w2_sources["MAGERR_BEST_W1"] = w1_sources["MAGERR_BEST"]
    w2_sources["MAG_PSF_W1"] = w1_sources["MAG_PSF"]
    w2_sources["MAGERR_PSF_W1"] = w1_sources["MAGERR_PSF"]
    w2_sources["CHI2_PSF_W1"] = w1_sources["CHI2_PSF"]
    return w2_sources


def work_(coadd_id,coadds,out):
    """
    coadd_id: ID (name/tile name) of the coadd, such as 1194p212
    coadds: Pandas DataFrame describing coadd metadata, such as epochs, bands, dates
    out: directory where to write results
    """
    coadds = coadds.copy()

    w2_meta = coadds.loc[(coadd_id,slice(None),2),:]

    # Plan synthetic planets
    # Pick positions and fluxes
    # TODO: Need to re-add epochal offsets from orbit simulation
    mov_x = []
    mov_y = []
    mov_flx = []
    mjdstart = None
    for i,meta in w2_meta.iterrows():
        if mjdstart is None:
            mjddiff = 0
        else:
            mjddiff = meta["MJDMEAN"]-mjdstart
        xs = []
        ys = []
        flxs = []
        random.seed(0)
        for j in xrange(10): # Number of synths to create
            # Pick pixel positions x, y, and pick flux
            px = random.uniform(0+10,2048-10)
            py = random.uniform(0+10,2048-10)
            flx = (random.uniform(50,400)/6.)
            # TODO: re-add orbit/parallax motion
            xs.append(px)
            ys.append(py)
            flxs.append(flx)
        mov_x.append(xs)
        mov_y.append(ys)
        mov_flx.append(flxs)

    # Store x,y,flx for use when coadding
    coadds.loc[w2_meta.index,"SYNTH_X"] = pd.Series(mov_x,w2_meta.index)
    coadds.loc[w2_meta.index,"SYNTH_Y"] = pd.Series(mov_y,w2_meta.index)
    coadds.loc[w2_meta.index,"SYNTH_FLUX"] = pd.Series(mov_flx,w2_meta.index)

    # Subset metadata for this coadd, epoch 0 only
    w1_meta = coadds.loc[(coadd_id,0,1),:]
    w2_meta = coadds.loc[(coadd_id,0,2),:]

    # Download coadds and add synths
    # (only doing W2, but using same func to grab w1 coadd)
    w1_hdr,w1_im = add_synths(w1_meta)
    w2_hdr,w2_im = add_synths(w2_meta)

    # Convert to fits format
    w1_fits = make_fits_image(w1_hdr,w1_im)
    w2_fits = make_fits_image(w2_hdr,w2_im)

    # Write out resulting files
    # (only doing W2)
    open("%s/%s_w2_synths.fits"%(out,coadd_id),"wb").write(w2_fits)

    # Extract sources & check missed fakes
    def check_f_b(meta,srcs):
        xs = meta["SYNTH_X"][:]
        ys = meta["SYNTH_Y"][:]
        flxs = meta["SYNTH_FLUX"][:]
        sources = []
        miss_xs = []
        miss_ys = []
        miss_flxs = []
        for i in xrange(len(xs)):
            x,y,flx = xs[i],ys[i],flxs[i]

            # Find detection within 3 pixel width square centered on synth
            tmp = srcs[((abs(srcs["XWIN_IMAGE"]-(x+1)) < 1) &
                        (abs(srcs["YWIN_IMAGE"]-(y+1)) < 1))].to_pandas()

            if tmp.shape[0] == 0:
                miss_xs.append(x)
                miss_ys.append(y)
                miss_flxs.append(flx)
                continue
            
            tmp = tmp.iloc[[0]].copy()
            tmp["SYNTH_X"] = x
            tmp["SYNTH_Y"] = y
            tmp["SYNTH_FLUX"] = flx
            sources.append(tmp)
            
        tps = None
        if len(sources) > 0:
            tps = pd.concat(sources,ignore_index=True)
        return tps,pd.DataFrame({"SYNTH_X":miss_xs,
                                 "SYNTH_Y":miss_ys,
                                 "SYNTH_FLUX":miss_flxs,
        })

    # Get PSFs from unwise_psf
    # TODO: Still need to fiure out whether to pass W2 PSF
    # for both W1 and W2
    w1_psf = unwise_psf.get_unwise_psf(1,coadd_id)
    w2_psf = unwise_psf.get_unwise_psf(2,coadd_id)

    # Try crowdsource
    
    # First, make a PSF
    #w2_psf_vp = cs_psf.VariablePixelizedPSF(cs_psf.central_stamp(w2_psf),normalize=-1)
    # TODO: need to decide whether this is the right mechanism for creating a PSF, and
    #       whether to twiddle any knobs (normalize?)
    w2_psf_s = cs_psf.SimplePSF(w2_psf)

    # Then, call fit_im with the image and PSF
    
    # Need to pass in a weight or crowdsource will crash.
    
    # Should recommend considering a weight default = 1
    
    # Also edited code to ignore a None "dq" when writing the flags column
    # need to figure out what that's all about
    
    # These results (x, y, flux, model, psf) have changed since the crowdsource
    # code that described fit_im's usage was written. Now, "x" contains the
    # source table, along with x, y, and so on
    
    # TODO: Figure out correct weight to use. coverage maps?
    x, y, flux, model, psf = crowdsource.fit_im(w2_im,w2_psf_s,weight=1)
    
    # Convert that source table to a dataframe for easy csv-ing
    df = pd.DataFrame(x)
    
    # Write to disk
    df.to_csv("%s/%s_cs_sources.csv"%(out,coadd_id),index=False)

    # Convert PSFs to fits for sextractor
    # Following sewpy's code, but there's some parameters we've left empty. Those may
    # be important.
    # TODO: what should they be?
    w1_psf_fits = make_fits_psf(w1_psf)
    w2_psf_fits = make_fits_psf(w2_psf)

    # Run sextractor
    sources = run_extract_sources(w1_fits,w2_fits,w1_psf=w1_psf_fits,w2_psf=w2_psf_fits)

    # Identify detected vs undetected synths from sextractor results
    # TODO: do this for crowdsource, too.
    tps,fns = check_f_b(w2_meta,sources)

    if tps is not None:
        tps.to_csv("%s/%s_se_tps.csv"%(out,coadd_id),index=False)
    fns.to_csv("%s/%s_se_fns.csv"%(out,coadd_id),index=False)
    
    sources.write("%s/%s_se_sources"%(out,coadd_id),overwrite=True,format="ascii") 


def work(args):
    coadds,out = args
    for coadd_id,coadds_ in coadds.groupby(level=[0]):
        try:
            print coadd_id
            work_(coadd_id,coadds_,out)
        except Exception,e:
            import traceback
            traceback.print_exc()
            print e
            open("%s/%s.failed"%(out,coadd_id),"wb")

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("out",type=str,help="Directory for storing results")
    ap.add_argument("--atlas",default="tr_neo3_index.fits")
    args = ap.parse_args()

    # Prepare output directory
    if not os.path.exists(args.out):
        # Make if it doesn't exist
        os.mkdir(args.out)
    elif not os.path.isdir(args.out):
        # It exists, but it's not a directory!
        raise Exception("'out' parameter is not a directory")

    # Read atlas into pandas
    atlas = aif.open(args.atlas)[1].data
    atlas = pd.DataFrame(atlas.tolist(),columns=atlas.names)
    atlas.set_index(["COADD_ID","EPOCH","BAND"],inplace=True)

    # Build work map
    num_workers = mp.cpu_count()-1
    #split_keys = np.array_split(atlas.index.levels[0],num_workers)
    # Set to only run one coadd for testing purposes
    split_keys = ["1194p212"]

    pool = mp.Pool(num_workers)

    pool.map(work,[(atlas.loc[(keys,slice(None),slice(None)),:],args.out)
                   for keys in split_keys])
    


if __name__ == "__main__": main()

