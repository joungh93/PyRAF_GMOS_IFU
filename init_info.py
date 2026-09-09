# ----- Imports ----- #
import numpy as np
import glob, os
from astropy.io import fits


current_dir = os.getcwd()
dir_raw = 'raw/'


# ----- Functions ----- #
def clean_string(value):
    """
    Make a header value safe for whitespace-separated info.txt.
    """
    return str(value).strip().replace(" ", "_")
    
    
def get_image_headers(hdul):
    """
    Return image-extension headers only.
    """
    headers = []

    for hdu in hdul[1:]:
        if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
            headers.append(hdu.header)

    return headers
    
    
# ----- Collecting & sorting raw files ----- #
os.chdir(dir_raw)
rawfile = sorted(glob.glob('*.fits'))


# ----- Reading FITS headers ----- #
f = open(current_dir+'/'+'info.txt','w')

f.write(
    "# FILENAME OBJTYPE OBSCLASS CENTWAVE DATALAB EXPTIME "
    "MASKNAME GRATING AIRMASS DATE-OBS MJD-OBS "
    "INSTRUME DETECTOR CCDSUM AMPINTEG ROI GAIN_SIG\n"
)


# ----- Writing the 'info.txt' ----- #
for filename in rawfile:

    with fits.open(filename) as hdul:
    
        h0 = hdul[0].header
        
        image_headers = get_image_headers(hdul)
        if len(image_headers) == 0:
            raise RuntimeError(f"No image extensions found: {filename}")

        h1 = image_headers[0]
        
        ### Basic information
        objtype = clean_string(h0.get("OBSTYPE", "NA"))
        objclass = clean_string(h0.get("OBSCLASS", "NA"))
        instrument = clean_string(h0.get("INSTRUME", "NA"))
        
        if instrument.startswith("GMOS"):
            centwave = str(h0.get("CENTWAVE", -999.0))
            grating  = clean_string(h0.get("GRATING", "NA"))
        elif instrument == "F2":
            centwave = str(h0.get("WAVELENG", -999.0))
            grating  = clean_string(h0.get("GRISM", "NA"))
        else:
            centwave = "-999"
            grating  = "NA"
        
        datalabel = clean_string(h0.get("DATALAB", "NA"))
        exptime = f"{float(h0.get('EXPTIME', 0.0)):.1f}"
        mask = clean_string(h0.get("MASKNAME", "NA"))

        try:
            airmass = f"{float(h0.get('AIRMASS', np.nan)):.4f}"
        except Exception:
            airmass = "nan"
        
        date = clean_string(h0.get("DATE-OBS", "NA"))
        mjd  = h0.get("MJD-OBS", h1.get("MJD-OBS", np.nan))
        mjd  = f"{float(mjd):.8f}"


        ### Detector configuration
        detector = clean_string(
            h0.get("DETECTOR", h0.get("DETTYPE", "NA"))
        )

        ccdsum = clean_string(
            h1.get("CCDSUM", h0.get("CCDSUM", "NA"))
        )

        ampinteg = clean_string(
            h0.get("AMPINTEG", h1.get("AMPINTEG", "NA"))
        )

        # ROI signature:
        roi_list = []
        for hdr in image_headers:
            sec = hdr.get("DETSEC", hdr.get("CCDSEC", "NA"))
            roi_list.append(
                str(sec).replace(" ", "")
            )
        roi_sig = "|".join(roi_list)
        
        # Gain signature
        gain_list = []
        for hdr in image_headers:
            gain = hdr.get("GAIN", np.nan)
            try:
                gain_list.append(f"{float(gain):.2f}")
            except Exception:
                gain_list.append("NA")
        gain_sig = ",".join(gain_list)        

        ### Writing the output files
        name = os.path.splitext(filename)[0]
        f.write(
            f"{name} "
            f"{objtype} "
            f"{objclass} "
            f"{centwave} "
            f"{datalabel} "
            f"{exptime} "
            f"{mask} "
            f"{grating} "
            f"{airmass} "
            f"{date} "
            f"{mjd} "
            f"{instrument} "
            f"{detector} "
            f"{ccdsum} "
            f"{ampinteg} "
            f"{roi_sig} "
            f"{gain_sig}\n"
        )

f.close()

os.chdir(current_dir)
