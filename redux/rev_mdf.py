#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:25:20 2020

@author: jlee
"""


import numpy as np
import glob, os
import copy
import g0_init_cfg as ic
from astropy.io import fits


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)

iraf.unlearn('gfextract')


# ----- Revising the MDF (copy) ----- #

###########################
########## w6900 ##########
###########################
d = ic.dir_wav[0]    # Run every time per central wavelength
dir_sci = sorted(glob.glob(d+"/*"))
for j in np.arange(len(dir_sci)):
    # Moving each science directory
    name_sci = dir_sci[j].split("/")[-1]
    print("Moving path for "+name_sci+"...")
    os.chdir(current_dir+"/"+dir_sci[j])
    iraf.chdir(current_dir+"/"+dir_sci[j])

    # Reading database file
    flat = np.loadtxt(ic.lst_flat, dtype=str)
    flat0 = flat.item(0)

    # Name of aperture file in the database
    if (ic.nslit == 1):
        apfile = ['aperg'+flat0+'_1']
    if (ic.nslit == 2):
        apfile = ['aperg'+flat0+'_1', 'aperg'+flat0+'_2']

    # Reading MDF file
    mdfdata, hdr = fits.getdata(ic.nmdf, ext=1, header=True)
    idx_apr_eff = []

    # ----- slit-1 ----- #
    f = open(ic.dir_db+apfile[0]+'_old','r')
    dbfile = f.readlines()
    dbfile = np.array(dbfile)
    f.close()
    os.system('rm -rfv '+ic.dir_db+apfile[0])

    g = open(ic.dir_db+apfile[0],'w')
    N_apr = 750
    idx_lines = np.arange(len(dbfile))
    for i in np.arange(N_apr):
        apr_num = i+1
        apr_lines = (dbfile == '\taperture\t{0:d}\n'.format(apr_num))
        if (np.sum(apr_lines) == 1):
            apr_idx = idx_lines[apr_lines][0]
            n_curve = int(dbfile[apr_idx+17].split('\t')[-1].split('\n')[0])
            apr_info = copy.deepcopy(dbfile[apr_idx-4:apr_idx-4+23+n_curve])

            # ----- START ----- #
            g.writelines(apr_info)
            idx_apr_eff.append(int(apr_info[1].split()[3])-1)                  
            # ----- END ----- #     

    g.close()

    # ----- slit-2 ----- #
    if (ic.nslit == 2):
        f = open(ic.dir_db+apfile[1]+'_old','r')
        dbfile = f.readlines()
        dbfile = np.array(dbfile)
        f.close()
        os.system('rm -rfv '+ic.dir_db+apfile[1])

        g = open(ic.dir_db+apfile[1],'w')
        N_apr = 750
        idx_lines = np.arange(len(dbfile))
        for i in np.arange(N_apr):
            apr_num = 750+i+1
            apr_lines = (dbfile == '\taperture\t{0:d}\n'.format(apr_num))
            if (np.sum(apr_lines) == 1):
                apr_idx = idx_lines[apr_lines][0]
                n_curve = int(dbfile[apr_idx+17].split('\t')[-1].split('\n')[0])
                apr_info = copy.deepcopy(dbfile[apr_idx-4:apr_idx-4+23+n_curve])

            # ----- START ----- #
            if (((apr_num >= 795) & (apr_num <= 1150)) | \
                ((apr_num >= 1156) & (apr_num <= 1499)) | \
                (apr_num == 1153)):
                apr_info[1] = apr_info[1].replace(apr_info[1].split()[3], '{0:d}'.format(apr_num+1))
                apr_info[2] = '\ttitle\t{0:.3f}   {1:.3f} '.format(mdfdata['XINST'][apr_num+1-1], mdfdata['YINST'][apr_num+1-1])+mdfdata['BLOCK'][apr_num+1-1]+'\n'
                apr_info[4] = '\taperture\t{0:d}\n'.format(apr_num+1)                       
                g.writelines(apr_info)
            elif ((apr_num == 1151) | \
                  (apr_num == 1154)):
                apr_info[1] = apr_info[1].replace(apr_info[1].split()[3], '{0:d}'.format(apr_num+2))
                apr_info[2] = '\ttitle\t{0:.3f}   {1:.3f} '.format(mdfdata['XINST'][apr_num+2-1], mdfdata['YINST'][apr_num+2-1])+mdfdata['BLOCK'][apr_num+2-1]+'\n'
                apr_info[4] = '\taperture\t{0:d}\n'.format(apr_num+2)                       
                g.writelines(apr_info)            
            else:
                g.writelines(apr_info)
            idx_apr_eff.append(int(apr_info[1].split()[3])-1) 
            # ----- END ----- #

        g.close()

    # Overwriting new MDF file
    newmdfdata = copy.deepcopy(mdfdata)
    bool_apr_eff = np.zeros(len(newmdfdata), dtype=bool)
    bool_apr_eff[idx_apr_eff] = True
    newmdfdata['BEAM'][bool_apr_eff] = 1
    newmdfdata['BEAM'][~bool_apr_eff] = -1

    hdu0 = fits.PrimaryHDU()
    hdu1 = fits.BinTableHDU()
    hdu1.data = newmdfdata
    hdu1.header = hdr
    hdul = fits.HDUList([hdu0, hdu1])
    hdul.writeto(ic.nmdf, overwrite=True)

    # Interative tasks for the first science data for each central wavelength
    if (j == 0):
        dir_db0 = current_dir+"/"+dir_sci[j]+"/"+ic.dir_db
        apfile0 = apfile
        # Verify the MDF again
        iraf.imdelete('erg@'+ic.lst_flat)
        iraf.gfextract('rg'+flat0, fl_inter='yes', line=ic.pk_line, exslits=ic.eslit)
    else:
        for k in np.arange(len(apfile)):
            os.system('cp -rpv '+dir_db0+apfile0[k]+' '+ic.dir_db)

    # Coming back to current path
    os.chdir(current_dir)
    iraf.chdir(current_dir) 


###########################
########## w7000 ##########
###########################
d = ic.dir_wav[1]    # Run every time per central wavelength
dir_sci = sorted(glob.glob(d+"/*"))
for j in np.arange(len(dir_sci)):
    # Moving each science directory
    name_sci = dir_sci[j].split("/")[-1]
    print("Moving path for "+name_sci+"...")
    os.chdir(current_dir+"/"+dir_sci[j])
    iraf.chdir(current_dir+"/"+dir_sci[j])

    # Reading database file
    flat = np.loadtxt(ic.lst_flat, dtype=str)
    flat0 = flat.item(0)

    # Name of aperture file in the database
    if (ic.nslit == 1):
        apfile = ['aperg'+flat0+'_1']
    if (ic.nslit == 2):
        apfile = ['aperg'+flat0+'_1', 'aperg'+flat0+'_2']

    # Reading MDF file
    mdfdata, hdr = fits.getdata(ic.nmdf, ext=1, header=True)
    idx_apr_eff = []

    # ----- slit-1 ----- #
    f = open(ic.dir_db+apfile[0]+'_old','r')
    dbfile = f.readlines()
    dbfile = np.array(dbfile)
    f.close()
    os.system('rm -rfv '+ic.dir_db+apfile[0])

    g = open(ic.dir_db+apfile[0],'w')
    N_apr = 750
    idx_lines = np.arange(len(dbfile))
    for i in np.arange(N_apr):
        apr_num = i+1
        apr_lines = (dbfile == '\taperture\t{0:d}\n'.format(apr_num))
        if (np.sum(apr_lines) == 1):
            apr_idx = idx_lines[apr_lines][0]
            n_curve = int(dbfile[apr_idx+17].split('\t')[-1].split('\n')[0])
            apr_info = copy.deepcopy(dbfile[apr_idx-4:apr_idx-4+23+n_curve])

            # ----- START ----- #
            g.writelines(apr_info)
            idx_apr_eff.append(int(apr_info[1].split()[3])-1)                   
            # ----- END ----- #     

    g.close()

    # ----- slit-2 ----- #
    if (ic.nslit == 2):
        f = open(ic.dir_db+apfile[1]+'_old','r')
        dbfile = f.readlines()
        dbfile = np.array(dbfile)
        f.close()
        os.system('rm -rfv '+ic.dir_db+apfile[1])

        g = open(ic.dir_db+apfile[1],'w')
        N_apr = 750
        idx_lines = np.arange(len(dbfile))
        for i in np.arange(N_apr):
            apr_num = 750+i+1
            apr_lines = (dbfile == '\taperture\t{0:d}\n'.format(apr_num))
            if (np.sum(apr_lines) == 1):
                apr_idx = idx_lines[apr_lines][0]
                n_curve = int(dbfile[apr_idx+17].split('\t')[-1].split('\n')[0])
                apr_info = copy.deepcopy(dbfile[apr_idx-4:apr_idx-4+23+n_curve])

            # ----- START ----- #
            g.writelines(apr_info)
            idx_apr_eff.append(int(apr_info[1].split()[3])-1)  
            # ----- END ----- #

        g.close()

    # Overwriting new MDF file
    newmdfdata = copy.deepcopy(mdfdata)
    bool_apr_eff = np.zeros(len(newmdfdata), dtype=bool)
    bool_apr_eff[idx_apr_eff] = True
    newmdfdata['BEAM'][bool_apr_eff] = 1
    newmdfdata['BEAM'][~bool_apr_eff] = -1

    hdu0 = fits.PrimaryHDU()
    hdu1 = fits.BinTableHDU()
    hdu1.data = newmdfdata
    hdu1.header = hdr
    hdul = fits.HDUList([hdu0, hdu1])
    hdul.writeto(ic.nmdf, overwrite=True)

    # Interative tasks for the first science data for each central wavelength
    if (j == 0):
        dir_db0 = current_dir+"/"+dir_sci[j]+"/"+ic.dir_db
        apfile0 = apfile
        # Verify the MDF again
        iraf.imdelete('erg@'+ic.lst_flat)
        iraf.gfextract('rg'+flat0, fl_inter='yes', line=ic.pk_line, exslits=ic.eslit)
    else:
        for k in np.arange(len(apfile)):
            os.system('cp -rpv '+dir_db0+apfile0[k]+' '+ic.dir_db)

    # Coming back to current path
    os.chdir(current_dir)
    iraf.chdir(current_dir) 

