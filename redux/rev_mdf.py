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


# ----- Reading database file ----- #
flat0 = iraf.type(ic.lst_flat, Stdout=1)[0]
apfile = ['aperg'+flat0+'_1', 'aperg'+flat0+'_2']


# ----- Reading MDF file ----- #
mdfdata, hdr = fits.getdata(ic.mdf, extn=1, header=True)

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
		apr_info = copy.deepcopy(dbfile[apr_idx-4:apr_idx-4+28])

		# Revise the MDF database
		if (apr_num >= 50):
			apr_info[1] = apr_info[1].replace(apr_info[1].split()[3], '{0:d}'.format(apr_num+1))
			apr_info[2] = '\ttitle\t{0:.3f}   {1:.3f} '.format(mdfdata['XINST'][apr_num+1-1], mdfdata['YINST'][apr_num+1-1])+mdfdata['BLOCK'][apr_num+1-1]+'\n'
			apr_info[4] = '\taperture\t{0:d}\n'.format(apr_num+1)
			if ((apr_num == 137) | \
				(apr_num == 204) | \
				(apr_num == 449) | \
				(apr_num == 689) | \
				(apr_num == 691)):
				apr_info[1] = apr_info[1].replace(apr_info[1].split()[3], '{0:d}'.format(apr_num+2))
				apr_info[2] = '\ttitle\t{0:.3f}   {1:.3f} '.format(mdfdata['XINST'][apr_num+2-1], mdfdata['YINST'][apr_num+2-1])+mdfdata['BLOCK'][apr_num+2-1]+'\n'
				apr_info[4] = '\taperture\t{0:d}\n'.format(apr_num+2)
			g.writelines(apr_info)
		else:
			g.writelines(apr_info)

		idx_apr_eff.append(int(apr_info[1].split()[3])-1)

g.close()


# ----- slit-2 ----- #
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
		apr_info = copy.deepcopy(dbfile[apr_idx-4:apr_idx-4+28])

		# Revise the MDF database
		g.writelines(apr_info)

		idx_apr_eff.append(int(apr_info[1].split()[3])-1)

g.close()


# ----- Overwriting new MDF file ----- #
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


# Verify the MDF again
flat0 = iraf.type(ic.lst_flat, Stdout=1)[0]
iraf.imdelete('erg@'+ic.lst_flat)
iraf.gfextract('rg'+flat0, fl_inter='yes', line=1400)
