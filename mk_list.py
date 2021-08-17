#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019
@author: jlee
"""


import numpy as np
import glob, os
import pandas as pd


# ----- Directories ----- #
dir_bias = 'bias/'
dir_std = 'standard/'
dir_red = 'redux/'

for di in [dir_bias, dir_std, dir_red]:
	if (glob.glob(di) == []):
		os.system('mkdir '+di)

current_dir = os.getcwd()


# ----- Reading info.txt ----- #
df = np.genfromtxt('info.txt', dtype=None, encoding='ascii', comments='#',
	               names=('name','obstype','obsclass','wav0','datalab',
	               	      'exptime','mask','grating','airmass','mjd'))
dlab_id = pd.Series(df['datalab']).str[:-4].values
seq_num = pd.Series(df['datalab']).str[-3:].values.astype('int')


# ----- Data classification ----- #
bias = (df['obstype'] == 'BIAS')
arc = (df['obstype'] == 'ARC')
flat = (df['obstype'] == 'FLAT')
obj = (df['obstype'] == 'OBJECT')


# ----- OBJECT & FLAT list ----- #
sci = (obj & (df['obsclass'] == 'science'))
cen_wav, cnt_wav = np.unique(df['wav0'][sci], return_counts=True)
for i in np.arange(len(cen_wav)):
	dir_wav = f"w{10*cen_wav[i]:.0f}/"
	sci_wav = (sci & (df['wav0'] == cen_wav[i]))
	sci_idx = np.where(sci_wav)[0]

	if (glob.glob(dir_red+dir_wav) == []):
		os.system('mkdir '+dir_red+dir_wav)
		for j in np.arange(cnt_wav[i]):
			os.system('mkdir '+dir_red+dir_wav+df['name'][sci_wav][j]+'/')

	for j in np.arange(cnt_wav[i]):
		dir_sci = df['name'][sci_wav][j]+'/'

		f = open(dir_red+dir_wav+dir_sci+'sci.lis','w')
		f.write(df['name'][sci_wav][j]+'\n')
		f.close()

		f = open(dir_red+dir_wav+dir_sci+'flat.lis','w')
		flat_wav = (flat & (df['wav0'] == cen_wav[i]))
		flat_idx = np.abs(df['mjd'][flat_wav]-df['mjd'][sci_wav][j]).argmin()
		f.write(df['name'][flat_wav][flat_idx]+'\n')
		f.close()

		f = open(dir_red+dir_wav+dir_sci+'arc.lis','w')
		arc_wav = (arc & (df['wav0'] == cen_wav[i]))
		arc_idx = np.abs(df['mjd'][arc_wav]-df['mjd'][sci_wav][j]).argmin()
		f.write(df['name'][arc_wav][arc_idx]+'\n')
		f.close()

		# g = open(dir_red+dir_wav+dir_sci+'flat.lis','w')

		# for j in np.arange(cnt_wav[i]):
		# 	f.write(df['name'][sci_wav][j]+'\n')
		# 	if ((df['obstype'][sci_idx[j]-1] == 'FLAT') & \
		# 		(dlab_id[sci_idx[j]-1] == dlab_id[sci_idx[j]])):
		# 		flat_idx = sci_idx[j]-1
		# 	elif ((df['obstype'][sci_idx[j]+1] == 'FLAT') & \
		# 		  (dlab_id[sci_idx[j]+1] == dlab_id[sci_idx[j]])):
		# 		flat_idx = sci_idx[j]+1
		# 	g.write(df['name'][flat_idx]+'\n')
		# f.close()
		# g.close()


# ----- STARNDARD & FLAT list ----- #
std = (obj & (df['obsclass'] == 'partnerCal'))
std_idx = np.where(std)[0][0]

if (np.sum(std) != 1):
	raise ValueError("Please check if the number of standard star frame is one!")

f = open(dir_std+'std.lis','w')
f.write(df['name'][std_idx]+'\n')
f.close()

std_flat_wav = (flat & (df['wav0'] == df['wav0'][std_idx]))
std_flat_idx = np.abs(df['mjd'][std_flat_wav]-df['mjd'][std_idx]).argmin()
f = open(dir_std+'std_flat.lis','w')
f.write(df['name'][std_flat_wav][std_flat_idx]+'\n')
f.close()

std_arc_wav = (arc & (df['wav0'] == df['wav0'][std_idx]))
std_arc_idx = np.abs(df['mjd'][std_arc_wav]-df['mjd'][std_idx]).argmin()
f = open(dir_std+'std_arc.lis','w')
f.write(df['name'][std_arc_wav][std_arc_idx]+'\n')
f.close()


# ----- BIAS list ----- #
f = open(dir_bias+'bias.lis','w')
for i in np.arange(np.sum(bias)):
	f.write(df['name'][bias][i]+'\n')
f.close()


# # ----- ARC list ----- #
# for i in np.arange(np.sum(arc)):
# 	arc_wav0 = df['wav0'][arc][i]
# 	dir_wav = f"w{10*arc_wav0:.0f}/"
# 	f = open(dir_red+dir_wav+'arc.lis','w')
# 	f.write(df['name'][arc][i]+'\n')
# 	f.close()

	# if (arc_wav0 == df['wav0'][std_idx]):
	# 	std_arc_idx = np.abs(df['mjd'][arc] - df['mjd'][std_idx]).argmin()
	# 	g = open(dir_std+'std_arc.lis','w')
	# 	g.write(df['name'][arc][i]+'\n')
	# 	g.close()
