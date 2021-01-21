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
	               	      'exptime','mask','grating'))
dlab_id = pd.Series(df['datalab']).str[:-4].values
seq_num = pd.Series(df['datalab']).str[-3:].values.astype(np.int)


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
	if (glob.glob(dir_red+dir_wav) == []):
		os.system('mkdir '+dir_red+dir_wav)

	sci_wav = (sci & (df['wav0'] == cen_wav[i]))
	sci_idx = np.where(sci_wav)[0]

	f = open(dir_red+dir_wav+'sci.lis','w')
	g = open(dir_red+dir_wav+'flat.lis','w')
	for j in np.arange(cnt_wav[i]):
		f.write(df['name'][sci_wav][j]+'\n')
		if ((df['obstype'][sci_idx[j]-1] == 'FLAT') & \
			(dlab_id[sci_idx[j]-1] == dlab_id[sci_idx[j]])):
			flat_idx = sci_idx[j]-1
		elif ((df['obstype'][sci_idx[j]+1] == 'FLAT') & \
			  (dlab_id[sci_idx[j]+1] == dlab_id[sci_idx[j]])):
			flat_idx = sci_idx[j]+1
		g.write(df['name'][flat_idx]+'\n')
	f.close()
	g.close()


# ----- STARNDARD & FLAT list ----- #
std = (obj & (df['obsclass'] == 'partnerCal'))
std_idx = np.where(std)[0][0]

if (np.sum(std) != 1):
	raise ValueError("Please check if the number of standard star frame is one!")

f = open(dir_std+'std.lis','w')
f.write(df['name'][std_idx]+'\n')
f.close()

if ((df['obstype'][std_idx-1] == 'FLAT') & \
	(dlab_id[std_idx-1] == dlab_id[std_idx])):
	std_flat_idx = std_idx-1
elif ((df['obstype'][std_idx+1] == 'FLAT') & \
	  (dlab_id[std_idx+1] == dlab_id[std_idx])):
	std_flat_idx = std_idx+1

g = open(dir_std+'std_flat.lis','w')
g.write(df['name'][std_flat_idx]+'\n')
g.close()


# ----- BIAS list ----- #
f = open(dir_bias+'bias.lis','w')
for i in np.arange(np.sum(bias)):
	f.write(df['name'][bias][i]+'\n')
f.close()


# ----- ARC list ----- #
for i in np.arange(np.sum(arc)):
	arc_wav0 = df['wav0'][arc][i]
	dir_wav = f"w{10*arc_wav0:.0f}/"
	f = open(dir_red+dir_wav+'arc.lis','w')
	f.write(df['name'][arc][i]+'\n')
	f.close()

	if (arc_wav0 == df['wav0'][std_idx]):
		g = open(dir_std+'std_arc.lis','w')
		g.write(df['name'][arc][i]+'\n')
		g.close()

