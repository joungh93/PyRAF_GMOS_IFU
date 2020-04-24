#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019

@author: jlee
"""


import numpy as np
import glob, os


# ----- Manual setting for standard flat ----- #
std_flat = 'S20190228S0015'    # GCAL flat for the standard star
std_arc_w = '700'    # ARC wav0 for the standard star


# ----- Directories ----- #
dir_bias = 'bias'
dir_std = 'standard'
dir_red = 'redux6'

df = np.genfromtxt('info.txt', dtype=None, encoding='ascii',
	               names=('name','obstype','obsclass','wav0'))

obj = ((df['obstype'] == 'OBJECT') & (df['obsclass'] == 'science'))
cen_wav, cnt_wav = np.unique(df['wav0'][obj], return_counts=True)

os.system('mkdir '+dir_bias)
os.system('mkdir '+dir_std)
for j in np.arange(len(cen_wav)):
	os.system('mkdir '+dir_red+'_'+str(int(cen_wav[j])))


# ----- Writing list files ----- #

# BIAS
bias = (df['obstype'] == 'BIAS')
os.system('mkdir '+dir_bias)
f = open(dir_bias+'/'+'bias.lis','w')
for i in np.arange(np.sum(bias)):
	f.write(df['name'][bias][i]+'\n')
f.close()

# ARC
arc = (df['obstype'] == 'ARC')
cen_wav, cnt_wav = np.unique(df['wav0'][arc], return_counts=True)

for j in np.arange(len(cen_wav)):
	di = dir_red+'_'+str(int(cen_wav[j]))
	f = open(di+'/'+'arc.lis','w')
	arc_w = (df['wav0'][arc] == cen_wav[j])
	for i in np.arange(cnt_wav[j]):
		f.write(df['name'][arc][arc_w][i]+'\n')
	f.close()
	if (cen_wav[j] == float(std_arc_w)):
		os.system('cp -rpv '+di+'/'+'arc.lis '+dir_std+'/'+'std_arc.lis')

# FLAT
flat = ((df['obstype'] == 'FLAT') & (df['obsclass'] == 'partnerCal'))
cen_wav, cnt_wav = np.unique(df['wav0'][flat], return_counts=True)

for j in np.arange(len(cen_wav)):
	di = dir_red+'_'+str(int(cen_wav[j]))
	f = open(di+'/'+'flat.lis','w')
	flat_w = (df['wav0'][flat] == cen_wav[j])
	for i in np.arange(cnt_wav[j]):
		if (df['name'][flat][flat_w][i] != std_flat):
			f.write(df['name'][flat][flat_w][i]+'\n')
	f.close()	

f = open(dir_std+'/'+'std_flat.lis','w')
f.write(std_flat+'\n')
f.close()

# STANDARD
std = ((df['obstype'] == 'OBJECT') & (df['obsclass'] == 'partnerCal'))
f = open(dir_std+'/'+'std.lis','w')
for i in np.arange(np.sum(std)):
	f.write(df['name'][std][i]+'\n')
f.close()

# OBJECT
obj = ((df['obstype'] == 'OBJECT') & (df['obsclass'] == 'science'))
cen_wav, cnt_wav = np.unique(df['wav0'][obj], return_counts=True)

for j in np.arange(len(cen_wav)):
	di = dir_red+'_'+str(int(cen_wav[j]))
	os.system('mkdir '+di)
	f = open(di+'/'+'sci.lis','w')
	obj_w = (df['wav0'][obj] == cen_wav[j])
	for i in np.arange(cnt_wav[j]):
		f.write(df['name'][obj][obj_w][i]+'\n')
	f.close()

