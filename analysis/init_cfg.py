#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019

@author: jlee
"""


import numpy as np
import glob, os


centwave = ['700','680']
nredux = 'redux4'
cube_list, cube_name = [], []
for i in centwave:
	dir_redux = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/'+nredux+'_'+i+'/'
	cube_list += glob.glob(dir_redux+'*_3D.fits')

for i in np.arange(len(cube_list)):
	cube_name.append(cube_list[i].split('/')[-1].split('cstxeqxbrg')[-1].split('_3D.fits')[0])

cube_list = np.array(cube_list)[np.argsort(cube_name)]
# cube_list = sorted(cube_list)

dir_iraf = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/'

cube_spa_off = ['N20190611S0265', 'N20190612S0125', 'N20190612S0128',
                'N20190612S0129', 'N20190613S0229', 'N20190613S0230',
                'N20190613S0233', 'N20190613S0234', 'N20190613S0237']

cube_ref = 'N20190611S0257'
