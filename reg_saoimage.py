#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 18:13:15 2020

@author: jlee
"""


import numpy as np

def read_region(regfile, regtype='circle'):
	'''
	regfile - name of the region file (dtype: string)
	regtype - shape of the region (default: 'circle')
		'circle' - x_center, y_center, radius
		'ellipse' - x_center, y_center, semi-major axis, semi-minor axis, angle(deg)
		'box' - x_center, y_center, x_size, y_size, angle(deg)
	'''
	f = open(regfile, "r")
	ll = f.readlines()
	f.close()

	if (regtype == 'circle'):
		x0, y0, rad = [], [], []
		for i in np.arange(len(ll)):
			line_split = ll[i][len(regtype)+1:-2].split(",")
			x0.append(np.float32(line_split[0]))
			y0.append(np.float32(line_split[1]))
			rad.append(np.float32(line_split[2]))
		return_list = [x0, y0, rad]
	
	if (regtype == 'ellipse'):
		x0, y0, a, b, angle = [], [], [], [], []
		for i in np.arange(len(ll)):
			line_split = ll[i][len(regtype)+1:-2].split(",")
			x0.append(np.float32(line_split[0]))
			y0.append(np.float32(line_split[1]))
			hrad, vrad, ang0 = np.array(line_split[2:]).astype('float')
			if (hrad >= vrad):
				a.append(hrad)
				b.append(vrad)
				angle.append(ang0)
			else:
				a.append(vrad)
				b.append(hrad)
				if (ang0-90. < 0.):
					angle.append(ang0+270.)
				else:
					angle.append(ang0-90.)
		return_list = [x0, y0, a, b, angle]

	if (regtype == 'box'):
		x0, y0, xsize, ysize, angle = [], [], [], [], []
		for i in np.arange(len(ll)):
			line_split = ll[i][len(regtype)+1:-2].split(",")
			x0.append(np.float32(line_split[0]))
			y0.append(np.float32(line_split[1]))
			xsize.append(np.float32(line_split[2]))
			ysize.append(np.float32(line_split[3]))
			angle.append(np.float32(line_split[4]))
		return_list = [x0, y0, xsize, ysize, angle]

	return return_list
