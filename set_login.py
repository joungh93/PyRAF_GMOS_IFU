#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 8 11:48:02 2023
@author: jlee
"""


import numpy as np
import glob, os


current_dir = os.getcwd()+"/"
# dir_root    = os.path.abspath("../")+"/"


# ----- Reading the initial 'login.cl' ----- #
with open("login.cl", "r") as f:
    lg = f.readlines()


# ----- Username ----- #
os.system("logname > username.txt")
with open("username.txt", "r") as f:
    usr = f.readline()
usr_id = usr.strip()


# ----- Re-writing 'login.cl' file ----- #
idx06_0 = lg[6].index('"')+1
idx06_1 = idx06_0 + lg[6][idx06_0:].index('"')
lg[6]  = lg[6].replace( lg[6][idx06_0:idx06_1], current_dir)

idx07_0 = lg[7].index('"')+1
idx07_1 = idx07_0 + lg[7][idx07_0:].index('"')
lg[7]  = lg[7].replace( lg[7][idx07_0:idx07_1], "/tmp/"+usr_id+"/")

idx10_0 = lg[10].index('"')+1
idx10_1 = idx10_0 + lg[10][idx10_0:].index('"')
lg[10] = lg[10].replace(lg[10][idx10_0:idx10_1], usr_id)

idx74_0 = lg[74].index('=')+1
lg[74] = lg[74].replace(lg[74][idx74_0:], current_dir+"login.cl\n")

with open("login.cl", "w") as f:
    f.writelines(lg)

