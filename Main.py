'''
Group A26 SVV
'''
import numpy as np

aircraft = "CRJ700" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
Ca = 0.605  # m
la = 2.661  # m
x1 = 0.172  # m
x2 = 1.211  # m
x3 = 2.591  # m
xa = 0.35   # m
ha = 0.205  # m
tsk = 1.1/1000  # m
tsp = 2.8/1000  # m
tst = 1.2/1000  # m
hst = 16./1000   # m
wst = 19./1000   # m
nst = 15  # -
d1 = 0.01154  # m
d3 = 0.01840  # m
theta = m.radians(28)  # rad
P = 97.4*1000  # N