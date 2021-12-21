# -*- coding: utf-8 -*-

import numpy as np
import hyperclip
from matplotlib import pyplot as plt
from time import time

n = 2
m = 3

hyperplanes = [hyperclip.Hyperplane().set_by_points(np.random.random((n,n))) for i_m in range(m)]

st = time()
X = np.random.random((10**6,n))

id_pos_side = np.ones(X.shape[0])
for hyp in hyperplanes:
    id_pos_side = np.all((id_pos_side, hyp.side(X)), axis=0)

mc_time = time()-st

st = time()
hc = hyperclip.Hyperclip().set_hyperplanes(hyperplanes)
vol = hc.volume()

hc_time = time() - st

print('10**6 MonteCarlo :', id_pos_side.mean(), 'Hyperclip :', vol)
print('10**6 MonteCarlo time :', mc_time, 'Hyperclip time :', hc_time)