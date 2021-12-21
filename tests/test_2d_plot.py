# -*- coding: utf-8 -*-

import numpy as np
import hyperclip
from matplotlib import pyplot as plt

n = 2
m = 3

np.random.seed(29)
hyperplanes = [hyperclip.Hyperplane().set_by_points(np.random.random((n,n))) for i_m in range(m)]
np.random.seed(None)

X = np.random.random((10**6,n))

id_pos_side = np.ones(X.shape[0])
for hyp in hyperplanes:
    id_pos_side = np.all((id_pos_side, hyp.side(X)), axis=0)

fig, axs = plt.subplots()
axs.set_aspect('equal', 'box')
plt.scatter(X[id_pos_side, 0], X[id_pos_side, 1], s=2, color='gray')

for hyp in hyperplanes:
    sol = hyp.compute_n_solutions()
    x_a, y_a, x_b, y_b = sol.flat
    
    a = (y_b-y_a)/(x_b-x_a)
    b = y_a - x_a * a
    
    y_0 = b
    y_1 = a * 1 + b
    
    plt.plot([0, 1], [y_0, y_1])

plt.xlim([0,1])
plt.ylim([0,1])


hc = hyperclip.Hyperclip().set_hyperplanes(hyperplanes)

vol = hc.volume()
# print('10**6 MonteCarlo : ', id_pos_side.mean(), 'Hyperclip :', vol)

plt.text(0.25,0.2, "10**6 MonteCarlo : "+str(round(id_pos_side.mean(),4)))
plt.text(0.25,0.1, "Hyperclip : "+str(round(vol,4)))