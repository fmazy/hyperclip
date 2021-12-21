# -*- coding: utf-8 -*-

import numpy as np
import hyperclip
from matplotlib import pyplot as plt

n = 2
m = 3

hyperplanes = [hyperclip.Hyperplane().set_by_points(np.random.random((n,n))) for i_m in range(m)]

X = np.random.random((10**6,n))

id_pos_side = np.ones(X.shape[0])
for hyp in hyperplanes:
    id_pos_side = np.all((id_pos_side, hyp.side(X)), axis=0)

fig, axs = plt.subplots()
axs.set_aspect('equal', 'box')
plt.scatter(X[id_pos_side, 0], X[id_pos_side, 1], s=2, color='gray')

for hyp in hyperplanes:
    sol = hyp.compute_n_solutions()
    x_a = sol[0,0]
    x_b = sol[1,0]
    y_a = sol[0,1]
    y_b = sol[1,1]
    
    a = (y_b-y_a)/(x_b-x_a)
    b = y_a - x_a * a
    
    y_0 = b
    y_1 = a * 1 + b
    
    plt.plot([0, 1], [y_0, y_1])

plt.xlim([-0.1,1])
plt.ylim([0,1])


hc = hyperclip.Hyperclip().set_hyperplanes(hyperplanes)

vol = hc.volume()
print('10**6 MonteCarlo : ', id_pos_side.mean(), 'Hyperclip :', vol)