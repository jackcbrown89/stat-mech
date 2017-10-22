import numpy as np
import matplotlib.pyplot as plt

hist_points = [5, 50, 500]
f = open('snapshots_cgs.txt', 'r').read().split(',')[2:-1][::3]
print(len(f))

runs = [np.matrix(x).reshape((9,9)) for x in f]
print(runs[0])

snapshot_steps = [x*len(f)//512 for x in range(1, 512, 16)]

# for step in snapshot_steps:
#     plt.matshow(runs[step])
#     print(step)
#     plt.show()


plt.plot(runs)
print(runs[0].shape)
plt.show()
