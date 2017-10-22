import numpy as np
import matplotlib.pyplot as plt

hist_points = [5, 50, 500]
f = open('snapshots_spins_up.txt', 'r').read().split(',')[:-1]
runs = [np.matrix(x).reshape(10, 10) for x in f]
print(runs[0])

snapshot_steps = [0, 10, 50, 100, 250, 500-1]

for step in snapshot_steps:
    plt.matshow(runs[step])
    plt.show()


# for i in runs[::2000]:
#     plt.plot(i)
#     plt.show()
