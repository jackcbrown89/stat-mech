import numpy as np
import matplotlib.pyplot as plt

hist_points = [5, 50, 488]
f = open('distr_prob.txt', 'r').readlines()
runs = np.zeros((len(f), len(f[1].split(',')[1:-1])))

for i in range(0, len(f)-1):
    runs[i] = np.array(f[i].split(',')[1:-1])

print(runs[0][0:25])
hist_steps = [[x[hist_points[0]-1] for x in runs], [x[hist_points[1]-1] for x in runs], [x[hist_points[2]-1] for x in runs]]

for step in hist_steps:
    plt.hist(step)
    plt.show()


# for i in runs[::2000]:
#     plt.plot(i)
#     plt.show()
