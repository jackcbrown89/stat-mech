import numpy as np
import matplotlib.pyplot as plt

hist_points = [5, 50, 500]
f = open('snapshots_mag_cg_double.txt', 'r').readlines()
# runs = [float(x) for x in f][int(len(f)//2):]
runs = [float(x) for x in f]
print(f)
# hist_steps = [[x[hist_points[0]-1] for x in runs], [x[hist_points[1]-1] for x in runs], [x[hist_points[2]-1] for x in runs]]
plt.plot(runs[20:])
# plt.errorbar(x=range(0, len(runs[20:]), 5),y=runs[20::5], yerr=[.001 for x in range(0, len(runs[20::5]))], fmt=',')
plt.title('Magnetization vs. Iteration')
plt.xlabel('Transition Iteration')
plt.ylabel('< m^2 >')
plt.show()

# for step in hist_steps:
#     plt.hist(step)
#     plt.show()

#
# for i in runs[::2000]:
#     plt.plot(i)
#     plt.show()
