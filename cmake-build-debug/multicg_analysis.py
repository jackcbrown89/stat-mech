import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize

hist_points = [5, 50, 500]
f = open('snapshots/snapshots_mag_nocg_cg1_cg2_b.25.3.6.txt', 'r').readlines()
# runs = [float(x) for x in f][int(len(f)//2):]
stops = [len(f)//3, 2*len(f)//3]
no_cg_runs = [[float(x) for x in f[:stops[0]][::3]], [float(x) for x in f[stops[0]:stops[1]][::3]],
              [float(x) for x in f[stops[1]:][::3]]]
cg1_runs = [[float(x) for x in f[:stops[0]][1::3]], [float(x) for x in f[stops[0]:stops[1]][1::3]],
              [float(x) for x in f[stops[1]:][1::3]]]
cg2_runs = [[float(x) for x in f[:stops[0]][2::3]], [float(x) for x in f[stops[0]:stops[1]][2::3]],
              [float(x) for x in f[stops[1]:][2::3]]]

runs = normalize([no_cg_runs[0][25:], cg1_runs[0][25:], cg2_runs[0][25:],
                  no_cg_runs[1][25:], cg1_runs[1][25:], cg2_runs[1][25:],
                  no_cg_runs[2][25:], cg1_runs[2][25:], cg2_runs[2][25:]])
runs_30 = normalize([no_cg_runs[1][25:], cg1_runs[1][25:], cg2_runs[1][25:]])
runs_60 = normalize([no_cg_runs[2][25:], cg1_runs[2][25:], cg2_runs[2][25:]])

# hist_steps = [[x[hist_points[0]-1] for x in runs], [x[hist_points[1]-1] for x in runs], [x[hist_points[2]-1] for x in runs]]
i=0
for run in runs[6:9]:
    if i == 0:
        plt.plot(run, label='No CG')
    if i == 1:
        plt.plot(run, label='Single CG')
    if i == 2:
        plt.plot(run, label='Double CG')
    i += 1
# plt.errorbar(x=range(0, len(runs[20:]), 5),y=runs[20::5], yerr=[.001 for x in range(0, len(runs[20::5]))], fmt=',')
plt.title('Magnetization vs. Iteration (βJ = 0.6)')
plt.xlabel('Transition Iteration')
plt.ylabel('< m^2 >')
plt.legend()
plt.show()

# for step in hist_steps:
#     plt.hist(step)
#     plt.show()

#
# for i in runs[::2000]:
#     plt.plot(i)
#     plt.show()
