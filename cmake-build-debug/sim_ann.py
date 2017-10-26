import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
import binascii

lin_sched = [float(x) for x in open('snapshots/sim_ann_bs100_512_lin.txt', 'r').readlines()]
log_sched = [float(x) for x in open('snapshots/sim_ann_bs100_512_log.txt', 'r').readlines()]
exp_sched = [float(x) for x in open('snapshots/sim_ann_bs100_512_exp.txt', 'r').readlines()]

# for i in range(0, len(M)):
#     PMs = [(x[0], x[1]) for x in zip(M, configs)]
#     l = [[x[0] for x in PMs if x[0] == 0], [x[0] for x in PMs if x[0] == 2], [x[0] for x in PMs if x[0] == 4]]

for i in [lin_sched, exp_sched[:-600], log_sched[:80]]:
    plt.plot(i)
    plt.title('Simulated Annealing Logarithmic Scheduling')
    plt.xlabel('Iterations')
    plt.ylabel('Energy')
    plt.show()

