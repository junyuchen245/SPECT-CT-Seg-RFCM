#! /usr/bin/python3
from os.path import exists
from sys import exit
from os import system
from runcmd3 import runcmd, waitall

numnoise = 50
beta_end = 8
beta_start = 0
gamma_end = 10
gamma_start = 0
gamma_val  = [0,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128, 0.256, 0.512]
beta_val   = [0,0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064,0.0128]
#beta_val = [2.56, 5.12, 10.24, 20.48, 40.96, 81.92, 163.84, 327.68, 655.36]
for beta_i in range(beta_start, beta_end + 1):
    for gamma_i in range(gamma_start, gamma_end + 1):
        for n in range(1, numnoise + 1):
            beta = beta_val[beta_i];
            gamma= gamma_val[gamma_i];
            python_cmd = "eval_seg_3D_cluster_v1.py %1.4f %1.3f %d" % (beta, gamma, n)
            log_filename = 'log_files/run_eval_seg_%1.4f_%1.3f_%d.log' % (beta, gamma, n)
            cmd = "python3 %s >& %s" % (python_cmd, log_filename)
            print(cmd)
            #exit(1)
            runcmd(cmd, waittime = 5, maxruns = 40)
waitall(waittime=1)
