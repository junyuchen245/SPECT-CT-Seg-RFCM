#! /usr/bin/env python3
from sys import argv, exit, stdout

pids=[]
from subprocess import Popen
from os import waitpid, WNOHANG
from time import sleep

def runcmd(cmd,maxruns=4,waittime=5, debug=False,rundir=None):
	"""Run a limited number of jobs at a time in the background
	runcmd(cmd, maxruns=4, waittime=5, debug=False)

	  Runs the command 'cmd'in the background under the shell. A maximum of
	  'maxruns' jobs at a time can be run. If this function is called and that
	  many jobs are already active then the subroutine waits and checks for
	  one to exist before starting the requested job. It waits 'waittime' seconds
	  between checks for a completed job. If 'debug' is true then lots of
	  information is printed.
	"""
	global pids
	while len(pids) >= maxruns:
		if (debug):
			print("waiting on ",pids)
			stdout.flush()
		running=[]
		for p in pids:
			try:
				n=waitpid(p,WNOHANG)
			except OSError:
				continue
			if len(n) == 2 and n[0] == 0 and n[1] == 0 :
				running.append(p)
		pids=running
		print("running=",pids)
		sleep(waittime)
	if (debug):
		print("run ",cmd)
		stdout.flush()
	
	p=Popen(cmd,cwd=rundir,shell=True)
	pids.append(p.pid)

def waitall(debug=False,waittime=5):
	global pids
	while len(pids) > 0:
		if (debug):
			print("waiting on ",pids)
			stdout.flush()
		running=[]
		for p in pids:
			try:
				n=waitpid(p,WNOHANG)
			except OSError:
				continue
			if len(n) == 2 and n[0] == 0 and n[1] == 0 :
				running.append(p)
		pids=running
		print("still running=",pids)
		if len(pids) > 0:
			sleep(waittime)
	print("all done")
	
	
if __name__ == '__main__':
	for i in range(10):
		runcmd('sleep '+str(3+i*0.1), maxruns=4, waittime=1)
	waitall(waittime=1)
