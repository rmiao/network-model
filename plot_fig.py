#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

DATA_PATH = "/tmp/dump/data/"

x = np.loadtxt(DATA_PATH + 't.dat')
# y = np.loadtxt(DATA_PATH + 'u.dat')
y = np.loadtxt(DATA_PATH + 'egress.dat')

# print x,y
# print len(y)

numflow = len(y[0])


print "plot figure"
print "numflow: ", numflow



plt.figure(1)

# Rate figure
plt.subplot(321)
for ii in range(numflow):
    plt.plot(x, y[:,ii])

plt.xlabel('t (Sec)')
plt.ylabel('TX (Gbps)')
plt.ylim(0,60)
plt.title('Egress rate per host')
plt.legend()

# DCQCN Rate figure
y_dcqcn = np.loadtxt(DATA_PATH + 'dcqcn_egress.dat')
plt.subplot(322)
for ii in range(numflow):
    plt.plot(x, y_dcqcn[:,ii])

plt.xlabel('t (Sec)')
plt.ylabel('DCQCN TX (Gbps)')
plt.ylim(0,1000)
plt.title('DCQCN Egress rate per host')
plt.legend()

# Queue Figure
queue = np.loadtxt(DATA_PATH + 'queue.dat')

numqueue = len(queue[0])
plt.subplot(323)
for ii in range(numqueue):
    plt.plot(x, queue[:,ii])

plt.xlabel('t (Sec)')
plt.ylabel('Queue (KB)')
plt.title('Switch queue length per port')
plt.legend()



# CNP Figure
cnp = np.loadtxt(DATA_PATH + 'cnp.dat')
numhost = len(cnp[0])
plt.subplot(324)

print len(cnp[0])
for ii in range(numhost):
    plt.plot(x, cnp[:,ii])
plt.xlabel('t (Sec)')
plt.ylabel('CNP Sent (Kpps)')
plt.title('CNP Sent per Host')
plt.legend()




# Pause count
pause_count = np.loadtxt(DATA_PATH + 'pause_count.dat')
numhost = len(pause_count[0])
plt.subplot(325)

print len(pause_count[0])
for ii in range(numhost):
    plt.plot(x, pause_count[:,ii])
plt.xlabel('t (Sec)')
plt.ylabel('Pause count')
plt.title('Pause count per Host')
plt.legend()


# Pause duration
pause_duration = np.loadtxt(DATA_PATH + 'pause_duration.dat')
numhost = len(pause_duration[0])
plt.subplot(326)

print len(pause_duration[0])
for ii in range(numhost):
    plt.plot(x, pause_duration[:,ii])
plt.xlabel('t (Sec)')
plt.ylabel('Pause duration (us)')
plt.title('Pause duration per Host')
plt.legend()


plt.show()
