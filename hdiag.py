import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import os
import myfun as mf

plt.rcdefaults()
plt.rc('lines', linewidth=1.5)
plt.rc('axes', labelsize=26)
plt.rc('xtick', labelsize=26)
plt.rc('ytick', labelsize=26)
plt.rc("legend", fontsize=15)  

##### Cut parameters
x1_min, x1_max, chisqred_max, pmax_g, pmax_r, nmin_g, nmin_r = mf.hdiag_params()

### Histogram Bin size
dz = 0.01
dphase = 1
dx1 = 0.5
dc = 0.05

### Load most recent fit 
files = []
path = "runs/"
for which in os.listdir(path):
	if os.path.isfile(os.path.join(path,which)) and "fitall" in which:
		files.append(str(which))
filename = files[-1]

a = genfromtxt(path+filename,dtype=str)
SN = a[:,0]

a = genfromtxt(path+filename)
z = a[:,1]
t0 = a[:,2]
x0 = a[:,3]
x1 = a[:,4]
c = a[:,5]
ebv_mw = a[:,6]
rv_mw = a[:,7]
mu = a[:,8]
sig_mu = a[:,9]
firstg = a[:,10]
firstr = a[:,11]
ng = a[:,12]
nr = a[:,13]
chisq = a[:,14]
ndof = a[:,15]

a = open(path+filename)
for i, line in enumerate(a):
	if i==0:
		Nmarsh = float(line.split(' ')[1])

chisq_red = np.array(chisq) / ndof

print("%s are missing g-band photometry" %len(np.array(firstg)[np.where(np.array(firstg)==1e6)]))
print("%s are missing r-band photometry" %len(np.array(firstr)[np.where(np.array(firstr)==1e6)]))

######################################################
####################### HISTOGRAMS ###################
######################################################

fig1=plt.figure(figsize=(12,8))   

ax1=fig1.add_subplot(2,2,1)
ax1.set_xlabel("Redshift")
ax2=fig1.add_subplot(2,2,2)
ax2.set_xlabel("Phase")
ax3=fig1.add_subplot(2,2,3)
ax3.set_xlabel("x1")
ax4=fig1.add_subplot(2,2,4)
ax4.set_xlabel("c")

###### Redhsift
minz = 0
maxz = max(z)
nz = int(maxz/dz)+1
ax1.hist(z,bins=arange(-dz/2.,dz*(nz+1/2.),dz),label="%d of %d SNe"%(len(z),Nmarsh))

##### Apply cuts
w=list(zip(SN,x1,c,firstg,firstr,ng,nr,z,mu,sig_mu,chisq_red))
print(type(chisqred_max), chisqred_max)
chisqred_max = np.float(chisqred_max)
ww=list(filter(lambda r:r[1]>=x1_min and r[1]<=x1_max 
					    and r[3]<=pmax_g and r[4]<=pmax_r 
					    and r[5]>=nmin_g and r[6]>=nmin_r
					    and r[10]<=chisqred_max,w))
SN,x1,c,firstg,firstr,ng,nr,z,mu,sig_mu,chisqnew = list(zip(*ww))

firstep = min(firstg,firstr)

#print np.array(SN)[np.where(np.array(firstep)<-14)]

###### First epoch
absmax = max(-min(firstep),max(firstep))
nep = int(absmax/dphase)+1
ax2.hist(firstep,bins=arange(-dphase*(nep+1/2.),dphase*(nep+1/2.),dphase),label="%d of %d SNe"%(len(firstep),Nmarsh))

###### x1
absmax = max(-min(x1),max(x1))
nx1 = int(absmax/dx1)+1
ax3.hist(x1,bins=arange(-dx1*(nx1+1/2.),dx1*(nx1+1/2.),dx1),label="%d of %d SNe"%(len(x1),Nmarsh))

###### c
absmax = max(-min(c),max(c))
nc = int(absmax/dc)+1
ax4.hist(c,bins=arange(-dc*(nc+1/2.),dc*(nc+1/2.),dc),label="%d of %d SNe"%(len(c),Nmarsh))

###### PLOT
ax1.set_xlim(0,max(z)+dz/2.)
ax2.set_xlim(-22,max(firstep))
ax3.set_xlim(x1_min,x1_max)
ax4.set_xlim(-max(-min(c),max(c)),max(-min(c),max(c)))

subplots_adjust(hspace=0.3)

ax1.legend(frameon=False)
ax2.legend(frameon=False)
ax3.legend(frameon=False)
ax4.legend(frameon=False)

fig1.show()

#fig1.savefig("hist.pdf")

######################################################
################### HUBBLE DIAGRAM ###################
######################################################

fig2=plt.figure(figsize=(12,8))   
ax=fig2.add_subplot(1,1,1)
ax.set_xlabel("z")
ax.set_ylabel("m$_B$")

ax.scatter(z,mu,color="white",edgecolor="red",zorder=3)
for i, txt in enumerate(SN):
	ax.errorbar(z[i],mu[i],yerr=sig_mu[i],color="red",zorder=2)
	#ax.annotate(txt[:],(z[i]+0.001,mu[i]),size=8)

#fig2.show()

fig2.savefig("hubble1.pdf")	
