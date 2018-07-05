import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import os

plt.rcdefaults()
plt.rc('lines', linewidth=1.5)
plt.rc('axes', labelsize=26)
plt.rc('xtick', labelsize=26)
plt.rc('ytick', labelsize=26)
plt.rc("legend", fontsize=15)  

##### Cut parameters

x1_min = -5
x1_max = 5

epming = 0
epminr = 0

chisqred_max = 5

ngmin = 1
nrmin = 1

### Histograms

dz = 0.01
dep = 1
dx1 = 0.5
dc = 0.05

# Upload parameters from the most recent fit 
files = []
path = "runs/"
for which in os.listdir(path):
	if os.path.isfile(os.path.join(path,which)) and "fitall" in which:
		files.append(str(which))
filename = files[-1]


fig=plt.figure(figsize=(18,10))   
ax1=fig.add_subplot(2,2,1)
ax1.set_xlabel("Redshift")
ax2=fig.add_subplot(2,2,2)
ax2.set_xlabel("Phase")
ax3=fig.add_subplot(2,2,3)
ax3.set_xlabel("x1")
ax4=fig.add_subplot(2,2,4)
ax4.set_xlabel("c")

SN = []
z = []
t0 = []
x0 = []
x1 = []
c = []
ebv_mw = []
rv_mw = []
mu = []
sig_mu = []
epg = []
epr = []
ng = []
nr = []
chisq = []
ndof = []

a = open(path+filename)
for i, line in enumerate(a):

	if i==0:
		Nmarsh = float(line.split(' ')[1])
	if not line.startswith("#"):
		SN.append(str(line.split(' ')[0]))
		z.append(float(line.split(' ')[1]))
		t0.append(float(line.split(' ')[2]))
		x0.append(float(line.split(' ')[3]))
		x1.append(float(line.split(' ')[4]))
		c.append(float(line.split(' ')[5]))
		ebv_mw.append(float(line.split(' ')[6]))
		rv_mw.append(float(line.split(' ')[7]))
		mu.append(float(line.split(' ')[8]))
		sig_mu.append(float(line.split(' ')[9]))
		epg.append(float(line.split(' ')[10]))
		epr.append(float(line.split(' ')[11]))
		ng.append(float(line.split(' ')[12]))
		nr.append(float(line.split(' ')[13]))
		chisq.append(float(line.split(' ')[14]))
		ndof.append(float(line.split(' ')[15]))

chisq_red = np.array(chisq) / ndof

print "%s are missing g-band photometry" %len(np.array(epg)[np.where(np.array(epg)==1e6)])
print "%s are missing r-band photometry" %len(np.array(epr)[np.where(np.array(epr)==1e6)])

### Plot histograms

minz = 0
maxz = max(z)
nz = int(maxz/dz)+1
ax1.hist(z,bins=arange(-dz/2.,dz*(nz+1/2.),dz),label="%d of %d SNe"%(len(z),Nmarsh))

### Select good data
w=list(zip(SN,x1,c,epg,epr,ng,nr,x0,z,mu,sig_mu,chisq_red))
ww=list(filter(lambda r:r[1]>=x1_min and r[1]<=x1_max 
					    and r[3]<=epming and r[4]<=epminr 
					    and r[5]>=ngmin and r[6]>=nrmin
					    and r[11]<=chisqred_max,w))
SN,x1,c,epg,epr,ng,nr,x0,z,mu,sig_mu,chisqnew = list(zip(*ww))

firstep = min(epg,epr)

#print np.array(SN)[np.where(np.array(firstep)<-14)]

##########

absmax = max(-min(firstep),max(firstep))
nep = int(absmax/dep)+1
ax2.hist(firstep,bins=arange(-dep*(nep+1/2.),dep*(nep+1/2.),dep),label="%d of %d SNe"%(len(firstep),Nmarsh))

absmax = max(-min(x1),max(x1))
nx1 = int(absmax/dx1)+1
ax3.hist(x1,bins=arange(-dx1*(nx1+1/2.),dx1*(nx1+1/2.),dx1),label="%d of %d SNe"%(len(x1),Nmarsh))

absmax = max(-min(c),max(c))
nc = int(absmax/dc)+1
ax4.hist(c,bins=arange(-dc*(nc+1/2.),dc*(nc+1/2.),dc),label="%d of %d SNe"%(len(c),Nmarsh))


ax1.set_xlim(0,max(z)+dz/2.)
ax2.set_xlim(-22,max(firstep))
ax3.set_xlim(x1_min,x1_max)
ax4.set_xlim(-max(-min(c),max(c)),max(-min(c),max(c)))

subplots_adjust(hspace=0.3)

ax1.legend(frameon=False)
ax2.legend(frameon=False)
ax3.legend(frameon=False)
ax4.legend(frameon=False)

#fig.show()

#fig.savefig("hist.pdf")

######### HUBBLE DIAGRAM

fig2=plt.figure(figsize=(12,8))   
ax=fig2.add_subplot(1,1,1)
ax.set_xlabel("z")
ax.set_ylabel("m$_B$")
ax.set_xlim(0.01,0.17)
ax.set_ylim(15.1,20.9)

ax.scatter(z,mu,color="white",edgecolor="red",zorder=3)
for i, txt in enumerate(SN):
	ax.errorbar(z[i],mu[i],yerr=sig_mu[i],color="red",zorder=2)

	#ax.annotate(txt[:],(z[i]+0.001,mu[i]),size=8)

fig2.show()

#fig2.savefig("hubble1.pdf")

