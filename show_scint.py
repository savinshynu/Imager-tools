import numpy as np
import sys, glob
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.stats import ks_2samp as ks
from mad import median_absolute_deviation as mad


filename = sys.argv[1]

hf = h5py.File(filename,'r')

im_dat =  hf.get('image')
t_dat = hf.get('time')
head_dat = hf.get('header')
el_dat = hf.get('elev')

#print t_dat[:10,:]


cyg, cas, vir, tau = [0,1,2,3]
stokes = 0
source = vir
ind  = np.where((el_dat[:,vir] > 55) & (el_dat[:,vir] < 56) & (np.int_(head_dat[:,0])==30))[0]
#ind2  = np.where((el_dat[:,cyg] > 55) & (el_dat[:,cyg] < 56) & (np.int_(head_dat[:,0])==30))[0]
#ind  = np.where( (abs(el_dat[:,cyg] - el_dat[:,cas]) < 0.5) & (np.int_(head_dat[:,0])==30))[0]
#print ind
#ind = np.arange(14300,14400)
#for x in ind:
#    print el_dat[x,cyg],el_dat[x,cas]
#ydat = np.mean(im_dat[ind,source,:,stokes],axis = 1)
#xdat = el_dat[ind,cas]

#print t_dat[ind,:]


#plt.scatter(xdat,ydat)
#plt.xlabel("Elevation in degrees")
#plt.ylabel("Power in a.u.")
#plt.show()

#plt.plot(range(im_dat.shape[2]),np.mean(im_dat[ind,source,:,stokes],axis=0))
#plt.show()
dat_source = im_dat[ind,:,:,stokes]
dev_med_cyg = np.zeros((dat_source.shape[0],dat_source.shape[2]))
dev_med_cas = np.zeros((dat_source.shape[0],dat_source.shape[2]))

med_dat_cyg = np.median(im_dat[ind,cyg,:,stokes],axis = 0)
med_dat_cas = np.median(im_dat[ind,cas,:,stokes],axis = 0)

def on_press(event):
    if event.key.isspace():
        if anim.running:
            anim.event_source.stop()
        else:
            anim.event_source.start()
        anim.running ^= True

fig,ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', on_press)
listresp = []
for i in range(dat_source.shape[0]):
    ax.set_xlabel('Channels')
    ax.set_ylabel('Response of Source')
    #ax.set_ylim(0.7,1.4)
    dev_med_cyg[i,:] = dat_source[i,cyg,:] / med_dat_cyg
    dev_med_cas[i,:] = dat_source[i,cas,:] / med_dat_cas
    std1 = np.std(dev_med_cyg[i,:])
    std2 = np.std(dev_med_cas[i,:])
    #std2 = np.std(dev_med/np.mean(dev_med))
    #stat, pval = ks(dat_source[i,:],med_dat)

    #if pval > 0.05 :#std1 < 2 :
    #if std2 <  0.02: #std_med :
    #if check(dev_med) < 10:     #== 'good':
    im1, = ax.plot(range(dat_source.shape[2]),dat_source[i,vir,:])
    #im1, = ax.plot(range(dev_med_cas.shape[1]),dev_med_cas[i,:])
    listresp.append([im1])

anim = animation.ArtistAnimation(fig, listresp, interval=500, repeat_delay=1000)
anim.running = True
anim.save('movie_vir_scint.mp4')
plt.show()




#med_dat_cyg = np.median(im_dat[ind,cyg,:,stokes],axis = 0)
#med_dat_cas = np.median(im_dat[ind,cas,:,stokes],axis = 0)
#mean_dat = np.mean(im_dat[ind,source,:,stokes],axis = 0)
plt.plot(range(med_dat_cyg.shape[0]),med_dat_cyg)
plt.plot(range(med_dat_cas.shape[0]),med_dat_cas)
plt.ylabel('Source response a.u.')
plt.xlabel('Channels')
plt.title('Median responses')
plt.show()




sig_dev_cyg = np.zeros(ind.shape)
sig_dev_cas = np.zeros(ind.shape)
#area_dat = np.zeros(ind.shape)
#error_dat = np.zeros(ind.shape)
for i,x in enumerate(ind):
    dev_cyg = im_dat[x,cyg,:,stokes] / med_dat_cyg
    dev_cas = im_dat[x,cas,:,stokes] / med_dat_cas
    #stat, pval = ks(im_dat[x,source,:,stokes],med_dat)
    #area = np.trapz(dev,range(dev.shape[0]))
    #area_dat[i] = area/198
    #cum_error =  np.sum(np.absolute(dev)) #pval, np.std(dev), area
    #error_dat[i] = cum_error
    sig_dev_cyg[i] = np.std(dev_cyg)
    sig_dev_cas[i] = np.std(dev_cas)
    #print check(dev)  #sig_dev[i]

print np.median(sig_dev_cyg),np.mean(sig_dev_cas)
std_med_cyg =  np.median(sig_dev_cyg)
std_med_cas =  np.median(sig_dev_cas)
plt.plot(range(len(sig_dev_cyg)),sig_dev_cyg,label='Cyg')
plt.plot(range(len(sig_dev_cas)),sig_dev_cas,label = 'Cas')
plt.ylabel("Standard deviation of response in each divided by median")
plt.xlabel('No. of integrations')
#plt.plot(range(len(area_dat)),area_dat)
#plt.plot(range(len(error_dat)),error_dat)
plt.legend()
plt.show()



#dat_source = im_dat[ind,:,:,stokes]
#dev_med_cyg = np.zeros((dat_source.shape[0],dat_source.shape[2]))
#dev_med_cas = np.zeros((dat_source.shape[0],dat_source.shape[2]))

fig,ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', on_press)
listresp = []
for i in range(dat_source.shape[0]):
    ax.set_xlabel('Channels')
    ax.set_ylabel('Response of Source')
    ax.set_ylim(0.7,1.4)
    dev_med_cyg[i,:] = dat_source[i,cyg,:] / med_dat_cyg
    dev_med_cas[i,:] = dat_source[i,cas,:] / med_dat_cas
    std1 = np.std(dev_med_cyg[i,:])
    std2 = np.std(dev_med_cas[i,:])
    #std2 = np.std(dev_med/np.mean(dev_med))
    #stat, pval = ks(dat_source[i,:],med_dat)

    #if pval > 0.05 :#std1 < 2 :
    if std1 <  0.02: #std_med :
    #if check(dev_med) < 10:     #== 'good':
       #im1, = ax.plot(range(dat_source.shape[1]),dat_source[i,:])
       im1, = ax.plot(range(dev_med_cyg.shape[1]),dev_med_cyg[i,:])
       listresp.append([im1])

anim = animation.ArtistAnimation(fig, listresp, interval=500, repeat_delay=1000)
anim.running = True
#ani.save('movie_cas_ut23.mp4')
plt.show()


good_cyg = np.where((sig_dev_cyg < 0.02 ))[0]
good_cas = np.where((sig_dev_cas < 0.02 ))[0]
print len(good_cyg),len(good_cas)
plt.plot(range(dat_source.shape[2]),np.mean(dat_source[good_cyg,cyg,:],axis = 0),label = 'cyg')
plt.plot(range(dat_source.shape[2]),np.mean(dat_source[good_cas,cas,:],axis = 0),label = 'cas')
#plt.ylim(0.7,1.4)
plt.ylabel('Source response a.u.')
plt.xlabel('Channels')
plt.title('Mean of responses with std < 2%')
plt.legend()
plt.show()

ratio = np.mean(dat_source[good_cas,cas,:],axis = 0)/np.mean(dat_source[good_cyg,cyg,:],axis = 0)
plt.plot(range(len(ratio)),ratio)
plt.ylim(0.7,1.4)
plt.ylabel('Ratio 0f Cas with Cyg in a.u.')
plt.xlabel('Channels')
plt.show()


par_jy = [[4.695,0.085,-0.178],[5.745,-0.770,0],[5.023,-0.856,0],[3.915,-0.299,0]]

def flux_jy_baars(v,source): # From Baars et al. paper values

    a,b,c = par_jy[source]

    log_s = a + b*np.log10(v) + c*(np.log10(v)**2)

    return 10**log_s


f1, f2,  bw = head_dat[ind[0],:3]
f2 += bw
f_chan_hz = np.arange(f1,f2,bw)

flux_cas_jy = flux_jy_baars(f_chan_hz,cas)
flux_cyg_jy = flux_jy_baars(f_chan_hz,cyg)

act_cas = ratio*flux_cyg_jy

plt.plot(f_chan_hz,act_cas)
plt.ylabel('Cas flux density in Jy')
plt.xlabel('Frequency MHz')
plt.show()
