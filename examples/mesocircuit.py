#!/usr/bin/env python
'''Calculation of population rates and rate spectra'''

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import os
from pymeanfield import Circuit

# plot tuning
plt.rcParams.update({
    'axes.titlesize': u'small'
})

def colorbar(fig, ax, im):
    rect = np.array(ax.get_position().bounds)
    rect[0] += rect[2] + 0.01
    rect[2] = 0.01
    cax = fig.add_axes(rect)
    cbar = plt.colorbar(im, cax=cax)
    cbar.locator = MaxNLocator(nbins=4)
    cbar.update_ticks()


def get_mesocircuit_params():
    '''returns parameters adapted from the microcircuit parameterset'''
    params = {}
    
    params['populations'] = ['23E', '23I', '4E', '4I', 
                             '5E', '5I', '6E', '6I']
    # number of neurons in populations
    params['N'] = np.array([330928, 93344, 350640, 87664, 77600, 17040, 230320, 47168])
        
    ### Neurons
    params['C'] = 250.0    # membrane capacitance in pF
    params['taum'] = 10.0  # membrane time constant in ms
    params['taur'] = 2.0   # refractory time in ms
    params['V0'] = -65.0   # reset potential in mV
    params['Vth'] = -50.0  # threshold of membrane potential in mV
    
    ### Synapses
    params['tauf'] = 0.5  # synaptic time constant in ms
    params['de'] = 1.5    # delay of excitatory connections in ms
    params['di'] = 0.75   # delay of inhibitory connections in ms
    # standard deviation of delay of excitatory connections in ms
    params['de_sd'] = params['de']*0.5 
    # standard deviation of delay of inhibitory connections in ms
    params['di_sd'] = params['di']*0.5
    # delay distribution, options: 'none', 'gaussian' (standard deviation 
    # is defined above), 'truncated gaussian' (standard deviation is 
    # defined above, truncation at zero)
    params['delay_dist'] = 'truncated_gaussian' #'none'
    # PSC amplitude in pA
    params['w'] = 87.8*0.5 # 0.5 being the default time constant of the microcircuit
    
    ### Connectivity 
    # indegrees
    params['I'] = (np.array([
        [908052706, 429956872, 415727848, 194859668, 67971857, 0, 47420816, 0],
        [342156175, 98313800, 84752610, 34561135, 44894844, 0, 7390917, 0],
        [73143801, 15807409, 501095081, 341493367, 14923606, 146694, 299954115, 0],
        [164520432, 1942287, 200388899, 101101409, 1837434, 0, 175415933, 0],
        [211873205, 36965998, 112686078, 3174062, 41092415, 40948453, 29863645, 0],
        [25346349, 3505443, 12579906, 268973, 6508893, 7602929, 2763177, 0],
        [97369358, 11615440, 139575131, 27448088, 83861841, 6332481, 172194137, 202086287],
        [46567837, 360336, 4602613, 169201, 8306879, 526387, 58658806, 26429834]]).T / params['N'].astype(float)).T
    # ratio of inhibitory to excitatory weights
    params['g']=4.0
    
    ### External input
    params['v_ext'] = 8.0 # in Hz
    # number of external neurons
    params['Next'] = np.array([1796, 1686, 2200, 2064, 1938, 1941, 3296, 2382])
    
    ### Neural response
    # Transfer function is either calculated analytically ('analytical')
    # or approximated by an exponential ('empirical'). In the latter case
    # the time constants in response to an incoming impulse ('tau_impulse'), 
    # as well as the instantaneous rate jumps ('delta_f') have to be
    # specified.
    params['tf_mode'] = 'analytical'   
    # number of modes used when fast response time constants are calculated
    params['num_modes'] = 1
           
    # file storing results
    params['datafile'] = 'results_mesocircuit.h5'

    return params

params = get_mesocircuit_params()

# get rate predictions
circ = Circuit('microcircuit', params, analysis_type='stationary')
print 'firing rates', circ.th_rates

# get power spectra and sensitivity measures
circ = Circuit('microcircuit', params, analysis_type='dynamical')

# prep folder for figures output
figdir = os.path.join('mesocircuit_figures', circ.param_hash)
if not os.path.isdir(os.path.split(figdir)[0]):
    os.mkdir(os.path.split(figdir)[0])
if not os.path.isdir(figdir):
    os.mkdir(figdir)
    

# plot power spectra
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
freqs, power = circ.create_power_spectra()
for i in range(8):
    ax.semilogy(freqs, np.sqrt(power[i]), label=circ.populations[i])
ax.axis(ax.axis('tight'))
ax.set_xlabel(r'$f$ (Hz)')
ax.set_ylabel(r'rate PSD (Hz$^2$/s)')
ax.legend(loc='best', fontsize='small')
ax.grid('on')
fig.savefig(os.path.join(figdir, 'mesocircuit_PSD.pdf'), bbox_inches='tight')
# plt.close(fig)


# compute and plot sensitivity measure as in Bos et al for a range of frequencies
# using the dominant mode
freqs = np.arange(51)*10
nrows = int(np.round(np.sqrt(freqs.size)))
ncols = freqs.size // nrows
if freqs.size % ncols > 0:
    nrows += 1
fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Re(\mathbf{Z}(f))$')
for i, freq in enumerate(freqs):
    ax = axes.flatten()[i]
    im = ax.matshow(circ.get_sensitivity_measure(freq=freq).real, cmap='RdBu_r', vmin=-2, vmax=2)
    # ax.set_title(r'$f={}$ Hz'.format(freq))
    ax.set_title('{} Hz'.format(freq), va='top')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
colorbar(fig, ax, im)
# trash unused axes
for i in range(freqs.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])
fig.savefig(os.path.join(figdir, 'mesocircuit_Re(Z(f)).pdf'), bbox_inches='tight')    
# plt.close(fig)


fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Im(\mathbf{Z}(f))$')
for i, freq in enumerate(freqs):
    ax = axes.flatten()[i]
    im = ax.matshow(circ.get_sensitivity_measure(freq=freq).imag, cmap='RdBu_r', vmin=-2, vmax=2)
    # ax.set_title(r'$f={}$ Hz'.format(freq))
    ax.set_title('{} Hz'.format(freq), va='top')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
colorbar(fig, ax, im)
# trash unused axes
for i in range(freqs.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])
fig.savefig(os.path.join(figdir, 'mesocircuit_Im(Z(f)).pdf'), bbox_inches='tight')    
# plt.close(fig)


fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Im(\mathbf{Z}(f))$')
fig.suptitle(r'$\mathbf{Z}^\mathrm{amp}(f)$')
_, eigs = circ.create_eigenvalue_spectra('MH') # first returned arg is freqs
for i, freq in enumerate(freqs):
    index = 0
    # get eigenvalue closest to instablilty
    eigc = eigs[index][np.argmin(abs(eigs[index]-1))]
    
    Z = circ.get_sensitivity_measure(freq)
    k = np.asarray([1, 0])-np.asarray([eigc.real, eigc.imag])
    k /= np.sqrt(np.dot(k, k))
    k_per = np.asarray([-k[1], k[0]])
    k_per /= np.sqrt(np.dot(k_per, k_per))
    Z_amp = Z.real*k[0]+Z.imag*k[1]
    Z_freq = Z.real*k_per[0]+Z.imag*k_per[1]

    ax = axes.flatten()[i]
    im = ax.matshow(Z_amp, cmap='RdBu_r', vmin=-2, vmax=2)
    # ax.set_title(r'$f={}$ Hz'.format(freq))
    ax.set_title('{} Hz'.format(freq), va='top')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
colorbar(fig, ax, im)
# trash unused axes
for i in range(freqs.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])
fig.savefig(os.path.join(figdir, 'mesocircuit_Z_amp(f)).pdf'), bbox_inches='tight')    
# plt.close(fig)


fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Im(\mathbf{Z}(f))$')
fig.suptitle(r'$\mathbf{Z}^\mathrm{freq}(f)$')
for i, freq in enumerate(freqs):
    index = 0
    # get eigenvalue closest to instablilty
    eigc = eigs[index][np.argmin(abs(eigs[index]-1))]
        
    Z = circ.get_sensitivity_measure(freq)
    k = np.asarray([1, 0])-np.asarray([eigc.real, eigc.imag])
    k /= np.sqrt(np.dot(k, k))
    k_per = np.asarray([-k[1], k[0]])
    k_per /= np.sqrt(np.dot(k_per, k_per))
    Z_amp = Z.real*k[0]+Z.imag*k[1]
    Z_freq = Z.real*k_per[0]+Z.imag*k_per[1]

    ax = axes.flatten()[i]
    im = ax.matshow(Z_freq, cmap='RdBu_r', vmin=-2, vmax=2)
    # ax.set_title(r'$f={}$ Hz'.format(freq))
    ax.set_title('{} Hz'.format(freq), va='top')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
colorbar(fig, ax, im)
# trash unused axes
for i in range(freqs.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])
fig.savefig(os.path.join(figdir, 'mesocircuit_Z_freq(f)).pdf'), bbox_inches='tight')    
# plt.close(fig)


# compute sensitivity measura for f=0 for all modes,
# real part. 
modes = np.arange(8)
nrows = int(np.round(np.sqrt(modes.size)))
ncols = modes.size // nrows
if modes.size % nrows > 0:
    ncols += 1
fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Re(\mathbf{Z}(f=0))$')
for i, mode in enumerate(modes):
    ax = axes.flatten()[i]
    im = ax.matshow(circ.get_sensitivity_measure(freq=0, index=mode).real, cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_title(r'mode {}'.format(mode), va='top')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
colorbar(fig, ax, im)
# trash unused axes
for i in range(modes.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])
fig.savefig(os.path.join(figdir, 'mesocircuit_Re(Z(f=0))_all_modes.pdf'), bbox_inches='tight')    
# plt.close(fig)


# compute sensitivity measura for f=0 for all modes,
# imaginary part
modes = np.arange(8)
nrows = int(np.round(np.sqrt(modes.size)))
ncols = modes.size // nrows
if modes.size % nrows > 0:
    ncols += 1
fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Im(\mathbf{Z}(f=0))$')
for i, mode in enumerate(modes):
    ax = axes.flatten()[i]
    im = ax.matshow(circ.get_sensitivity_measure(freq=0, index=mode).imag, cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_title(r'mode {}'.format(mode), va='top')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
colorbar(fig, ax, im)
# trash unused axes
for i in range(modes.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])
fig.savefig(os.path.join(figdir, 'mesocircuit_Im(Z(f=0))_all_modes.pdf'), bbox_inches='tight')    
# plt.close(fig)




_, eigs = circ.create_eigenvalue_spectra('MH') # first returned arg is freqs
freqs = [0, 85, 280, 330]

for f in freqs:
    eig_index = np.arange(8)
    nrows = int(np.round(np.sqrt(eig_index.size)))
    ncols = eig_index.size // nrows
    if eig_index.size % nrows > 0:
        ncols += 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
    fig.suptitle(r'$\mathbf{{Z}}^\mathrm{{amp}}(f={})$'.format(f))
    
    for i in eig_index:
        # get eigenvalue closest to instablilty
        eigc = eigs[i][np.argmin(abs(eigs[i]-1))]
        
        
        Z = circ.get_sensitivity_measure(f)
        k = np.asarray([1, 0])-np.asarray([eigc.real, eigc.imag])
        k /= np.sqrt(np.dot(k, k))
        k_per = np.asarray([-k[1], k[0]])
        k_per /= np.sqrt(np.dot(k_per, k_per))
        Z_amp = Z.real*k[0]+Z.imag*k[1]
        Z_freq = Z.real*k_per[0]+Z.imag*k_per[1]
    
        ax = axes.flatten()[i]
        im = ax.matshow(Z_amp, cmap='RdBu_r', vmin=-2, vmax=2)
        ax.set_title(r'eig_index {}'.format(i))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    colorbar(fig, ax, im)
    # trash unused axes
    for i in range(modes.size, nrows*ncols):
        plt.delaxes(axes.flatten()[i])
    fig.savefig(os.path.join(figdir, 'mesocircuit_Z_amp(f={})_all_modes).pdf').format(f), bbox_inches='tight')    
    # plt.close(fig)
    

    # eig_index = np.arange(8)
    # nrows = int(np.round(np.sqrt(eig_index.size)))
    # ncols = eig_index.size // nrows
    # if eig_index.size % nrows > 0:
    #     ncols += 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
    fig.suptitle(r'$\mathbf{{Z}}^\mathrm{{freq}}(f={})$'.format(f))
    
    for i in eig_index:
        # get eigenvalue closest to instablilty
        eigc = eigs[i][np.argmin(abs(eigs[i]-1))]
        
        
        Z = circ.get_sensitivity_measure(f)
        k = np.asarray([1, 0])-np.asarray([eigc.real, eigc.imag])
        k /= np.sqrt(np.dot(k, k))
        k_per = np.asarray([-k[1], k[0]])
        k_per /= np.sqrt(np.dot(k_per, k_per))
        Z_amp = Z.real*k[0]+Z.imag*k[1]
        Z_freq = Z.real*k_per[0]+Z.imag*k_per[1]
    
        ax = axes.flatten()[i]
        im = ax.matshow(Z_freq, cmap='RdBu_r', vmin=-2, vmax=2)
        ax.set_title(r'eig_index {}'.format(i))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    colorbar(fig, ax, im)
    # trash unused axes
    for i in range(modes.size, nrows*ncols):
        plt.delaxes(axes.flatten()[i])
    fig.savefig(os.path.join(figdir, 'mesocircuit_Z_freq(f={})_all_modes).pdf').format(f), bbox_inches='tight')    
    # plt.close(fig)
        


plt.show()
