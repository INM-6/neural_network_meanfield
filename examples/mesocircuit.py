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

# parameters
dic = {'dsd': 1.0,          # not sure what this does, not used anywhere that I can find
       'delay_dist': 'truncated_gaussian',
       # 'Next' : np.array([1796, 1686, 2000, 2064, 1938, 1941, 3296, 2382]), # noise indegrees
       'Next' : np.array([1796, 1686, 2200, 2064, 1938, 1941, 3296, 2382]), # noise indegrees
       # 'Next' : np.array([1796, 1686, 2374, 2064, 1938, 1941, 3296, 2382]),  # noise indegrees
       'w' : 87.8*0.5,      # PSC amplitude in pA, mul by time constant of Potjans model
       'tauf' : 0.5         # synapse time constant in ms
       }
# get rate predictions
circ = Circuit('mesocircuit', dic, analysis_type='stationary')
print 'firing rates', circ.th_rates

# get power spectra and sensitivity measures
circ = Circuit('mesocircuit', dic, analysis_type='dynamical')

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
