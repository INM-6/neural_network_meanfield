#!/usr/bin/env python

# Calculation of firing rates

import circuit

circ = circuit.Circuit('mesocircuit', analysis_type='stationary')
print 'firing rates', circ.th_rates

# Calculation of population rate spectra

import matplotlib.pyplot as plt
import numpy as np
import circuit

dic = {'dsd': 1.0, 'delay_dist': 'truncated_gaussian'}
circ = circuit.Circuit('mesocircuit', dic)
freqs, power = circ.create_power_spectra()

plt.figure(figsize=(10,10))
for i in range(8):
    plt.plot(freqs, np.sqrt(power[i]), label=circ.populations[i])
plt.yscale('log')
plt.axis('tight')
plt.legend()


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
    ax.matshow(circ.get_sensitivity_measure(freq=freq).real, cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_title(r'$f={}$ Hz'.format(freq))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
# trash unused axes
for i in range(freqs.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])

fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Im(\mathbf{Z}(f))$')
for i, freq in enumerate(freqs):
    ax = axes.flatten()[i]
    ax.matshow(circ.get_sensitivity_measure(freq=freq).imag, cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_title(r'$f={}$ Hz'.format(freq))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
# trash unused axes
for i in range(freqs.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])


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
    ax.matshow(circ.get_sensitivity_measure(freq=0, index=mode).real, cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_title(r'mode {}'.format(mode))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
# trash unused axes
for i in range(modes.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])


# compute sensitivity measura for f=0 for all modes,
# imaginary part
modes = np.arange(8)
nrows = int(np.round(np.sqrt(modes.size)))
ncols = modes.size // nrows
if modes.size % nrows > 0:
    ncols += 1
fig, axes = plt.subplots(nrows, ncols, figsize=(10,10))
fig.suptitle(r'$\Re(\mathbf{Z}(f=0))$')
for i, mode in enumerate(modes):
    ax = axes.flatten()[i]
    ax.matshow(circ.get_sensitivity_measure(freq=0, index=mode).imag, cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_title(r'mode {}'.format(mode))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
# trash unused axes
for i in range(modes.size, nrows*ncols):
    plt.delaxes(axes.flatten()[i])




freqs, eigs = circ.create_eigenvalue_spectra('MH')
# fmax = freqs[0] # look at 0-frequency
fmax = [0, 85, 330]

for f in fmax:
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
        ax.matshow(Z_amp, cmap='RdBu_r', vmin=-2, vmax=2)
        ax.set_title(r'eig_index {}'.format(i))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    # trash unused axes
    for i in range(modes.size, nrows*ncols):
        plt.delaxes(axes.flatten()[i])


for f in fmax:
    eig_index = np.arange(8)
    nrows = int(np.round(np.sqrt(eig_index.size)))
    ncols = eig_index.size // nrows
    if eig_index.size % nrows > 0:
        ncols += 1
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
        ax.matshow(Z_freq, cmap='RdBu_r', vmin=-2, vmax=2)
        ax.set_title(r'eig_index {}'.format(i))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    # trash unused axes
    for i in range(modes.size, nrows*ncols):
        plt.delaxes(axes.flatten()[i])


plt.show()
