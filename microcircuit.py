#!/usr/bin/env python

# Calculation of firing rates

import circuit
import cProfile
import pstats
import StringIO

circ = circuit.Circuit('microcircuit', analysis_type='stationary')
print 'firing rates', circ.th_rates

# Calculation of population rate spectra

import matplotlib.pyplot as plt
import numpy as np
import circuit

dic = {'dsd': 1.0, 'delay_dist': 'truncated_gaussian'}
circ = circuit.Circuit('microcircuit', dic)


# save_file_name = 'profile.dat'
# cProfile.run("circuit.Circuit('microcircuit', dic)", save_file_name)
# stream = StringIO.StringIO()
# p = pstats.Stats(save_file_name, stream=stream)
# 
# p.strip_dirs().sort_stats('cumtime').print_stats()
# print stream.getvalue()

freqs, power = circ.create_power_spectra()

plt.figure()
for i in range(8):
    plt.plot(freqs, np.sqrt(power[i]), label=circ.populations[i])
plt.yscale('log')
plt.legend()
plt.show()
