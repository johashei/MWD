import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as npoly
import os
import seaborn as sns
import sys
import warnings

import MWD as daq

def main():
    
    try:
        trace_directory = sys.argv[1]
        energy_outfile = sys.argv[2]
    except IndexError:
        help_message = "Usage: python3 MWD_spectrum.py <directory containing traces> <file to which the energy values will be written>"
        print(help_message)
        sys.exit(0)

    # Trigger parameters
    cfd_delay = 8
    cfd_threshold = 150
    glitch_filter_threshold = 75
    # MWD parameters
    trapezoid_length = 600
    rise_time = 50
    decay_time = 40e3
    # Baseline parameters
    baseline_fit_window = 100
    baseline_extra_length = 110
    baseline_length = trapezoid_length + rise_time + baseline_extra_length
    # Parameters for fitting the decay time
    full_energy_lower_bound = 1500 # read from output histogram
    decay_length = 600
    plot = 3 # number of fits of the decay time to plot

    energies = []
    fitted_decay_time = []
    for file in os.scandir(trace_directory):
        try:
            sample, pulse = np.loadtxt(file, delimiter='\t', unpack=True)
        except ValueError:
            message = f"File {file} could not be read."
            warnings.warn(message, UserWarning)
            continue
        except Exception as ex:
            message = f"The following exception occurred when reading file {file}:"
            print(message)
            raise ex 

        cfd = daq.Trigger(pulse)
        cfd.set_parameters(cfd_delay, cfd_threshold, glitch_filter_threshold)
        cfd.find_triggers()
        trigger_sample_numbers = cfd.get_trigger_sample_numbers()
        if not trigger_sample_numbers.size:
            continue

        mwd = daq.MWD(pulse)
        mwd.set_parameters(rise_time, decay_time, trapezoid_length)
        mwd.do_mwd()
        trapezoid = mwd.get_trapezoid()

        energy_sample_numbers = trigger_sample_numbers + trapezoid_length
        baseline = daq.find_baseline(trapezoid, trigger_sample_numbers, baseline_fit_window, baseline_length)
        for sample in energy_sample_numbers[energy_sample_numbers < len(pulse)]:
            energy = trapezoid[sample] - baseline[sample]
            energies.append(energy)
            
            # Fit decay time 
            decay = pulse[sample - trapezoid_length+rise_time:sample]
            x = np.asarray(range(sample - trapezoid_length+rise_time, sample))
            fit = npoly.Polynomial.fit(x=x, y=np.log(decay), deg=1, w=np.sqrt(decay))
            fitted_decay_time.append(-1/fit.convert().coef[1])
            if plot:
                plt.plot(range(len(pulse)), pulse, label="data")
                plt.plot(x, np.exp(fit.convert().coef[0])*np.exp(fit.convert().coef[1]*x), label="fit")
                plt.legend()
                plt.show()
                plot -= 1

    np.savetxt(energy_outfile, energies)

    fig, ax = plt.subplots(1,2, figsize=(10,4))
    sns.histplot(data=energies, bins=2000, ax=ax[0])
    ax[0].set_title("spectrum")
    ax[0].set_xlabel("energy channel")
    ax[0].set_ylabel("counts", ha='right', rotation=0)
    ax[0].yaxis.set_label_coords(0,1)
    ax[0].set_xlim([200,2000])
    sns.scatterplot(x=energies, y=fitted_decay_time,  ax=ax[1])
    ax[1].set_title("decay times")
    ax[1].set_xlabel("energy channel")
    ax[1].set_ylabel("time sample", ha='right', rotation=0)
    ax[1].yaxis.set_label_coords(0,1)
    ax[1].set_xlim([-10,3500])
    ax[1].set_ylim([-0.8e5, 2.2e5])
    fig.tight_layout()
    fig.savefig("spec_and_decay.pdf", format='pdf')

    plt.show()

if __name__ == '__main__':
    main()

