import matplotlib.pyplot as plt
import numpy as np

def main():
    names = ["single_pulse", "multiple_pulses"]
    for name in names:
        fig, ax = plt.subplots(1,2, sharey=True)
        signal = np.loadtxt(f"./test_traces/{name}.dat")
        [deconvoluted, trapezoid, filtered_trigger_pulse, baseline, is_trigger_sample, sample_times] = np.loadtxt(f"./test_traces/{name}_reference.tsv", delimiter='\t', usecols=[2,3,4,5,6,7], unpack=True)
        ch = np.asarray(range(len(signal)))
        
        ax[0].set_title("Moving Window Deconvolution")
        ax[0].set_xlabel("time channel")
        ax[0].plot(ch, signal, label="signal")
        ax[0].plot(ch, deconvoluted, label="deconvoluted")
        ax[0].plot(ch, trapezoid, label="trapezoid")
        ax[0].plot(ch, baseline, label="baseline")

        ax[1].set_title("Constant Fraction Discriminator")
        ax[1].set_xlabel("time channel")
        ax[1].plot(ch, signal, label="signal")
        ax[1].plot(ch[:-8], signal[8:], label="delayed signal")
        ax[1].plot(ch[:-8], signal[8:] - signal[:-8], label="trigger pulse") 
        ax[1].plot(ch, filtered_trigger_pulse, label="filtered trigger pulse")
        
        ax[1].set_ylim(ax[1].get_ylim())
        ax[1].vlines(np.nonzero(is_trigger_sample), *ax[0].get_ylim(), ls=':', label="trigger")
        ax[0].vlines(np.nonzero(sample_times), *ax[0].get_ylim(), ls=':', label="readout")
        ax[0].legend()
        ax[1].legend()

        fig.tight_layout()

    plt.show()

if __name__ == '__main__':
    main()
