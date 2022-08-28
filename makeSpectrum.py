import MWD as daq
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as npoly
import os
import seaborn as sns
import sys
import warnings


def main():
    
    try:
        traceDirectory = sys.argv[1]
        energyOutFile = sys.argv[2]
    except IndexError:
        helpMessage = "Usage: python3 MWD_spectrum.py <directory containing traces> <file to which the energy values will be written>"
        print(helpMessage)
        sys.exit(0)

    # Trigger parameters
    cfdDelay = 8
    cfdThreshold = 150
    glitchFilterThreshold = 75
    # MWD parameters
    trapezoidLength = 600
    riseTime = 50
    decayTime = 40e3
    # Baseline parameters
    BLFL = 100
    baselineExtraLength = 110
    baselineLength = trapezoidLength + riseTime + baselineExtraLength
    # Parameters for fitting the decay time
    fullEnergyLowerBound = 1500 # read from output histogram
    decayLength = 600
    plot = 3

    energies = []
    fittedDecayTime = []
    for file in os.scandir(traceDirectory):
        # TODO: change so it can read from binary file format
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
        cfd.Set_parameters(cfdDelay, cfdThreshold, glitchFilterThreshold)
        cfd.findTriggers()
        triggerSampleNumbers = cfd.Get_triggerSampleNumbers()
        if not triggerSampleNumbers.size:
            continue

        mwd = daq.MWD(pulse)
        mwd.Set_parameters(riseTime, decayTime, trapezoidLength)
        mwd.doMWD()
        trapezoid = mwd.Get_trapezoid()

        energySampleNumbers = triggerSampleNumbers + trapezoidLength
        baseline = daq.findBaseline(trapezoid, triggerSampleNumbers, BLFL, baselineLength)
        for sample in energySampleNumbers[energySampleNumbers < len(pulse)]:
            energy = trapezoid[sample] - baseline[sample]
            energies.append(energy)
            
            # Fit decay time 
            decay = pulse[sample - trapezoidLength+riseTime:sample]
            x = np.asarray(range(sample - trapezoidLength+riseTime, sample))
            fit = npoly.Polynomial.fit(x=x, y=np.log(decay), deg=1, w=np.sqrt(decay))
            fittedDecayTime.append(-1/fit.convert().coef[1])
            if plot:
                plt.plot(range(len(pulse)), pulse, label="data")
                plt.plot(x, np.exp(fit.convert().coef[0])*np.exp(fit.convert().coef[1]*x), label="fit")
                plt.legend()
                plt.show()
                plot -= 1

    np.savetxt(energyOutFile, energies)

    fig, ax = plt.subplots(1,2)
    sns.histplot(data=energies, bins=2000, ax=ax[0])
    ax[0].set_title("spectrum")
    sns.scatterplot(x=energies, y=fittedDecayTime,  ax=ax[1])
    ax[1].set_title("decay times")
    plt.show()

if __name__ == '__main__':
    main()

