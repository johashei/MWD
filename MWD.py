import matplotlib.pyplot as plt
from numba import jit
import numpy as np
import sys

def main():
    pulse = np.loadtxt("./test_traces/single_pulse.txt")
    
    # Set parameters
    ## CFD parameters
    cfd_delay = 8
    cfd_threshold = 150
    glitch_filter_threshold = 75
    ## MWD parameters 
    trapezoid_length = 600
    rise_time = 450
    decay_time = 5e3
    ## Baseline estimation parameters
    baseline_length = trapezoid_length + rise_time + 110
    baseline_fit_window = 100 

    cfd = Trigger(pulse, debug=True)
    cfd.set_parameters(cfd_delay, cfd_threshold, glitch_filter_threshold)
    cfd.find_triggers()
    trigger_sample_numbers = cfd.get_trigger_sample_numbers()
    energy_sample_numbers = trigger_sample_numbers + trapezoid_length

    mwd = MWD(pulse)
    mwd.set_parameters(rise_time, decay_time, trapezoid_length)
    mwd.do_mwd()
    trapezoid = mwd.get_trapezoid()

    baseline = find_baseline(trapezoid, trigger_sample_numbers, baseline_fit_window, baseline_length)
    
    energies = []
    for sample in energy_sample_numbers[energy_sample_numbers < len(pulse)]:
        energies.append(trapezoid[sample] - baseline[sample])
    
    print(f"{energies = }")
    
    x = range(len(trapezoid))
    plt.plot(x, trapezoid)
    print(trigger_sample_numbers[0,0])
    plt.axvline(trigger_sample_numbers[0,0])
    plt.axvline(trigger_sample_numbers[0,0] + rise_time)
    plt.axvline(trigger_sample_numbers[0,0] + rise_time + trapezoid_length)
    plt.axvline(energy_sample_numbers[0,0], ls=':')
    plt.show()

    # These two functions are used to check that this code gives the same result as MWD.m by James Lawson
    # See README.md for more information.
    _write_local_variables(mwd, cfd, baseline)
    _test_local_variables()

    return 0

def _write_local_variables(mwd: 'MWD', cfd: 'Trigger', baseline : np.ndarray, filename="python_local_variables.tsv"):
    test_array = np.zeros([len(mwd.pulse), 8])
    test_array[mwd.trapezoid_length:, 0] = mwd.pulse[mwd.trapezoid_length:] \
                                         - mwd.pulse[:-mwd.trapezoid_length]
    deconvoluted = mwd.mwd(mwd.pulse, mwd.decay_time, mwd.trapezoid_length)
    test_array[:, 1] = moving_window_sum(deconvoluted, window=mwd.rise_time)/mwd.decay_time
    test_array[:, 2] = deconvoluted
    test_array[:, 3] = mwd.trapezoid
    test_array[:cfd.cfd_delay, 4] = 0
    test_array[cfd.cfd_delay:, 4] = cfd.filtered_trigger_pulse
    test_array[:, 5] = baseline
    test_array[:, 6] = cfd.get_is_trigger_sample()*111
    sample_times = np.roll(cfd.get_is_trigger_sample(), mwd.trapezoid_length)
    sample_times[:mwd.trapezoid_length] = False
    test_array[:, 7] = sample_times*111
    np.savetxt(filename, np.c_[test_array], 
               fmt="%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%3d\t%3d")
    return 

def _test_local_variables(reference_file="matlab_local_variables.tsv", file_to_test="python_local_variables.tsv"):
    reference = np.loadtxt(reference_file, delimiter='\t')
    rewritten = np.loadtxt(file_to_test, delimiter='\t')
    errors = rewritten-reference
    print(f"the largest errors in each variable are:")
    sys.stdout.write("  ")
    for error in np.max(errors,0):
        sys.stdout.write(f" {error:8.2e} ") 
    sys.stdout.write("\nin")
    for i,j in enumerate(np.argmax(errors,0)):
        sys.stdout.write(f" {rewritten[j, i]:8.3f} ")
    sys.stdout.write("\nat")
    for idx in np.argmax(errors,0):
        sys.stdout.write(f" {idx:8d} ")
    return 


def find_baseline(trapezoid, trigger_sample_numbers, baseline_fit_window, baseline_length):
    baseline = np.copy(trapezoid)
    for trigger_sample_number in trigger_sample_numbers[trigger_sample_numbers >= baseline_fit_window]:
        baseline[trigger_sample_number:trigger_sample_number + baseline_length] = trapezoid[trigger_sample_number]
    return baseline


class Trigger(object):
    def __init__(self, pulse, debug=False):
        self.pulse = pulse
        self.debug = debug
        self.is_trigger_sample = np.zeros_like(self.pulse, dtype=bool)
        self.filtered_trigger_pulse = None # this is only filled if debug is True

    def set_parameters(self, cfd_delay, cfd_threshold, glitch_filter_threshold): 
        self.trigger_pulse = abs(self.pulse[cfd_delay:] - self.pulse[:-cfd_delay])
        self.cfd_delay = cfd_delay
        self.cfd_threshold = cfd_threshold
        self.glitch_filter_threshold = glitch_filter_threshold
      
    def get_trigger_sample_numbers(self):
        return np.array(self.is_trigger_sample.nonzero())
    
    def get_is_trigger_sample(self):
        return self.is_trigger_sample
    
    def find_triggers(self):
        filtered_trigger_pulse = Trigger.glitch_filter(self.trigger_pulse, self.glitch_filter_threshold)
        is_trigger_event = Trigger.trigger(filtered_trigger_pulse, self.cfd_threshold)
        self.is_trigger_sample[self.cfd_delay:] = is_trigger_event
        if self.debug:
            self.filtered_trigger_pulse = filtered_trigger_pulse

    @staticmethod
    @jit(nopython=True)
    def glitch_filter(unfiltered_pulse, filter_threshold):
        filtered_pulse = np.empty_like(unfiltered_pulse)
        filtered_pulse[0] = unfiltered_pulse[0]
        for i, amplitude in enumerate(unfiltered_pulse[1:]):
            if abs(amplitude - filtered_pulse[i]) > filter_threshold:
                filtered_pulse[i+1] = amplitude
            else:
                filtered_pulse[i+1] = filtered_pulse[i]
        return(filtered_pulse)

    @staticmethod
    @jit(nopython=True)
    def trigger(trigger_pulse, trigger_threshold):
        is_trigger_event = np.empty_like(trigger_pulse, dtype=np.bool_) # np type for jit to work 
        is_trigger_event[0] = False
        for i, (previous_element, current_element) in enumerate(zip(trigger_pulse[:-1], \
                                                                  trigger_pulse[1:])):
            is_trigger_event[i+1] = current_element >= trigger_threshold and \
                                  previous_element < trigger_threshold
        return(is_trigger_event)


class MWD(object):
    def __init__(self, pulse):
        self.pulse = pulse

    def read_input_parameters(self, input_file_name):
        message = "Not implemented yet. Please call the 'set_parameters' method."
        raise NotImplementedError(message)

    def set_parameters(self, rise_time, decay_time, trapezoid_length):
        self.rise_time = rise_time
        self.decay_time = decay_time
        self.trapezoid_length = trapezoid_length

    def get_trapezoid(self):
        return self.trapezoid

    def do_mwd(self):
        deconvoluted = self.mwd(self.pulse, self.decay_time, self.trapezoid_length)
        self.trapezoid = moving_window_sum(deconvoluted, window=self.rise_time)/self.rise_time 

    @staticmethod
    @jit(nopython=True)
    def mwd(pulse, decay_time, trapezoid_length):
        deconvoluted = np.zeros_like(pulse)
        deconvoluted[trapezoid_length:] = pulse[trapezoid_length:] \
                                        + moving_window_sum(pulse, trapezoid_length)[trapezoid_length:]/decay_time \
                                        - pulse[:-trapezoid_length]
        return (deconvoluted)


# must be outside class to be called from another jitted function
@jit(nopython=True)
def moving_window_sum(array, window):
    output = np.zeros_like(array)
    output[window] = np.sum(array[:window])
    for i in range(window,len(array)-1):
        # Instead of recomputing the sum in each iteration, 
        # subtract the first element and add the new one.
        output[i+1] = output[i] - array[i-window] + array[i]
    return output



if __name__ == '__main__':
    main()








