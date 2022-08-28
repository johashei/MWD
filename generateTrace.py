import numpy as np

class Trace(object):
    def __init__(self, riseTime, decayTime):
        self.riseTime = riseTime
        self.decayTime = decayTime
    
    def Generate(pulseStartTimes, traceLength):
        self.signal = np.zeros(traceLength)
        for time in pulseStartTimes:
            self.AddPulse(time)

    def AddPulse(self, start):
       pass 
