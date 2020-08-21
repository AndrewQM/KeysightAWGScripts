import matplotlib.pyplot as plt
import numpy as np
from collections import deque

def gaussian(sigma, center, time):
    # sigma = pulse_width / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-(time - center) ** 2 / (2 * sigma ** 2))




#each pulse is it's own object with correponding parameters like hightime, sigma, and amplitude
class pulse():
    def __init__(self, hightime, risetime, center, sampleTime, amplitude, Data1_, Data2_):
        self.hightime = hightime
        self.risetime = risetime
        self.center = center
        self.sampleTime = sampleTime
        self.amplitude = amplitude
        self.left_center = 0
        self.right_center = 0
        self.nearCenterSampleIndex = 0
        self.Data1 = Data1_
        self.Data2 = Data2_
                                                    ##########global?
        self.nearCenterSampleIndex = self.center // sampleTime
        self.index_left = int(self.nearCenterSampleIndex)
        self.index_right = int(self.nearCenterSampleIndex + 1)
        #print(self.nearCenterSampleIndex)

    def writePulseGaussian(self, time):
        # defines the theoretical (smooth, unsampled) pulse with risetime, hightime, and center.
        # writePulseFragment() calls this
        if (self.center - 0.5 * self.hightime) < time < (self.center + 0.5 * self.hightime):
            # print("one")
            amp = self.amplitude
        elif time <= (self.center - 0.5 * self.hightime):
            # print("two")
            self.left_center = round((self.center - 0.5 * self.hightime) / self.sampleTime) * self.sampleTime

            #print("left_center divided by sample: ", self.left_center / self.sampleTime)
            #print("time divided by sampleTIme: ", time / self.sampleTime)
            #print("sampled left center: ", self.left_center)
            amp = self.amplitude * gaussian(self.risetime, self.left_center, time)
        elif time >= (self.center + 0.5 * self.hightime):
            # print("three")
            self.right_center = round((self.center + 0.5 * self.hightime) / self.sampleTime) * self.sampleTime

            amp = self.amplitude * gaussian(self.risetime, self.right_center, time)
        else:
            # print("error")
            exit(1)
        return amp

    def writeSampledPulse(self, sequence_):
        #Writes several samples into the main sequence that form a pulse

        minLeftAmplitude = 1
        minRightAmplitude = 1
        req = 0
        default_index_left = self.index_left
        default_index_right= self.index_right
        # this should go inside the pulse class
        # write the pulse to the list iterating out from near the center
        while minLeftAmplitude > .002 or minRightAmplitude > .002:
            # print("location of pulse between indexes:", pulse_time/sampleTime)
            sequence_[self.index_left] = self.writePulseGaussian(self.index_left * self.sampleTime)
            if sequence_[self.index_left] <= minLeftAmplitude:
                minLeftAmplitude = sequence_[self.index_left]
            self.index_left = self.index_left - 1
            # print("self.index_left: ", self.index_left)

            sequence_[self.index_right] = self.writePulseGaussian(self.index_right * self.sampleTime)
            if sequence_[self.index_right] <= minRightAmplitude:
                minRightAmplitude = sequence_[self.index_right]
            self.index_right = self.index_right + 1
            # print("self.index_right: ", self.index_right)

            req = req + 1
            # print("cycle")
            if req > 1000:
                print("cycle error")
                exit(1)

        self.index_left = default_index_left
        self.index_right = default_index_right

    def plotSampledPulse(self, sequence_, count_, laserTime_, sampleTime_, figcount_):


        #be aware laserTime and sampleTime are still global
        x = []
        y = []
        for i in range(self.index_left - int(2 * laserTime_ // sampleTime_),
                       self.index_right + int(2 * laserTime_ // sampleTime_)): x.append(i * sampleTime_ * 1e9)

        for i in range(self.index_left - int(2 * laserTime_ // sampleTime_),
                       self.index_right + int(2 * laserTime_ // sampleTime_)): y.append(sequence_[i])


        plt.figure(figcount_)
        figcount_ = figcount_ + 1
        plt.axvline(self.left_center * 1e9, color="#d1d1d1")
        plt.axvline(self.right_center * 1e9, color="#d1d1d1")

        plt.plot(x, y, color="#4189c4")

        plt.axvline(self.center * 1e9, color="#ffc4b5")
        plt.axvline((self.center + laserTime_) * 1e9, color="#ffc4b5")
        plt.axvline((self.center - laserTime_) * 1e9, color="#ffc4b5")
        plt.axvline((self.center + 2 * laserTime_) * 1e9, color="#ffc4b5")
        plt.axvline((self.center - 2 * laserTime_) * 1e9, color="#ffc4b5")
        # print("pulse laser time: ", pulse_time*1e9)
        # print("pulse laser time: ", pulse_time*1e9)
        # print("pulse laser time plus 1: ", pulse_time*1e9)

        plt.xlabel("time (ns)")
        plt.ylabel("Amplitude")
        plt.title("pulse in cycle number " + str(count_ + 1))
        plt.show()




class PulseSequence():
    def __init__(self, Data1_, Data2_, sampleTime_, laserTime_, deadPulses_, sequenceLength_, usableSymbols_, padding_ = 0):
        assert len(Data1_) <= usableSymbols_
        assert len(Data2_) <= usableSymbols_
        self.pulseList = [] #this will be a list of pulse objects, each with their own center, hightime, risetime, etc. specified.
        self.times = [0] * len(Data1_)
        self.Data1 = Data1_
        self.Data2 = Data2_
        self.padding = padding_
        self.sampleTime = sampleTime_
        self.laserTime = laserTime_
        self.deadPulses = deadPulses_
        self.figcount = 1
        self.sequence1 = [0] * sequenceLength_
        self.sequence2 = [0] * sequenceLength_
        self.zerosequence = [0]*sequenceLength_
        self.clock_sequence = [0] * sequenceLength_

        for i in range(sequenceLength_ // 10):  # just sets the beginning of the sequence to 1, rest is zero.
            self.clock_sequence[i] = 1


    def generateTimesList(self):
        # generate at list of universal AWG pulse times from the start time of each cycle
        # these AWG pulse times align with pulses of the laser
        self.cycleStart = self.padding
        for i in range(len(self.times)):
            self.times[i] = self.cycleStart
            self.cycleStart = self.cycleStart + (1 + self.deadPulses)*self.laserTime
        #print(len(self.times))

    def genDefaultPulseObjects(self, hightime_, sigma_, amplitude_, Data1_, Data2_):
        #this writes a list of pulse objects
        for i in range(len(self.Data1)):
            self.pulseList.append(pulse(hightime_, sigma_, self.times[i], self.sampleTime, amplitude_, Data1_[i], Data2_[i]))
        #print(self.pulseList)

    # Where symbols get transformed into pulses in different channels
    def writePulses(self):
        for pulseObject in self.pulseList:

            if pulseObject.Data1 == 1 and pulseObject.Data2 == 0: # "early"
                pulseObject.writeSampledPulse(self.sequence1)
            if pulseObject.Data1 == 0 and pulseObject.Data2 == 1: # "late"
                pulseObject.writeSampledPulse(self.sequence2)
            if pulseObject.Data1 == 1 and pulseObject.Data2 == 1: # "phase basis"
                pulseObject.amplitude = 0.5
                pulseObject.writeSampledPulse(self.sequence1)
                pulseObject.writeSampledPulse(self.sequence2)

    def plotSomePulses(self, start, end):
        for count, pulseObject in enumerate(self.pulseList[start:end]):
            pulseObject.plotSampledPulse(self.sequence1, count, self.laserTime, self.sampleTime, self.figcount)
            self.figcount = self.figcount + 1
            pulseObject.plotSampledPulse(self.sequence2, count, self.laserTime, self.sampleTime, self.figcount)
            self.figcount = self.figcount + 1


    def delayChannel(self, delay, channel = 1):
        if channel == 1:
            delaySequence = deque(self.sequence1)
        if channel == 2:
            delaySequence = deque(self.sequence2)
        if channel == 3:
            delaySequence = deque(self.clock_sequence)
        delaySamples = int((delay/1e12)/self.sampleTime)
        delaySequence.rotate(delaySamples)
        if channel == 1:
            self.sequence1 = list(delaySequence)
            print("Channel 1 delayed by ", delaySamples*round(self.sampleTime*1e12,2), " picoseconds, or ", delaySamples, "samples.")
        if channel == 2:
            self.sequence2 = list(delaySequence)
            print("Channel 2 delayed by ", delaySamples*round(self.sampleTime*1e12,2), " picoseconds, or ", delaySamples, "samples.")
        if channel == 3:
            self.clock_sequence = list(delaySequence)
            print("Clock delayed by ", delaySamples*round(self.sampleTime*1e12,2), " picoseconds, or ", delaySamples, "samples.")

    def plotSequencePortion(self, start, end):

        figure_size = (8, 6)
        fig = plt.figure(figsize=figure_size)
        border = 0.05
        spacing = 0.08
        height = (1 - 2 * border - spacing) / 2
        width = 1 - 2 * border

        ax2_size = [border, border,
                    width, height]
        ax1_size = [border, border + height + spacing,
                    width, height]

        ax1 = fig.add_axes(ax1_size)
        ax2 = fig.add_axes(ax2_size, sharex=ax1)

        plt.figure(self.figcount)
        self.figcount = self.figcount + 1
        x = [i*self.sampleTime*1e9 for i in range(start, end)]
        initialTime = x[0]
        pulseTime = initialTime  ###### change!

        while pulseTime < x[-1]:
            ax1.axvline(pulseTime, color="#ffc4b5")
            ax2.axvline(pulseTime, color="#ffc4b5")
            pulseTime = pulseTime + self.laserTime*1e9


        y1 = self.sequence1[start:end]
        y2 = self.sequence2[start: end]
        ax1.plot(x,y1)
        ax2.plot(x, y2)
        ax1.set_title("Channel 1")
        ax2.set_title("Channel 2")
        ax1.set_ylim(top=1.2)
        ax2.set_ylim(top=1.2)