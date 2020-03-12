import numpy as np

from fractions import Fraction
import matplotlib.pyplot as plt
import math
import csv
import random

# Cast bytes to bytearray
mutable_bytes = bytearray(b'\x00\x0F')

# Bytearray allows modification
mutable_bytes[0] = 255
mutable_bytes.append(255)
print(mutable_bytes)

# Cast bytearray back to bytes
immutable_bytes = bytes(mutable_bytes)
print(immutable_bytes)


#########################################################
#pulse time, cycles, rep rate, etc.





#AWG params
sampleRate = 90 #Gsps
MAXsampleDepth = 512 #KSamples

#System params
laserRate = 5 #GHz
amplitude = .7
activeTime = 12 #ns

MAXtcspcCycleRate = 8 #MHz

filesave = False




sampleTime = 1/(sampleRate*1e9)
laserTime = 1/(laserRate*1e9)
hightime = laserTime*.1


sigma = laserTime*.2
MINcycleTime = 1/(MAXtcspcCycleRate*1e6)
MINsamples_per_cycle = ((sampleRate*1e9)//(MAXtcspcCycleRate*1e6)) - 1


F = Fraction(int(laserRate*100), int(sampleRate*100))
mult = (MINsamples_per_cycle)//F.denominator
cycleDepth = (mult + 1)*F.denominator
cycleTime = cycleDepth*sampleTime
pulses_per_cycle = cycleTime/laserTime

###### should throw an error if not an int
cycles_per_sequence = int((MAXsampleDepth*1e3)//cycleDepth)

sequenceLength = int(cycleDepth*cycles_per_sequence)

figcount = 0





print("cycles per sequence: ", cycles_per_sequence)
print("sequence length: ", sequenceLength)
print("cycle time: ", cycleTime)
print("cycle depth is: ", cycleDepth)

print()
print()

print(F)


print("laser pulses per cycle: ", pulses_per_cycle)
print("TCSPC cycle rate: ", 1/(cycleTime*1e6), "MHz")


#how many pulses is about a microsecond or a little more?

#print("There are ", round(pulses_per_sequence,6), "pulses per sequence")
print("each laser pulse is", round(laserTime/sampleTime,6), "samples")
print("pulse time: ", laserTime*1e12, "ps")
print("number of pulses that fit in", activeTime ,"ns: ", 1e-9*activeTime//laserTime)




#########################################################
#shape and width of pulse



def gaussian(sigma, center, time):
    # sigma = pulse_width / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-(time - center) ** 2 / (2 * sigma ** 2))



class datagen():
    def __init__(self, max_bins, total_cycles):
        self.max_bins = max_bins
        self.total_cycles = int(total_cycles)
        self.bins = 2**math.floor(math.log2(max_bins))
        self.peaks = 16
        self.pulses_per_value = int(self.total_cycles//self.peaks) #2 for 32-ary alphabet
        self.datatype = "NAN"
        print("data bins available: ", self.bins)


    def testPattern1(self):
        self.datatype = "testpattern1"
        self.prob = [0]*16
        for i in range(0,8):
            self.prob[i] = 8-i
        for i in range(8,16):
            self.prob[i] = i - 7
        for item in self.prob: item = item/sum(self.prob)
        self.prob = [math.floor(self.bins*(i/sum(self.prob))) for i in self.prob]
        print("prob: ", self.prob)

        self.left_over = self.total_cycles - sum(self.prob)
        print("summed prob: ", sum(self.prob))
        print("left over: ", self.left_over)

        self.prob7 = self.left_over//2
        self.prob[8] = int(self.left_over - self.prob7)
        self.prob[7] = int(self.prob7)


        print("prob: ", self.prob)
        print("sum of prob", sum(self.prob))
        print()


        for item in self.prob:
            print(item)


        plt.figure(1)
        plt.bar(range(0,16),self.prob)

    def testPattern2(self):
        global figcount
        self.datatype = "testpattern2"
        self.prob = [0] * self.bins
        bandleft = 1
        bandright = 1
        cyclecapacity = self.total_cycles

        while cyclecapacity > 0:
            #print(bandleft)
            #print(bandright)
            #print("cyclecapacity is: ", cyclecapacity)
            for i in range(self.bins - 1, (self.bins-1) - bandright, -1):
                if cyclecapacity > 0:
                    self.prob[i] = self.prob[i] + 1
                    cyclecapacity = cyclecapacity - 1
                else:
                    break

            for i in range(0, bandleft):
                if cyclecapacity > 0:
                    self.prob[i] = self.prob[i] + 1
                    cyclecapacity = cyclecapacity - 1
                else:
                    break

            bandleft = bandleft + 1
            bandright = bandright + 1

        plt.figure(figcount)
        figcount = figcount + 1
        plt.bar(range(0, len(self.prob)), self.prob)
        print(self.prob)

    def testPulse(self, bin):
        self.datatype = "testpulse_bin_" + str(bin)
        self.prob = [0] * self.bins
        self.prob[bin] = self.total_cycles

    def generatePulseSequence(self):
        # generates list of numbers that correspond to PPM slots. For example,
        # if the PPM alphabet is 32 slots long, then pulseSequnces will hold values ranging
        # from 0 to 31.
        self.pulseSequence = [0]*self.total_cycles
        i = 0
        for index in range(len(self.prob)):
            for repeat in range(self.prob[index]):
                self.pulseSequence[i] = index
                i = i+1
        random.Random(4).shuffle(self.pulseSequence)
        print("this is the pulseSequence: ", self.pulseSequence)

        return self.pulseSequence





#in the future I might want to have different hightimes for pulses that are in the center of samples. Or vice versa.
#then every pulse should have its own object. So that a special method could be written that
class pulse():
    def __init__(self, hightime, risetime, center, sampleTime, amplitude):
        self.hightime = hightime
        self.risetime = risetime
        self.center = center
        self.sampleTime = sampleTime
        self.amplitude = amplitude
        self.left_center = 0
        self.right_center = 0
        self.nearCenterSampleIndex = 0
                                                    ##########global?
        self.nearCenterSampleIndex = self.center // sampleTime
        self.index_left = int(self.nearCenterSampleIndex)
        self.index_right = int(self.nearCenterSampleIndex + 1)

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

        # this should go inside the pulse class
        # write the pulse to the list iterating out from near the center
        while minLeftAmplitude > .002 or minRightAmplitude > .002:
            # print("location of pulse between indexes:", pulse_time/sampleTime)
            sequence_[self.index_left] = self.writePulseGaussian(self.index_left * sampleTime)
            if sequence_[self.index_left] <= minLeftAmplitude:
                minLeftAmplitude = sequence_[self.index_left]
            self.index_left = self.index_left - 1
            # print("self.index_left: ", self.index_left)

            sequence_[self.index_right] = self.writePulseGaussian(self.index_right * sampleTime)
            if sequence_[self.index_right] <= minRightAmplitude:
                minRightAmplitude = sequence_[self.index_right]
            self.index_right = self.index_right + 1
            # print("self.index_right: ", self.index_right)

            req = req + 1
            # print("cycle")
            if req > 1000:
                print("cycle error")
                exit(1)

    def plotSampledPulse(self, sequence_, count_, slot_):
        global figcount

        #be aware laserTime and sampleTime are still global
        x = []
        y = []
        for i in range(self.index_left - int(2 * laserTime // sampleTime),
                       self.index_right + int(2 * laserTime // sampleTime)): x.append(i * sampleTime * 1e9)

        for i in range(self.index_left - int(2 * laserTime // sampleTime),
                       self.index_right + int(2 * laserTime // sampleTime)): y.append(sequence_[i])


        plt.figure(figcount)
        figcount = figcount + 1
        plt.axvline(self.left_center * 1e9, color="#d1d1d1")
        plt.axvline(self.right_center * 1e9, color="#d1d1d1")

        plt.plot(x, y, color="#4189c4")

        plt.axvline(self.center * 1e9, color="#ffc4b5")
        plt.axvline((self.center + laserTime) * 1e9, color="#ffc4b5")
        plt.axvline((self.center - laserTime) * 1e9, color="#ffc4b5")
        plt.axvline((self.center + 2 * laserTime) * 1e9, color="#ffc4b5")
        plt.axvline((self.center - 2 * laserTime) * 1e9, color="#ffc4b5")
        # print("pulse laser time: ", pulse_time*1e9)
        # print("pulse laser time: ", pulse_time*1e9)
        # print("pulse laser time plus 1: ", pulse_time*1e9)

        plt.xlabel("time (ns)")
        plt.ylabel("Amplitude")
        plt.title("pulse in cycle number " + str(count_ + 1) +
            " in slot " + str(slot_))
        plt.show()




class PulseSequence():
    def __init__(self, cycles_per_sequence_, pulseSequence_, sampleTime_, padding_ = 0):
        self.pulseList = [] #this will be a list of pulse objects, each with their own center, hightime, risetime, etc. specified.
        self.times = [0] * cycles_per_sequence_
        self.pulseSequence = pulseSequence_
        self.padding = padding_
        self.sampleTime = sampleTime_

        self.sequence = [0]*sequenceLength
        self.zerosequence = [0]*sequenceLength

    def generateTimesList(self):
        # generate at list of universal AWG pulse times from the start time of each cycle
        # these AWG pulse times align with pulses of the laser
        self.cycleStart = self.padding
        for i in range(len(self.times)):
            self.times[i] = self.cycleStart + pulseSequence[i] * laserTime
            self.cycleStart = self.cycleStart + cycleTime
        print(self.times)

    def genDefaultPulseObjects(self, hightime_, sigma_, amplitude_):
        for i in range(len(self.pulseSequence)):
            self.pulseList.append(pulse(hightime_, sigma_, self.times[i], self.sampleTime, amplitude_))
        #print(self.pulseList)


    def writePulses(self):
        for pulseObject in self.pulseList:
            pulseObject.writeSampledPulse(self.sequence)

    def plotSomePulses(self, number):
        for count, pulseObject in enumerate(self.pulseList[0:number]):
            pulseObject.plotSampledPulse(self.sequence, count, self.pulseSequence[count])



data = datagen(1e-9*activeTime//laserTime, cycles_per_sequence)
data.testPattern2()
#data.testPulse(1) #for putting a single pulse in a specified slot in all cycles.
pulseSequence = data.generatePulseSequence()


PulseSeq = PulseSequence(cycles_per_sequence, pulseSequence, sampleTime, 0)
PulseSeq.generateTimesList()
PulseSeq.genDefaultPulseObjects(hightime, sigma, amplitude)
PulseSeq.writePulses()
PulseSeq.plotSomePulses(3)



#will want to pad the beginning with some zeros and move them to the very end of the sequence.
    #why??
#laser_pulse = 97 # this is bin 25 of, say, 256 in the PPM
#pulse_time = laser_pulse*laserTime
#p1 = pulse(hightime, sigma, pulse_time, sampleTime, amplitude)




file_name = "AWG" + "_LRate_" + str(laserRate) + "GHz" + "_SRate_" + str(
    sampleRate) + "GHz" + "_dBins_" + str(data.bins) + "_TCycles_" + str(
    cycles_per_sequence) + "_SLength_" + str(sequenceLength) + "_dtype_" \
            + str(data.datatype)
if filesave:
    with open(file_name + '.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(zip(PulseSeq.sequence, PulseSeq.zerosequence))










