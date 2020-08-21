import numpy as np
from fractions import Fraction
import matplotlib.pyplot as plt
import math
import csv
import random
import yaml
import json
#import numpy as np
from random import *
from CodeGeneratorUtil import *






with open("symbol_file.json", "r") as readfile:
    dat = json.load(readfile)

clockgen = True
filesave = True

# AWG params
sampleRate = dat["ExpParams"]["SampleRate"]["val"]
MAXsampleDepth = 1e3 * dat["ExpParams"]["MAXsampleDepth"]["val"]  # memory limit of AWG

# System params
laserRate = dat["ExpParams"]["laserRepRate"]["val"]  # GHz
amplitude = dat["ExpParams"]["defaultAmplitude"]["val"]
deadPulses = dat["ExpParams"]["deadPulses"]["val"]
#Data = dat["Data"]
sampleTime = 1 / (sampleRate * 1e9)
laserTime = 1 / (laserRate * 1e9)
hightime = laserTime * .2

sigma = laserTime * .1
#MINcycleTime = 1 / (MAXsymbolRate * 1e6)
#MINsamples_per_cycle = ((sampleRate * 1e9) // (MAXsymbolRate * 1e6)) - 1

F = Fraction(int(laserRate * 1000000), int(sampleRate * 1000000))
# numerator is number of laser pulses
# denominator is number of samples tha take same amount to time

usableData = 700
padding = [0,0,0,0,0]
dataLength = usableData + len(padding) #data with 5 laserpulses of padding at beginning

usableData1 = [randint(0,1) for i in range(usableData)]
usableData2 = [randint(0,1) for i in range(usableData)]
Data1 = padding + usableData1
Data2 = padding + usableData2
assert len(Data1) == len(Data2) == dataLength
print(Data1[0:20])
print(Data2[0:20])


dataLength = len(Data1)
totalPulses_base = dataLength*(1+deadPulses)
print("totalPulses_base is: ", totalPulses_base)

print("Dead pusles: ", deadPulses)
if totalPulses_base > (F.numerator/F.denominator)*MAXsampleDepth:
    #too much data for memory in AWG. Round down to nearest set.
    sets = int(MAXsampleDepth//F.denominator)
    totalPulses_f = int(F.numerator * sets)
    totalSamples_f = int(F.denominator * sets)
    usableSymbols = int(totalPulses_f // (1 + deadPulses))
    usablePulses = int(usableSymbols * (1 + deadPulses)) #includes dead pulses
    print("Too much data for AWG memory")
    print(usableSymbols, "symbols will be written of the ", dataLength, "symbols given")
else:
    #all data fits in one cycle of AWG
    print("All symbols fit in AWG memory")
    sets = int(totalPulses_base//F.numerator)
    totalPulses_f = int(F.numerator*sets)
    #print(totalPulses_f)
    #print(F.numerator*sets)
    if totalPulses_f == totalPulses_base:
        #the number of pulses is a multiple of F.numerator
        totalSamples_f = int(F.denominator * sets)
        usableSymbols = totalPulses_base//(1+deadPulses) #which should be the same as dataLength
        assert usableSymbols == dataLength
    else:
        # add a certain number of extra sets so that there is at least deadPulses or more left over at the end.
        additional_sets = int((1 + deadPulses)//F.numerator) + 1
        totalPulses_f = int(F.numerator * (sets + additional_sets))
        totalSamples_f = int(F.denominator * (sets + additional_sets))
        usableSymbols = int(totalPulses_f // (1 + deadPulses))
        assert usableSymbols >= dataLength

print("usable symbols: ", usableSymbols)
print("totalSamples_f", totalSamples_f)


Data1 = Data1[:usableSymbols]
Data2 = Data2[:usableSymbols]

figcount = 0
PulseSeq = PulseSequence(Data1, Data2, sampleTime, laserTime, deadPulses, totalSamples_f, usableSymbols, padding_ = 0)
PulseSeq.generateTimesList()
PulseSeq.genDefaultPulseObjects(hightime, sigma, amplitude, Data1, Data2)
PulseSeq.writePulses()
PulseSeq.plotSomePulses(17,20) # plot the 17th pulse through the 20th pulse
PulseSeq.plotSequencePortion(0,10000)

if clockgen:
    clock_sequence = [0] * totalSamples_f

    for i in range(totalSamples_f//10):  # just sets the beginning of the sequence to 1, rest is zero.
        clock_sequence[i] = 1




#will want to pad the beginning with some zeros and move them to the very end of the sequence.
    #why??
#laser_pulse = 97 # this is bin 25 of, say, 256 in the PPM
#pulse_time = laser_pulse*laserTime
#p1 = pulse(hightime, sigma, pulse_time, sampleTime, amplitude)




file_name = "AWG" + "_LRate_" + str(laserRate) + "GHz" + "_SRate_" + str(
    sampleRate) + "GHz" + "_deadBins_" + str(deadPulses) +  "_SLength_" + str(totalSamples_f)


if filesave:
    with open(file_name + '_CHANNEL_1.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(zip(PulseSeq.sequence1, PulseSeq.zerosequence))

    with open(file_name + '_CHANNEL_2.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(zip(PulseSeq.sequence2, PulseSeq.zerosequence))


    if clockgen:
        with open(file_name + '_CLOCK' + '.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(zip(clock_sequence, PulseSeq.zerosequence))


# where saving to the database would be added
'''
file_name = "AWG" + "_LRate_" + str(laserRate) + "GHz" + "_SRate_" + str(
    sampleRate) + "GHz" + "_deadBins_" + str(deadPulses) +  "_SLength_" + str(sequenceLength)

file_name_yaml = "AWG" + "_LRate_" + str(laserRate) + "GHz" + "_SRate_" + str(
    sampleRate) + '.yaml'


users = [{"LaserRate": laserRate, 'laserTime': laserTime, 'deadBins': deadPulses, 'CyclesPerSequence': cycles_per_sequence, 'Pulses per Cycle': pulses_per_cycle}]

with open(file_name_yaml, 'w') as f:
    data = yaml.dump(users, f)

'''





