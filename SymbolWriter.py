import json

RootData = {}
SampleRate = {"val": 92.2375, "unit": "GHz"}
MAXsampleDepth = {"val": 256, "unit": "KSamples"}
laserRepRate = {"val": 4.7, "unit": "GHz"}
#ActiveTime = {"val": 120, "unit": "ns"}
#MAXsymbolRate = {"val": 4.6, "unit": "MHz"}
defaultAmplitude = {"val": 1, "unit": "None"}
deadPulses = {"val": 0, "unit":"Laser Pulses"}
RootData["ExpParams"] = {"SampleRate": SampleRate,
                                "MAXsampleDepth":MAXsampleDepth,
                                "laserRepRate": laserRepRate,
                                #"ActiveTime":ActiveTime,
                                #"MAXsymbolRate": MAXsymbolRate,
                                "defaultAmplitude": defaultAmplitude,
                                "deadPulses": deadPulses}

highTime = {"val": 0.3, "unit": "fraction of laser period"}
sigma = {"val": 0.15, "unit": "fraction of laser period"}
RootData["PulseParams"] = {"highTime": highTime, "sigma":sigma}

Data = {}

# or load it from a different file here...

#I don't think each symbol needs an amplitude... That's determined in the next script when symbols are interpreted as
# pusles on two or more channels
"""
amplitudes = ""
for character in Symbols:
    amplitude = 1
    if character is "A":
        #do something with amplitude
        amplitude = 1
    if character is "B":
        # do something with amplitude
        amplitude = 1
    #amplitudes.append(amplitude)
    amplitudes += str(amplitude)
#amplitudes = [str(i) for i in amplitudes]
RootData["Data"] = {"Symbols": Symbols, "Amplitudes": amplitudes}
"""
#RootData["Data"] = Symbols


with open("symbol_file.json", "w") as write_file:
    json.dump(RootData, write_file, indent = 4)


#for opening the file in python:
with open("symbol_file.json", "r") as readfile:
    dataloaded = json.load(readfile)

#example of reading the file
print(dataloaded["ExpParams"])
print(dataloaded["ExpParams"]["laserRepRate"]["val"])
#print(dataloaded["Data"])




