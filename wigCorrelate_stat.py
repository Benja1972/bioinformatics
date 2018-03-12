#!/usr/bin/env python



import os
import yaml
import subprocess

stream = file("config.yaml", 'r')
config = yaml.load(stream)
samples = config["samples"]

subprocess.call(['bash','-c', 'rm 05peak/wigCorrelation_stat.txt'])

for i in range(len(config["samples"])/2):
	inp1 = "05peak/"+samples[2*i]+"_FE.bw "
	inp2 = "05peak/"+samples[2*i+1]+"_FE.bw "
	bashCommand = "wigCorrelate "+ inp1 + inp2 + " >> 05peak/wigCorrelation_stat.txt"
	print bashCommand
	output=subprocess.check_output(['bash','-c', bashCommand])


subprocess.call(['bash','-c', 'rm 06broad_peak/wigCorrelation_stat.txt'])

for i in range(len(config["samples"])/2):
	inp1 = "06broad_peak/"+samples[2*i]+"_FE.bw "
	inp2 = "06broad_peak/"+samples[2*i+1]+"_FE.bw "
	bashCommand = "wigCorrelate "+ inp1 + inp2 + " >> 06broad_peak/wigCorrelation_stat.txt"
	print bashCommand
	output=subprocess.call(['bash','-c', bashCommand])


