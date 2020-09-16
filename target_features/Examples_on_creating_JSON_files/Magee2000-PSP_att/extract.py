import os
import pylab
import csv
import collections
import math
import numpy
import scipy
from itertools import cycle
import matplotlib.pyplot as plt
import json
import collections

amp_soma_vec_fig1 = []
amp_dend_vec_fig1 = []
dist_vec_fig1 = []
dist2_vec_fig1 = []

data = "Magee2000_Fig1_soma.csv"
csvreader = csv.reader(open(data, "r"))
for dist, amp in csvreader:
    dist_vec_fig1.append(float(dist))
    amp_soma_vec_fig1.append(float(amp))

dist_vec_fig1 = numpy.array(dist_vec_fig1)
amp_soma_vec_fig1 = numpy.array(amp_soma_vec_fig1)

data = "Magee2000_Fig1_dend.csv"
csvreader = csv.reader(open(data, "r"))
for dist, amp in csvreader:
    dist2_vec_fig1.append(float(dist))
    amp_dend_vec_fig1.append(float(amp))

dist2_vec_fig1 = numpy.array(dist2_vec_fig1)
amp_dend_vec_fig1 = numpy.array(amp_dend_vec_fig1)

plt.subplot(2,2,1)
plt.plot(dist_vec_fig1, amp_dend_vec_fig1, 'bo', dist_vec_fig1, amp_soma_vec_fig1, 'ro')
plt.subplot(2,2,2)
#plt.plot(dist_vec, 'b*', dist2_vec, 'r*')
plt.plot(dist_vec_fig1, amp_soma_vec_fig1/amp_dend_vec_fig1, 'bo')

output_file = open('Magee2000_Fig1_att.csv', "a")
for dist, att in zip(dist_vec_fig1, amp_soma_vec_fig1/amp_dend_vec_fig1):
    output_file.write("%f,%f\n" % (dist, att))
output_file.close()


amp_soma_vec_fig2 = []
amp_dend_vec_fig2 = []
dist_vec_fig2 = []
dist2_vec_fig2 = []

data = "Magee2000_Fig2_soma.csv"
csvreader = csv.reader(open(data, "r"))
for dist, amp in csvreader:
    dist_vec_fig2.append(float(dist))
    amp_soma_vec_fig2.append(float(amp))

dist_vec_fig2 = numpy.array(dist_vec_fig2)
amp_soma_vec_fig2 = numpy.array(amp_soma_vec_fig2)

data = "Magee2000_Fig2_dend.csv"
csvreader = csv.reader(open(data, "r"))
for dist, amp in csvreader:
    dist2_vec_fig2.append(float(dist))
    amp_dend_vec_fig2.append(float(amp))

dist2_vec_fig2 = numpy.array(dist2_vec_fig2)
amp_dend_vec_fig2 = numpy.array(amp_dend_vec_fig2)

plt.subplot(2,2,3)
plt.plot(dist_vec_fig2, amp_dend_vec_fig2, 'bo', dist_vec_fig2, amp_soma_vec_fig2, 'ro')
plt.subplot(2,2,4)
plt.plot(dist_vec_fig2, amp_soma_vec_fig2/amp_dend_vec_fig2, 'ro')

output_file = open('Magee2000_Fig2_att.csv', "a")
for dist, att in zip(dist_vec_fig2, amp_soma_vec_fig2/amp_dend_vec_fig2):
    output_file.write("%f,%f\n" % (dist, att))
output_file.close()

target_distances = [100,200,300]
tolerance = 50.0

plt.figure()
plt.plot(dist_vec_fig1, amp_soma_vec_fig1/amp_dend_vec_fig1, 'ko', dist_vec_fig2, amp_soma_vec_fig2/amp_dend_vec_fig2, 'go')

features = collections.OrderedDict()
stim = collections.OrderedDict()

for i, target_dist in enumerate(target_distances):
    target_fig1 = numpy.array(amp_soma_vec_fig1/amp_dend_vec_fig1)[((dist_vec_fig1 >= target_dist - tolerance) & (dist_vec_fig1 <= target_dist + tolerance))]
    target_fig2 = numpy.array(amp_soma_vec_fig2/amp_dend_vec_fig2)[((dist_vec_fig2 >= target_dist - tolerance) & (dist_vec_fig2 <= target_dist + tolerance))]

    #target = target_fig1
    target = numpy.concatenate((target_fig1, target_fig2), axis=0)
    mean = numpy.mean(target)
    std = numpy.std(target)

    print mean, std

    features['mean_attenuation_soma/dend_'+str(target_dist)+'_um'] = str(mean)+' mV'
    features['std_attenuation_soma/dend_'+str(target_dist)+'_um'] = str(std)+' mV'


    plt.errorbar(target_dist, mean, yerr=std, marker='o')

file_name= './feat_PSP_attenuation_target_data.json'
json.dump(features, open(file_name, "wb"), indent=4)

stim['target_distances'] = target_distances
stim['type'] = 'double exponential synapse'
stim['amplitude'] = None
stim['tau_rise'] = 0.1
stim['tau_decay'] = 3.0

file_name= './stim_PSP_attenuation_test.json'
json.dump(stim, open(file_name, "wb"), indent=4)


plt.show()
