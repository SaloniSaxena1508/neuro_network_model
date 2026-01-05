#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 14:29:16 2025

@author: salonisaxena
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os
import re
import time
from scipy.ndimage import gaussian_filter
from pathlib import Path

main_folder = "unstruc_no_overlap"
context_input_folder = "left_stimsound_active" #photostim or no photostim, active or passive
REPO_ROOT = Path.cwd()

file_path = REPO_ROOT/main_folder/context_input_folder

stimID = np.loadtxt(REPO_ROOT/main_folder/'stimID.txt').astype(int) 
contextID = np.loadtxt(REPO_ROOT/main_folder/'contextID.txt').astype(int)

ne1, ne2, ne_non, ni1, ni2, ni_non = 2040, 4080, 6800, 7160, 7520, 8000
N_neurons = 8000
N = np.arange(1, N_neurons+1)

def rate_v_time(spike_train): # spike train is an array of spike times and spiking neuron IDs
    """ average firing rate of each population versus time """
    s1_time, s1_id = spike_train[0], spike_train[1]
    bins_T = np.arange(0, T+db, db)
    epop1, epop2, e_rest, ipop1, ipop2, i_rest = [], [], [], [], [], []
    for i in range(len(bins_T)-1):
        tmin, tmax = bins_T[i], bins_T[i+1]
        t_ind = np.argwhere((s1_time >= tmin) & (s1_time < tmax))
        nu_ind = s1_id[t_ind]
        epop1.append(len(np.argwhere(nu_ind <= ne1)))
        epop2.append(len( np.argwhere((nu_ind - (ne1+1)) * (nu_ind - ne2) <= 0)))
        e_rest.append(len(np.argwhere((nu_ind - (ne2+1)) * (nu_ind - ne_non) <= 0)))
        ipop1.append(len(np.argwhere((nu_ind - (ne_non+1)) * (nu_ind - ni1) <= 0)))
        ipop2.append(len(np.argwhere((nu_ind - (ni1+1)) * (nu_ind - ni2) <= 0)))
        i_rest.append(len(np.argwhere((nu_ind - (ni2+1)) * (nu_ind - ni_non) <= 0)))
    return np.asarray(epop1)/db/2040, np.asarray(epop2)/db/2040, np.asarray(e_rest)/db/2720, np.asarray(ipop1)/db/360, np.asarray(ipop2)/db/360, np.asarray(i_rest)/db/480
    
    
def time_avg_pre(spike_train, ID):  
    spike_time, spike_id = spike_train[0], spike_train[1]
    pre = np.argwhere((spike_time >= 300) & (spike_time<T_on)) # pre-stimulus time period. discard transients upto 500
    spike_pre = spike_id[pre.astype(int)]   #neurons spike in pre interval
    
    rate_pre = np.zeros(len(ID))  #vector to store pre-stimulus firing rate of each neuron
    for i in range(0, len(ID)):
       rate_pre[i] = len(np.argwhere(spike_pre == ID[i]))/(T_on-300)       
    return rate_pre      

def time_avg_post(spike_train, ID):
    spike_time, spike_id = spike_train[0], spike_train[1]
    post = np.argwhere((spike_time <= T_off) & (spike_time >= T_on))
    spike_post = spike_id[post.astype(int)]
    
    rate_post = np.zeros(len(ID))
    for i in range(0, len(ID)):
       rate_post[i] = len(np.argwhere(spike_post == ID[i]))/(T_off-T_on)       
    return rate_post    


#%%INDICES OF EACH SUBPOPULATION, DEFINE PATH ADDRESSES 
db = 1
T = 4000 # 1500
dt = 0.1
T_on, T_off = 2000, 2300  #1000, 1300   #T_off=2300 for stim
num_trials = 100

N_neurons = 8000
N = np.arange(1, N_neurons+1)  #8001 range of neuron ids
ne1, ne2, ne_non, ni1, ni2, ni_non = 2040, 4080, 6800, 7160, 7520, 8000

E1 = np.arange(1, ne1+1) #epop1 1-2040
E2 = np.arange(ne1+1, ne2+1) #epop2  2041-4080
Erest = np.arange(ne2+1, ne_non+1) #epop nonsel 4081-6800
I1 = np.arange(ne_non+1, ni1+1) #ipop1 6801-7160
I2 = np.arange(ni1+1, ni2+1)  #ipop2 7161-7520
Irest = np.arange(ni2+1, N_neurons+1) #ipop nonsel 7520-8000

nonstimID = N[~np.isin(N, stimID)]  #unstimulated neurons

E1stim = stimID[np.argwhere(stimID < ne1+1)]   #stimulated neurons from epop1
E2stim = stimID[np.argwhere((stimID >= ne1+1) & (stimID < ne2+1) )] #stimulated neurons from epop2 
Ereststim = stimID[np.argwhere( (stimID >= ne2+1) & (stimID < ne_non+1) )] #stimulated neurons from enon-sel
I1stim = stimID[np.argwhere( (stimID >= ne_non+1) & (stimID < ni1+1) )]
I2stim = stimID[np.argwhere( (stimID >= ni1+1) & (stimID < ni2+1) )]
Ireststim = stimID[np.argwhere( (stimID >= ni2+1) & (stimID < N_neurons+1) )]

E1nonstim = nonstimID[np.argwhere(nonstimID < ne1+1)]   #stimulated neurons from epop1
E2nonstim = nonstimID[np.argwhere((nonstimID >= ne1+1) & (nonstimID < ne2+1) )] #stimulated neurons from epop2 
Erestnonstim = nonstimID[np.argwhere( (nonstimID >= ne2+1) & (nonstimID < ne_non+1) )] #stimulated neurons from enon-sel
I1nonstim = nonstimID[np.argwhere( (nonstimID >= ne_non+1) & (nonstimID < ni1+1) )]
I2nonstim = nonstimID[np.argwhere( (nonstimID >= ni1+1) & (nonstimID < ni2+1) )]
Irestnonstim = nonstimID[np.argwhere( (nonstimID >= ni2+1) & (nonstimID < N_neurons+1) <= 0)]

# context current
noncontID = N[~np.isin(N, contextID)]  #unstimulated neurons

E1_cont = contextID[np.argwhere(contextID < ne1+1)]   #stimulated neurons from epop1
E2_cont = contextID[np.argwhere((contextID >= ne1+1) & (contextID < ne2+1) )] #stimulated neurons from epop2 
Erest_cont = contextID[np.argwhere( (contextID >= ne2+1) & (contextID < ne_non+1) )] #stimulated neurons from enon-sel
I1_cont = contextID[np.argwhere( (contextID >= ne_non+1) & (contextID < ni1+1) )]
I2_cont = contextID[np.argwhere( (contextID >= ni1+1) & (contextID < ni2+1) )]
Irest_cont = contextID[np.argwhere( (contextID >= ni2+1) & (contextID < N_neurons+1) )]

E1_noncont = noncontID[np.argwhere(noncontID < ne1+1)]   #stimulated neurons from epop1
E2_noncont = noncontID[np.argwhere((noncontID >= ne1+1) & (noncontID < ne2+1) )] #stimulated neurons from epop2 
Erest_noncont = noncontID[np.argwhere( (noncontID >= ne2+1) & (noncontID < ne_non+1) )] #stimulated neurons from enon-sel
I1_noncont = noncontID[np.argwhere( (noncontID >= ne_non+1) & (noncontID < ni1+1) )]
I2_noncont = noncontID[np.argwhere( (noncontID >= ni1+1) & (noncontID < ni2+1) )]
Irest_noncont = noncontID[np.argwhere( (noncontID >= ni2+1) & (noncontID < N_neurons+1) <= 0)]

"""overlap_enon = np.loadtxt('/Users/salonisaxena/Downloads/code/model7/overlap_enon.txt').astype(int)
overlap_inon = np.loadtxt('/Users/salonisaxena/Downloads/code/model7/overlap_inon.txt').astype(int)
"""
""" RATES FOR INDIVIDUAL NEURONS PRE AND POST """
#np.random.seed(42)
ID = np.arange(1, 8001)
#ID = np.random.choice(nonstim, 1200)
folder_path = file_path

if context_input_folder == "left_sound_passive":
    pattern = r"2clusterLNP_LeftSoundI_ID\d+\.mat"  #SOUND, PASSIVE
    
if context_input_folder == "left_sound_active":
    pattern = r"2clusterLNP_Ibias_LeftSoundI_ID\d+\.mat"  #SOUND, ACTIVE
 
if context_input_folder == "left_stimsound_passive":
    pattern = r"2clusterLNP_LeftSSI_ID\d+\.mat"  #STIM+SOUND, PASSIVE

if context_input_folder == "left_stimsound_active":
    pattern = r"2clusterLNP_Ibias_LeftSSI_ID\d+\.mat"  #STIM+SOUND, ACTIVE

    
mat_files = [f for f in os.listdir(folder_path) if re.match(pattern, f)]

response_pre = np.zeros(len(ID))
for file in sorted(mat_files, key=lambda x: int(re.search(r'\d+', x).group())):  #average over trials
    file_path = os.path.join(folder_path, file)
    data = scipy.io.loadmat(file_path)  # Load .mat file
    spike_train = data['s1']
    response_pre += time_avg_pre(spike_train, ID)
   
response_pre = response_pre/num_trials
np.savetxt(file_path/'rate_pre.txt', response_pre)

response_post = np.zeros(len(ID))

for file in sorted(mat_files, key=lambda x: int(re.search(r'\d+', x).group())):  #average over trials
    file_path = os.path.join(folder_path, file)
    data = scipy.io.loadmat(file_path)  # Load .mat file
    spike_train = data['s1']
    response_post += time_avg_post(spike_train, ID)
 
response_post = response_post/num_trials   
np.savetxt(file_path/'rate_post.txt', response_post)

epop1_resp, epop2_resp, e_rest_resp, ipop1_resp, ipop2_resp, i_rest_resp = np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T)

for file in sorted(mat_files, key=lambda x: int(re.search(r'\d+', x).group())):  #average over populations
    file_path = os.path.join(folder_path, file)
    data = scipy.io.loadmat(file_path)  # Load .mat file
    spike_train = data['s1']
    response = rate_v_time(spike_train)
    epop1_resp += response[0]
    epop2_resp += response[1]
    e_rest_resp += response[2]
    ipop1_resp += response[3]
    ipop2_resp += response[4]
    i_rest_resp += response[5]

epop1_resp = epop1_resp/num_trials
epop2_resp = epop2_resp/num_trials
e_rest_resp = e_rest_resp/num_trials
ipop1_resp = ipop1_resp/num_trials
ipop2_resp = ipop2_resp/num_trials
i_rest_resp = i_rest_resp/num_trials

np.savetxt(file_path/'epop1_resp.txt', epop1_resp)
np.savetxt(file_path/'epop2_resp.txt', epop2_resp)
np.savetxt(file_path/'e_rest_resp.txt', e_rest_resp)
np.savetxt(file_path/'ipop1_resp.txt', ipop1_resp)
np.savetxt(file_path/'ipop2_resp.txt', ipop2_resp) 
np.savetxt(file_path/'i_rest_resp.txt', i_rest_resp)

