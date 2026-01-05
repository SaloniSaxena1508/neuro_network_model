#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 16:58:43 2025

@author: salonisaxena
This code imports the firing rates pre and post-stimulus and computes the modulation index for passive and active. 
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

main_folder = "unstruc_no_overlap"  # change to unstruc_overlap or struc_no_overlap
REPO_ROOT = Path.cwd()

file_path = REPO_ROOT/main_folder

stimID = np.loadtxt(file_path/'stimID.txt').astype(int)   # IDs of photostim neurons
contextID = np.loadtxt(file_path/'contextID.txt').astype(int)  # IDs of context neurons

ne1, ne2, ne_non, ni1, ni2, ni_non = 2040, 4080, 6800, 7160, 7520, 8000 # E1 IDs = 1:2041, E2 IDs = 2041:4081 etc
N_neurons = 8000
N = np.arange(1, N_neurons+1)

stim, context = 'yes', 'no' # stim = yes if photostim applied and no if only sound, context = yes for active and no for passive

""" stim modulation indices, active and passive """
if stim == 'yes' and context == 'yes':
    sound_post_active = np.loadtxt(file_path/'left_sound_active/rate_post.txt') * 1000
    stimsound_post_active = np.loadtxt(file_path/'left_stimsound_active/rate_post.txt') * 1000
    normalised_stim_active = (stimsound_post_active - sound_post_active)/(stimsound_post_active + sound_post_active)
    stim_active_no_denom = stimsound_post_active - sound_post_active
    """np.savetxt(file_path/'normalized_mi_stim_act.txt', normalised_stim_active)
    np.savetxt(file_path/'mi_stim_act_no_denom.txt', stim_active_no_denom)
    """

if stim=='yes' and context=='no':
    sound_post_passive = np.loadtxt(file_path/'left_sound_passive/rate_post.txt') * 1000 # post sound passive rates in Hz
    stimsound_post_passive = np.loadtxt(file_path/'left_stimsound_passive/rate_post.txt') * 1000
    normalised_stim_passive = (stimsound_post_passive - sound_post_passive)/(stimsound_post_passive + sound_post_passive)
    stim_passive_no_denom = stimsound_post_passive - sound_post_passive 
    """np.savetxt(file_path/'normalized_mi_stim_pass.txt', normalised_stim_passive)
    np.savetxt(file_path/'mi_stim_pass_no_denom.txt', stim_passive_no_denom)
    """


""" sound modulation index, passive and active"""
if stim == 'no' and context == 'no':
    sound_post_passive = np.loadtxt(file_path/'left_sound_passive/rate_post.txt') * 1000
    sound_pre_passive = np.loadtxt(file_path/'left_sound_passive/rate_pre.txt') * 1000
    normalised_sound_passive = (sound_post_passive - sound_pre_passive)/(sound_post_passive + sound_pre_passive)
    sound_passive_no_denom = sound_post_passive - sound_pre_passive 
    """np.savetxt(file_path/'normalized_mi_sound_pass.txt', normalised_sound_passive)
    np.savetxt(file_path/'mi_sound_pass_no_denom.txt', sound_passive_no_denom)
    """
    
if stim == 'no' and context =='yes':
    sound_post_active = np.loadtxt(file_path/'left_sound_active/rate_post.txt') * 1000
    sound_pre_active = np.loadtxt(file_path/'left_sound_active/rate_pre.txt') * 1000
    normalised_sound_active = (sound_post_active - sound_pre_active)/(sound_post_active + sound_pre_active)
    sound_active_no_denom = sound_post_active - sound_pre_active
    """np.savetxt(file_path/'normalized_mi_sound_act.txt', normalised_sound_active)
    np.savetxt(file_path/'mi_sound_act_no_denom.txt', sound_active_no_denom)
    """
