''' This file runs spike sorting on a signal segment of a compressed and an uncompressed signal, 
and compares the results with spike times found on the intracellular signal. These times are
found in the spikeSorting-HC1 intracellular.py file. Both that and this file are based on code 
available at https://spikeinterface.readthedocs.io/en/latest/index.html. 
'''

#%% Imports
import numpy as np
import spikeinterface.full as si
from pathlib import Path
from probeinterface import Probe, ProbeGroup, generate_linear_probe, generate_tetrode
from probeinterface.plotting import plot_probe, plot_probe_group
from spikeinterface.curation import CurationSorting

# import internal function
import sys
sys.path.append("C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data")
from src.data.saveFiltered import saveFiltered


# set new names for output files each time
rnd = np.random.randint(1,1000)
outputfolder = 'output' + str(rnd)
outputfolder_c = 'output_compressed' + str(rnd)
wavefolder = "waveforms_folder_hc_prefiltered" + str(rnd)
wavefolder_c = "waveforms_folder_hc_c" + str(rnd)


# parameters related to dataset
M = 600000 # length of dataset
sampling_frequency = 10_000.0  # sampling freq. of dataset
dtype = "float64"  # datatype of dataset
num_channels = 2  # number of channels in dataset

# added noise
noisy = True
sigma = 6


# choose which compressed signal to use
# CR, SNDR = 
# CR, SNDR = 11, 14.1
# CR, SNDR = 32, 10.1
# CR, SNDR = 35, 11.3
# # CR, SNDR = 53, 12.2 
# CR, SNDR = 61, 7.2
CR, SNDR = 73, 13.8
# # CR, SNDR = 98, 6.6
# CR, SNDR = 120, 6.2
# CR, SNDR = 131, 7.9 #ch=2
# CR, SNDR = 213, 5 #ch=4
# CR, SNDR = 261, 9.7 
# CR, SNDR = 361, 3
# CR, SNDR = 399, 3.6
# CR, SNDR = 624, 2.1


# multiprocesssing
global_job_kwargs = dict(n_jobs=1, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)
#%% This code filters the data, so that the uncompressed signal is filtered similarly to the compressed
saveFiltered('hc1_d533101_extracellular', fs=10000, f_low=300, max_samples=M, channel_max=num_channels, f_high=4000, matlab=True)
print("Remember to run matToBin('hc1_d533101_extracellular_f') in matlab")
#Now run matToBin('hc1_d533101_extracellular_f') in matlab

#%% import data
base_folder = Path('C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/bin')
if noisy == False:
    bin_str = "vq-lbg/hc1_cr-" + str(CR) +  "_sndr-" + str(SNDR) + "_c.bin"
else: 
        bin_str = "vq-lbg/hc1_cr-" + str(CR) +  "_sndr-" + str(SNDR) + "_c_noisy" + str(sigma) + ".bin"
        # bin_str = "vq-lbg/hc1_cr-" + str(CR) +  "_sndr-" + str(SNDR) + "_c_noisy.bin" # if file has unspecified noise level

# Load data 
recording = si.read_binary(file_paths=base_folder / "hc1_d533101_extracellular_f.bin", sampling_frequency=sampling_frequency,
                        num_channels=num_channels, dtype=dtype)
recording_c = si.read_binary(file_paths=base_folder / bin_str, sampling_frequency=sampling_frequency,
                        num_channels=num_channels, dtype=dtype)

print(recording)
print(recording_c)

#%% Make tetrodes
if num_channels == 4:
    tetrode = generate_tetrode()
    tetrode.set_device_channel_indices(np.arange(num_channels))

    plot_probe(tetrode, with_contact_id=True)

    recording_p = recording.set_probe(tetrode)
    recording_p_c = recording_c.set_probe(tetrode)
elif num_channels == 2: 
    
    linear_probe = generate_linear_probe(num_elec=2, ypitch=20)
    plot_probe(linear_probe, with_contact_id=True)
    linear_probe.set_device_channel_indices([0,1])

    recording_p = recording.set_probe(linear_probe)
    recording_p_c = recording_c.set_probe(linear_probe)
else: 
    raise ValueError('Probe needs to be defined for the specified num channels')

#%% Visualize
print(recording_p)
print(recording_p_c) 

w_ts = si.plot_traces(recording_p, time_range=(0, 5))
w_ts = si.plot_traces(recording_p_c, time_range=(0, 5))

#%% Pre-process
# filter (has been done beforehand instead)
recording_f = si.bandpass_filter(recording_p, freq_min=300, freq_max=4000)
recording_f_c = si.bandpass_filter(recording_p_c, freq_min=300, freq_max=4000)
# remove common refrence (has not been used)
recording_cmr = si.common_reference(recording_p, reference='global', operator='median')
recording_cmr_c = si.common_reference(recording_p_c, reference='global', operator='median')

# plot
# w = si.plot_timeseries(recording_f)
# w = si.plot_timeseries(recording_f_c)
# w = si.plot_timeseries(recording_cmr)
# w = si.plot_timeseries(recording_cmr_c)


# save recordings
recording_preprocessed = recording_p.save(format="binary")
recording_preprocessed_c = recording_p_c.save(format="binary")


#%% choose a sorter
print(si.available_sorters())
print(si.installed_sorters())

sorter = 'tridesclous'

print(si.get_default_sorter_params(sorter))

base_folder = Path('C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/clustering')
#%% Sort uncompressed data
sorting = si.run_sorter(sorter, recording_preprocessed, detect_threshold=4, output_folder=base_folder / outputfolder)  

#%% Sort compressed data
sorting_c = si.run_sorter(sorter, recording_preprocessed_c, detect_threshold=4, output_folder=base_folder / outputfolder_c)  

#%% remove excess spikes and extract waveforms
sorting_r = si.remove_excess_spikes(sorting, recording_preprocessed)
sorting_r_c = si.remove_excess_spikes(sorting_c, recording_preprocessed_c)
we = si.extract_waveforms(recording_preprocessed, sorting_r, folder=base_folder / wavefolder, overwrite=True, sparse=False, max_spikes_per_unit=2000, allow_unfiltered=True)
we_c = si.extract_waveforms(recording_preprocessed_c, sorting_r_c, folder=base_folder / wavefolder_c, overwrite=True, sparse=False, max_spikes_per_unit=2000, allow_unfiltered=True)
print(we)
print(we_c)

#%% Print the number of spikes found in the first unit
unit_id0 = sorting_r.unit_ids[0] # units are the different neurons we have found, here we look at one
wavefroms = we.get_waveforms(unit_id0) # we find the waveform and a template for the AP's of this neuron
print(wavefroms.shape)

template = we.get_template(unit_id0)
print(template.shape)

unit_id0 = sorting_r_c.unit_ids[0] # units are the different neurons we have found, here we look at one
wavefroms = we_c.get_waveforms(unit_id0) # we find the waveform and a template for the AP's of this neuron
print(wavefroms.shape)

template = we_c.get_template(unit_id0)
print(template.shape)


#%% plot waveforms and templates
w = si.plot_unit_waveforms(we, figsize=(16,16))
w = si.plot_unit_templates(we, figsize=(16,16))
si.plot_unit_waveforms_density_map(we, figsize=(16,16))

w = si.plot_unit_waveforms(we_c, figsize=(16,16))
w = si.plot_unit_templates(we_c, figsize=(16,16))
si.plot_unit_waveforms_density_map(we_c, figsize=(16,16))

# %% use compute parameters and make sorting viewer link (if needed)
# for uncompressed signal
# qm_params = si.get_default_qm_params()
# qm_params["presence_ratio"]["bin_duration_s"] = 1
# qm_params["amplitude_cutoff"]["num_histogram_bins"] = 5
# qm_params["drift"]["interval_s"] = 2
# qm_params["drift"]["min_spikes_per_interval"] = 2
# qm = si.compute_quality_metrics(we, qm_params=qm_params)
# print(qm)

# w1 = si.plot_quality_metrics(we, display=False, backend="sortingview")

# amplitudes         = si.compute_spike_amplitudes(   we)
# unit_locations     = si.compute_unit_locations(     we)
# spike_locations    = si.compute_spike_locations(    we)
# correlograms, bins = si.compute_correlograms(       we)
# similarity         = si.compute_template_similarity(we)
# w2 = si.plot_sorting_summary(we, display=False, curation=True, backend="sortingview")
#%% Merge units (if needed)

# merge compressed units
# cs = CurationSorting(parent_sorting=sorting_r_c)
# # cs.merge(units_to_merge=[0, 2])
# # cs.merge(units_to_merge=[1, 3])
# cs.merge(units_to_merge=np.arange(len(sorting_r_c.unit_ids)))
# sorting_r_c = cs.sorting

# merge uncompressed units
# cs = CurationSorting(parent_sorting=sorting_r)
# cs.merge(units_to_merge=[1, 2])
# sorting_r = cs.sorting

#%% Get compressed, uncompressed, and intracellular spike times
spike_times = sorting_r.to_spike_vector()['sample_index']
spike_class = sorting_r.to_spike_vector()['unit_index']

spike_times_c = sorting_r_c.to_spike_vector()['sample_index']
spike_class_c = sorting_r_c.to_spike_vector()['unit_index']


base_folder = Path('C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/npy')
if M == 240000:
    spike_times_true = np.load(  base_folder / 'spike_times_true_240000.npy')
    spike_times_true_c = np.load(base_folder / 'spike_times_true_240000.npy')
elif M == 600000:
    spike_times_true = np.load(  base_folder / 'spike_times_true_600000.npy')
    spike_times_true_c = np.load(base_folder / 'spike_times_true_600000.npy')
elif M == 1200000: 
    spike_times_true = np.load(  base_folder / 'spike_times_true.npy')
    spike_times_true_c = np.load(base_folder / 'spike_times_true.npy')
else: 
    raise ValueError('Intracellular spike times not found for this specified M; use ReadBinary.m, then run spikeSorting-HC1')

#%% Calculate TP, FP, and FN values. using a specified delta_t
delta = 10

# for uncompressed signal
true_positives = []
wrong_unit = []
false_positives = []

# for each extracellular spike find closest intracellular spike. Look at sample difference, if < delta; register as TP
for i in range(len(spike_times)):
    #find closest spike
    delta_spike = np.min(np.abs(spike_times[i] - spike_times_true))
    j = np.argmin(np.abs(spike_times[i] - spike_times_true))

    if delta_spike <= delta: #  check if spike close enough
        true_positives.append(spike_class[i]) # save as TP and register class
        spike_times_true = np.delete(spike_times_true, j)
    else: 
        false_positives.append(spike_class[i])


# print(f'True positives: {len(true_positives)}')
# print(f'False positives: {len(false_positives)}')
# print(f'Un detected: {len(spike_times_true)}')

# for compressed signal (same as for uncompressed)
true_positives_c = []
wrong_unit_c = []
false_positives_c = []
for i in range(len(spike_times_c)):
    #find closest spike
    delta_spike = np.min(np.abs(spike_times_c[i] - spike_times_true_c))
    j = np.argmin(np.abs(spike_times_c[i] - spike_times_true_c))

    if delta_spike <= delta: #  check if spike close enough
        true_positives_c.append(spike_class_c[i])
        spike_times_true_c = np.delete(spike_times_true_c, j)
    else: 
        false_positives_c.append(spike_class_c[i])

# print(f'True positives: {len(true_positives_c)}')
# print(f'False positives: {len(false_positives_c)}')
# print(f'Un detected: {len(spike_times_true_c)}')

#%%
# find unit with best match for signal and use this
true_positives = np.array(true_positives)
false_positives = np.array(false_positives)
unit_idx_focus = np.argmax( np.array([len(true_positives[ true_positives  == i]) for i in range(len(sorting_r.unit_ids))]) )
unit_idx_other = np.delete(sorting_r.unit_ids, unit_idx_focus)

print(f'TP: {len(  true_positives[ true_positives == unit_idx_focus ])}') # true spikes of correct class
print(f'FP: {len( false_positives[false_positives == unit_idx_focus])}') # spikes of that class that did not have intracellular match
print(f'FN: {len(spike_times_true) + len([tp for tp in true_positives if tp in unit_idx_other])}') # intracellular spikes not found, and true spikes wrongly classified

# find unit with best match for compressed signal and use this
true_positives_c = np.array(true_positives_c)
false_positives_c = np.array(false_positives_c)
unit_idx_focus_c = np.argmax( np.array([len(true_positives_c[ true_positives_c  == i]) for i in range(len(sorting_r_c.unit_ids))]) )
unit_idx_other_c = np.delete(sorting_r_c.unit_ids, unit_idx_focus_c)

print(f'TP compressed: {len(  true_positives_c[ true_positives_c == unit_idx_focus_c ])}')
print(f'FP compressed: {len( false_positives_c[false_positives_c == unit_idx_focus_c])}')
print(f'FN compressed: {len(spike_times_true_c) + len([tp for tp in true_positives_c if tp in unit_idx_other_c])}')

# Confusion matrix for comparison is implemented in SpikeInterface, and can be added if desired