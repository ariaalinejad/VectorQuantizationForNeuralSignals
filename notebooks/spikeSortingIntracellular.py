''' This file runs spike sorting on an intracellular signal and saves
a .npy array with the spike times. Code in this file is based on code 
available at https://spikeinterface.readthedocs.io/en/latest/index.html. 
'''

#%% Imports
import numpy as np
import spikeinterface.full as si
from pathlib import Path
from probeinterface import Probe, ProbeGroup
from probeinterface import generate_linear_probe
from probeinterface.plotting import plot_probe, plot_probe_group
from probeinterface import generate_tetrode

# make new output folder name each time
rnd = np.random.randint(1,1000)
outputfolder = 'output_intra' + str(rnd)
wavefolder = "waveforms_intra" + str(rnd)

# multiprocesssing
global_job_kwargs = dict(n_jobs=1, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)
#%% import data
M = 600000
sampling_frequency = 10_000.0 
num_channels = 1  
dtype = "float64" 

base_folder = Path('C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/matlab/src/read_hc1')

# Load data  
if M == 240000:
    recording = si.read_binary(file_paths=base_folder / 'hc1_d533101_intracellular_240000.bin', sampling_frequency=sampling_frequency,
                        num_channels=num_channels, dtype=dtype)
elif M == 600000:
    recording = si.read_binary(file_paths=base_folder / 'hc1_d533101_intracellular_600000.bin', sampling_frequency=sampling_frequency,
                        num_channels=num_channels, dtype=dtype)
else:
    raise ValueError('Binary not made for selected M, use ReadBinary.m')

print(recording)
#%% Make probe
linear_probe = generate_linear_probe(num_elec=1, ypitch=20)
plot_probe(linear_probe, with_contact_id=True)
linear_probe.set_device_channel_indices([0])

# set probe electrode
recording_p = recording.set_probe(linear_probe)

#%% Visualize
print(recording_p) 

w_ts = si.plot_traces(recording_p, time_range=(0, 0.25))

#%% Save recording 
# recording_f = si.bandpass_filter(recording_p, freq_min=300, freq_max=4000) # BP filter if desired
recording_preprocessed = recording_p.save(format="binary")
print(recording_preprocessed)

w = si.plot_timeseries(recording_preprocessed)

#%% choose a sorter
print(si.available_sorters())
print(si.installed_sorters())

sorter = 'tridesclous'

print(si.get_default_sorter_params(sorter))

#%% Run spike sorting
base_folder = Path('C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/clustering')
sorting = si.run_sorter(sorter, recording_preprocessed, detect_threshold=4, output_folder=base_folder / outputfolder)
print(sorting)


#%% extract waveforms
we = si.extract_waveforms(recording_preprocessed, sorting, folder=base_folder / wavefolder, overwrite=True, max_spikes_per_unit=2000, allow_unfiltered=True)
print(we)

#%% Plot spikes found for first unit, and template
unit_id0 = sorting.unit_ids[0] 
wavefroms = we.get_waveforms(unit_id0) 
print(wavefroms.shape)

template = we.get_template(unit_id0)
print(template.shape)


#%% plot waveforms
w = si.plot_unit_waveforms(we, figsize=(16,16))
w = si.plot_unit_templates(we, figsize=(16,16))
si.plot_unit_waveforms_density_map(we, figsize=(16,16))

#%% Print num spikes found in each unit
for i in range(len(sorting.unit_ids[:])):
    print(np.shape(we.get_waveforms( sorting.unit_ids[i])))


#%% Use peak detection instead of spike sorting
# from spikeinterface.sortingcomponents.peak_detection import detect_peaks
# import matplotlib.pyplot as plt

# noise_levels_int16 = si.get_noise_levels(recording_preprocessed, return_scaled=False)
# peaks = detect_peaks(recording_preprocessed)#, detect_threshold=10)

# from spikeinterface.sortingcomponents.peak_localization import localize_peaks

# peak_locations = localize_peaks(recording_preprocessed, peaks, method='center_of_mass')

# fs = recording_preprocessed.sampling_frequency
# fig, ax = plt.subplots(figsize=(10, 8))
# ax.scatter(peaks['sample_index'] / fs, peak_locations['y'], color='k', marker='.',  alpha=0.01)

#%% use sorting viewer to evaluate detected units
qm_params = si.get_default_qm_params()
qm_params["presence_ratio"]["bin_duration_s"] = 1
qm_params["amplitude_cutoff"]["num_histogram_bins"] = 5
qm_params["drift"]["interval_s"] = 2
qm_params["drift"]["min_spikes_per_interval"] = 2
qm = si.compute_quality_metrics(we, qm_params=qm_params)
print(qm)

amplitudes         = si.compute_spike_amplitudes(   we)
unit_locations     = si.compute_unit_locations(     we)
spike_locations    = si.compute_spike_locations(    we)
correlograms, bins = si.compute_correlograms(       we)
similarity         = si.compute_template_similarity(we)
w2 = si.plot_sorting_summary(we, display=False, curation=True, backend="sortingview")

# %% do chosen correction from sortingview
if M == 240000:
    uri = 'sha1://87f1bc503a6b69f529fa0a9768b4d8587afeb6bf'
elif M == 600000:
    uri = 'sha1://41d7b9d24de463f5fa1a49fed84f401f4989445a'
elif M == 1200000: # full signal
    uri = 'sha1://f8fa49c28beb589701e771f72c089b53290f59e2' 
else: 
    raise ValueError('Uri undefined for specified M, use sortingview link above')

# apply correction to sorting
sorting_curated_sv = si.apply_sortingview_curation(sorting, uri_or_json=uri)
print(sorting_curated_sv)
mask = sorting_curated_sv.get_property("accept")
print(mask)
sorting_curated_keep = sorting_curated_sv.select_units(sorting_curated_sv.unit_ids[mask])
print(sorting_curated_keep)

#%% extract waveforms from corrected sorting
we_2 = si.extract_waveforms(recording_preprocessed, sorting_curated_keep, folder="waveforms_folder_synth_currated_2", overwrite=True, max_spikes_per_unit=2000, allow_unfiltered=True)
print(we_2)

#%% Plot waveforms
w = si.plot_unit_waveforms(we_2, figsize=(16,16))
w = si.plot_unit_templates(we_2, figsize=(16,16))
si.plot_unit_waveforms_density_map(we_2, figsize=(16,16))

#%% get spike times
print(np.shape(sorting_curated_keep.to_spike_vector()['sample_index']))
spike_times_true = sorting_curated_keep.to_spike_vector()['sample_index']

#%% Save spike times as .npy
base_folder = Path('C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/npy')
if M == 240000:
    np.save(base_folder / 'spike_times_true_240000.npy', spike_times_true)  
elif M == 600000:
    np.save(base_folder / 'spike_times_true_600000.npy', spike_times_true)  
elif M ==1200000: # full signal
    np.save(base_folder / 'spike_times_true.npy', spike_times_true)  
else: 
    raise ValueError('Code for saving spike times not defined for given M')