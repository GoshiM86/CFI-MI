function [spike_train, states_binary] = binary_representation(spike_train, b_idle_threshold)

% binary_representation.m generates binary profile decomposing neural
% spiking activity into the states of quiescence '0' and state of firing '1' unfolded in time

% Input parameters:
% 1. spike_train: spike
% 2. b_idle_threshold - parameter needed to estimate periods of neural quiescence (state '0'); 
% (idle_threshold = b*mean(ISI_strema)); % b_idle_threshold corresponds to b

% Output parameters parameters:
% 1. spike train containing the same number of samples as array of states
% 2. states_binary: binary profile

%	This function is part of the SFI MI toolbox.

ISI = diff(spike_train); 
ISI = round(ISI.*100)./100;
idle_threshold = b_idle_threshold*mean(ISI); % defining idle threshold

states_binary = zeros(numel(ISI),1); % initilization of states
ind_ones = find(ISI < idle_threshold); %  state '1' of neural firing
states_binary(ind_ones,1) = ones(numel(ind_ones), 1); % binary states

spike_train = spike_train(1:end-1,:); % numel(spike_train) must be equal to numel(states_binary)


