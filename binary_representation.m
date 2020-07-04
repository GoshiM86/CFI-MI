function [spike_train, states_binary] = binary_representation(spike_train, b_idle_threshold)

% binary_representation.m generates binary profile of neural
% spiking activity by its decomposition into the states of quiescence '0' and state of firing '1' unfolded in time

% Input parameters:
% 1. spike_train: spiking activity
% 2. b_idle_threshold - parameter needed to be set for estimation periods of neural quiescence (state '0'); 
% (idle_threshold = b*mean(ISI_stream)); b_idle_threshold corresponds to b

% Output parameters parameters:
% 1. spike_train: spiking activity containing the same number of samples (timings) as array of states (total number of '1' and '0' states in time)
% 2. states_binary: binary profile of a input spiking activity

% This function is part of the CFI-MI toolbox.

%--------------------------------------------------------------------------

ISI = diff(spike_train); % inter-spike interval stream
ISI = round(ISI.*100)./100;
idle_threshold = b_idle_threshold*mean(ISI); % defining the idle threshold

states_binary = zeros(numel(ISI),1); % initilization of states
ind_ones = find(ISI < idle_threshold); %  states '1'
states_binary(ind_ones,1) = ones(numel(ind_ones), 1); % binary profile

spike_train = spike_train(1:end-1,:); % numel(spike_train) must be equal to numel(states_binary)


