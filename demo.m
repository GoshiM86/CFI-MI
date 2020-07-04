close all; clc;

% This implements the simulation of new proposed measure by G. Mijatovic, T. Loncar Turukalo, N.
% Bozanic and L. Faes: 'A measure of concurrent neural firing activity based on mutual information', 2020.

% CFI_MI index estimates the degree of concomitant firing between two neural units 
% based on a modified form of a mutual information (MI) applied on a coarse, binary representations 
% of firing activity (states of neural quiescence '0' and states of firing '1' unfolded in time).

% Script demo.m simulates a network of 1000 randomly coupled spiking neurons (see REF1)
% and estimates the CFI_MI between two selected cells with assessment of its 
% statistical significance (see REF2) and the corresponding visualization of spike trains and binary profiles. 

% This script uses next functions: binary_representation.m; function_CFI_MI.m; spiSeMe_surrogate_jodi.m as a part of the CFI-MI toolbox.

%--------------------------------------------------------------------------

%% Input parameters
b_idle_threshold = 3; % in order to estimate periods of neural quiescence (state '0');
%(idle_threshold = b*mean(ISI_strema)); b_idle_threshold corresponds to b

number_of_surrogates = 25; % mumber of surrogate sequences to be generated
% in order to assess statistical significance of the CFI_MI index

T_total = 2000; % duration of simulation in [ms]
alpha = 0.7; % scaling coefficient of synapsis matrix S
cell_ind1 = 10; % for example, 10th cell
cell_ind2 = 100; % for example, 100th cell

%% Realistic spiking model for producing coupled cortical dynamics, see REF1: Izhikevich EM (2003) Simple model of spiking neurons. IEEE Transactions onneural networks 14(6):1569–1572
% Simulation of a network of 1000 randomly coupled spiking neurons:
% 800 regular spiking (RS) - excitatory cells and 200 low-threshold spiking (LTS) - inhibitory cells

%--------------------------------------------------------------------------
% Part of code created by Eugene M. Izhikevich, February 25, 2003
Ne=800; % excitatory neurons
Ni=200; % inhibitory neurons
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri];
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];
S=[alpha*rand(Ne+Ni,Ne), -alpha*rand(Ne+Ni,Ni)];
v=-65*ones(Ne+Ni,1); % initial values of v
u=b.*v; % initial values of u
firings=[]; % spike timings


for t=1:T_total 
    I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
    fired=find(v>=30); % indices of spikes
    firings=[firings; t+0*fired,fired];
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    I=I+sum(S(:,fired),2);
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u); % stability
end;
% plot(firings(:,1),firings(:,2),'.');
%--------------------------------------------------------------------------

%% Extract cells, complete network of Ne+Ni cells
cells = {};
for j = 1:Ne+Ni
    y = firings(:,1); % timings
    cells{j} = y(find(firings(:,2) == j)); % timings in [ms] per jth cell
end

%% Estimate CFI_MI index between two exemplary cells
spike_train1 = cells{:,cell_ind1}; % spikes in [ms]
spike_train2 = cells{:,cell_ind2}; % spikes in [ms]

CFI_MI = function_CFI_MI(spike_train1, spike_train2, b_idle_threshold); % estimation of CFI_MI index

%% Statistical signficance, see REF2: The JOint DIstribution method (JODI) for generation of surrogate event sequences was originally proposed by L. Ricci et al. in Chaos 29 (2019),121102, <a href="matlab:web('https://doi.org/10.1063/1.5138250')">doi:10.1063/1.5138250</a>

% Generates surrogate sequences of Inter-Event-Intervals (IEI) corresponding to the original
% sequence stored in the ieiSequence (isiSurrogates_1 or isiSurrogates_2) array by means of the
% JOint DIstribution (JODI) algorithm.

ISI_1 = diff(spike_train1);
ISI_2 = diff(spike_train2);

isiSurrogates_1 = spiSeMe_surrogate_jodi(ISI_1, 'M', number_of_surrogates);
isiSurrogates_2 = spiSeMe_surrogate_jodi(ISI_2, 'M', number_of_surrogates);

CFI_MI_surrogates_ALL = [];
for k = 1 : size(isiSurrogates_1,2)
    %----------------------------------------------------------
    ISI_surrogates_1 =  isiSurrogates_1(:,k);
    spike_train_surrogates_1 = [0];
    for iii = 2 : numel(ISI_surrogates_1)
        spike_train_surrogates_1(iii) = ISI_surrogates_1(iii) + spike_train_surrogates_1(iii-1);
    end
    spike_train_surrogates_1 = spike_train_surrogates_1';
    %----------------------------------------------------------
    ISI_surrogates_2 =  isiSurrogates_2(:,k);
    spike_train_surrogates_2 = [0];
    for iii = 2 : numel(ISI_surrogates_2)
        spike_train_surrogates_2(iii) = ISI_surrogates_2(iii) + spike_train_surrogates_2(iii-1);
    end
    spike_train_surrogates_2 = spike_train_surrogates_2';
    %----------------------------------------------------------
    %  MI index on surrogates!
    CFI_MI_surrogates = function_CFI_MI(spike_train_surrogates_1, spike_train_surrogates_2, b_idle_threshold);
    CFI_MI_surrogates_ALL = [CFI_MI_surrogates_ALL; CFI_MI_surrogates];
end % k number of surrogates!

%% significant or not?
percentile_CFI_MI_all_surrogates_one_original_pair  = prctile(CFI_MI_surrogates_ALL,[2.5, 97.5]); % two-tailed hypothesis test with 5% significance

down_limit = percentile_CFI_MI_all_surrogates_one_original_pair(:, 1);
upper_limit = percentile_CFI_MI_all_surrogates_one_original_pair(:, 2);

if (CFI_MI >= down_limit && CFI_MI <= upper_limit)
    flags_CFI_MI =  0; % flag = 0 --> no correlation
elseif (CFI_MI <  down_limit)
    flags_CFI_MI =  -1; % flag = -1 --> anti-correlation
elseif (CFI_MI > upper_limit)
    flags_CFI_MI =  1; % flag = 1 --> correlation
else
end

if flags_CFI_MI == 1
    title1 = strcat('CFI_{MI} =', num2str(CFI_MI), '; statistically significant');
elseif flags_CFI_MI == - 1
    title1 = strcat('CFI_{MI} =', num2str(CFI_MI), '; is statistically significant')'
elseif flags_CFI_MI == 0
    title1 = strcat('CFI_{MI} =', num2str(CFI_MI), '; no statistically significant');
else
end

%% Visualize spikes trains and binary profiles with indication of CFI_MI value and its statistical significance

[spike_train1, states1] = binary_representation(spike_train1, b_idle_threshold);
[spike_train2, states2] = binary_representation(spike_train2, b_idle_threshold);

figure(1);
subplot(2,1,1);
plot(spike_train1, 0.5*ones(size(spike_train1,1),1),'k.'); hold all; % train 1
stairs(spike_train1, states1, 'Color','r','LineWidth',1.25,'LineStyle',':'); hold all;
xlim ([0 T_total]); ylim([min(states1)-0.5 max(states1)+0.5]);
ylabel('Spike train 1','FontSize',12,'FontWeight','bold');
subplot(2,1,2); xlabel('time [s]'); ylabel('spike__train1_binary');
plot(spike_train2, 0.5*ones(size(spike_train2,1),1),'k.'); hold all; % train 1
stairs(spike_train2, states2, 'Color','r','LineWidth',1.25,'LineStyle',':'); hold all;
xlim ([0 T_total]); ylim([min(states2)-0.5 max(states2)+0.5]);
ylabel('Spike train 2','FontSize',12,'FontWeight','bold'); 
xlabel('Time [ms]','FontSize',12,'FontWeight','bold');
title(title1, 'FontSize',12,'FontWeight','bold');





