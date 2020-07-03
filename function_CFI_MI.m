function CFI_MI = function_CFI_MI(spike_train1, spike_train2, b_idle_threshold)

% function_CFI_MI_ver2.m estimated concurrent firing activity based on
% modified mutual information applied on binary profiles (states of quiescence '0' and state of firing '1') of neural spiking activity

% Input parameters:
% 1. spike_train1
% 2. spike_train2
% 3. b_idle_threshold - parameter needed to estimate periods of neural quiescence (state '0'); 
% (idle_threshold = b*mean(ISI_strema)); % b_idle_threshold corresponds to b

% Output parameters:
% 1. CFI_MI index

% This function is part of the SFI MI toolbox.

%% Input parsing & validation
if (~iscolumn(spike_train1))
   spike_train1 = spike_train1';
end
if (~iscolumn(spike_train2))
   spike_train2 = spike_train2';
end

if (b_idle_threshold < 1)
    error('Function argument "b_idle_threshold" must be a positive integer.')
end
%% Two-state representation of both spike trains, binary profiles
[spike_train1, states1] = binary_representation(spike_train1, b_idle_threshold);
[spike_train2, states2] = binary_representation(spike_train2, b_idle_threshold);

%% Taking overlap time interal of two spike trains in order to mutually compare
t1 = max(min(spike_train1), min(spike_train2)); % begin
t2 = min(max(spike_train1), max(spike_train2)); % end

total_recording_time = t2 - t1; % overlapping interval

spike_train11 = spike_train1(find(spike_train1 >= t1 & spike_train1 <= t2));
spike_train22 = spike_train2(find(spike_train2 >= t1 & spike_train2 <= t2));
states11 = states1(find(spike_train1 >= t1 & spike_train1 <= t2));
states22 = states2(find(spike_train2 >= t1 & spike_train2 <= t2));

spike_train11_ALL = []; spike_train11_ALL = [spike_train11_ALL t1; spike_train11; t2];
states11_ALL = []; states11_ALL = [states11_ALL; max(states1(find(spike_train1 <= t1))); states11; min(states1(find(spike_train1 >= t2)))];
spike_train22_ALL = []; spike_train22_ALL = [spike_train22_ALL t1; spike_train22; t2];
states22_ALL = []; states22_ALL = [states22_ALL; max(states2(find(spike_train2 <= t1))); states22; min(states2(find(spike_train2 >= t2)))];

spike_train111 = spike_train11_ALL(find(fliplr([true, diff(fliplr(spike_train11_ALL'))~=0]))); % to exclude a possible repeating of first or last spike (identical samples beacuse of t1 and t2)
states111 = states11_ALL(find(fliplr([true, diff(fliplr(spike_train11_ALL'))~=0])));
spike_train222 = spike_train22_ALL(find(fliplr([true, diff(fliplr(spike_train22_ALL'))~=0]))); % to exclude a possible repeating of first or last spike (identical samples beacuse of t1 and t2)
states222 = states22_ALL(find(fliplr([true, diff(fliplr(spike_train22_ALL'))~=0])));

%% Final streams of spikes and binary states
spike_train1 = spike_train111; 
spike_train2 = spike_train222;
states1 = states111; 
states2 = states222;

possible_states = [0:1: max(max(states1), max(states2))]; % 0 and 1
spike_trains_labels = {'spike_train1', 'spike_train2'};

%%                          MARGINAL PROBABILITIES

marginal_prob_both_trains = [];
for i = 1 : numel(spike_trains_labels) % two spike trains
    
    if spike_trains_labels{:,i} == 'spike_train1' % train1
        spike_train = spike_train1;
        states = states1;
    else % train 2
        spike_train = spike_train2;
        states = states2;
    end
    
    marginal_prob_one_train = [];
    for ii = 1 : numel(possible_states)
        state_temp = possible_states(:, ii);
        xx = find(states' == state_temp); % temporary state of train 1
        c1 = diff([0 find([diff(find(states' == state_temp)) inf]>1)]); %length of the sequences
        d1 = cumsum(c1); %endpoints of the sequences
        spike_ranges = [];
        if ~isempty(xx)
            for f = 1 : numel(d1)
                upper_limit = xx(d1(f)) + 1;
                down_limit = xx(d1(f)) + 1 - c1(f);
                if upper_limit > numel(spike_train)
                    upper_limit = xx(d1(f));
                    spike_ranges = [spike_ranges; spike_train(down_limit) spike_train(upper_limit)]; % all time-intervals under temporary state
                else
                    spike_ranges = [spike_ranges; spike_train(down_limit) spike_train(upper_limit)]; % all time-intervals under temporary state
                end
            end
        end
        
        if ~isempty(spike_ranges)
            marginal_prob = sum(spike_ranges(:,2) - spike_ranges(:,1)) / total_recording_time;
        else
            marginal_prob = 0;
        end
        
        marginal_prob_one_train = [marginal_prob_one_train; marginal_prob];
        % save here 'spike_ranges' because of next step - estimating joint probabilities
        if spike_trains_labels{:,i} == 'spike_train1'
            S_struct.(['spike_ranges_train1_' num2str(ii-1)]) = spike_ranges;
        end
        
    end 
    
    marginal_prob_both_trains = [marginal_prob_both_trains marginal_prob_one_train];
    
end 

%%                          JOINT PROBABILITIES

joint_prob_ALL = [];

for j = 1 : numel(possible_states)  % possible states of spike train 1
    
    fns = fieldnames(S_struct);
    spike_ranges = S_struct.(fns{j});
    
    joint_prob_temp = [];
    for jj = 1 : numel(possible_states) % possible states of spike train 2
        state_temp_another_train = possible_states(:,jj); % temporary state of spike train 2
        
        interval_under_ramps = [];
        for jjj = 1 : size(spike_ranges,1)
            
            temp = spike_ranges(jjj,:);
            
            yy = states2(find(spike_train2 >= temp(1,1) & spike_train2 <= temp(1,2))); % all states of another train (spike train 2)
            zz = spike_train2(find(spike_train2 >= temp(1,1) & spike_train2 <= temp(1,2))); % all spikes of another train (spike train 2)
            
            %----------------------------------------------------------------------
            zz_ALL = [];  yy_ALL = [];
            zz_ALL = [zz_ALL; temp(1,1); zz; temp(1,2)];
            zz_ALL = sort(zz_ALL);
            
            for s = 1: numel(zz_ALL)
                xx = find(spike_train2 <= zz_ALL(s,:));
                yy_ALL = [yy_ALL; states2(xx(end, :))];
            end
            
            zz = zz_ALL;
            yy = yy_ALL;
            
            %======================================================================
            % another train is in state 0 or 1
            %======================================================================
            
            xx = find(yy' == state_temp_another_train);
            if ~isempty(xx)
                c1 = diff([0 find([diff(find(yy' == state_temp_another_train)) inf]>1)]); %length of the sequences
                d1 = cumsum(c1); %endpoints of the sequences
                for f = 1:numel(d1)
                    upper_limit = xx(d1(f)) + 1;
                    down_limit = xx(d1(f)) + 1 - c1(f);
                    if upper_limit > numel(zz)
                        upper_limit = xx(d1(f));
                        interval_under_ramps = [interval_under_ramps; zz(down_limit) zz(upper_limit)];
                    else
                        interval_under_ramps = [interval_under_ramps; zz(down_limit) zz(upper_limit)];
                    end
                end
            end 
        end 
        
        if ~isempty(interval_under_ramps)
            joint_prob = sum(interval_under_ramps(:,2) - interval_under_ramps(:,1)) / total_recording_time;
        else
            joint_prob = 0;
        end
        
        joint_prob_temp = [joint_prob_temp; joint_prob]; % first state of spike train 1 with all states of spike train 2
    end 
    
    joint_prob_ALL = [joint_prob_ALL  joint_prob_temp]; % all states of spike train 1 with all states of spike train 2
    
end 

% %% just to check
% just_to_check_joint = sum(sum(joint_prob_ALL)); % must be equal to 1
% just_to_check_marginals = sum(marginal_prob_both_trains); % must be equal to 1, 1

%%                          MUTUAL INFORMATION

marginal_prob_train1 = marginal_prob_both_trains(:,1);
marginal_prob_train2 = marginal_prob_both_trains(:,2);

MI_array = [];
for k = 1 : size(joint_prob_ALL, 2)
    part = joint_prob_ALL(:, k);
    for kk = 1 : size(part,1)
        joint_temp = part(kk, :);
        if (joint_temp ~=0 && marginal_prob_train1(k,:)~=0 && marginal_prob_train2(kk,:)~=0)
            MI_term = joint_temp*(log2((joint_temp /(marginal_prob_train1(k,:)*marginal_prob_train2(kk,:)))));
        else
            MI_term = 0;
        end
        MI_array = [MI_array; MI_term];
    end
end

MI = sum(MI_array);

%%                        ENTROPIES OF BOTH TRAINS

H1_array = []; % spike train 1
for kkk = 1 : size(marginal_prob_train1,1)
    part = marginal_prob_train1(kkk,:);
    if part~=0
        term = - part*log2(part);
        H1_array = [H1_array; term];
    end
end

entropy_1 = sum(H1_array);

H2_array = []; % spike train 2
for kkk = 1 : size(marginal_prob_train2,1)
    part = marginal_prob_train2(kkk,:);
    if part~=0
        term = - part*log2(part);
        H2_array = [H2_array; term];
    end
end

entropy_2 = sum(H2_array);

%% Average conditional probabilities

% marginal_idle1 = marginal_prob_train1(1,:);
marginal_work1 = marginal_prob_train1(2,:);
marginal_idle2 = marginal_prob_train2(1,:);
marginal_work2 = marginal_prob_train2(2,:);

joint_idle_idle = joint_prob_ALL(1,1);
joint_work_work = joint_prob_ALL(2,2);
joint_idle_work = joint_prob_ALL(2,1);
joint_work_idle = joint_prob_ALL(1,2);

p_c = 0; p_ac = 0;
if (marginal_idle2 ~=0 && marginal_work2~=0)
    p_c = 0.5*((joint_idle_idle / marginal_idle2) + (joint_work_work / marginal_work2));
    p_ac = 0.5*((joint_idle_work / marginal_work2) + (joint_work_idle / marginal_idle2));
end

%% Computation of Concurrent Firing activity based on normalized mutual information, CFI_MI index
if marginal_work1 == marginal_work2 % both trains are steadily in state 1
    CFI_MI = 1;
elseif joint_work_idle == 1 || joint_idle_work == 1 % both trains are steadily in opposite states
    CFI_MI = -1;
elseif MI == 0  % if MI = 0 --> CFI_MI = 0
    CFI_MI = MI;
elseif p_c > p_ac
    CFI_MI = MI / (min(entropy_1, entropy_2));
elseif p_c < p_ac
    CFI_MI = - MI / (min(entropy_1, entropy_2));
elseif p_c == p_ac
    CFI_MI = 0;
else
end

%% Ploting trains and binary profiles
figure;
subplot(2,1,1);
plot(spike_train1, 0.5*ones(size(spike_train1,1),1),'k.'); hold all; % train 1
stairs(spike_train1, states1, 'Color','r','LineWidth',1.25,'LineStyle',':'); hold all;
xlim([t1 t2]); ylim([min(states1)-1 max(states1)+1]);
ylabel('Spike train 1','FontSize',12,'FontWeight','bold');
subplot(2,1,2); xlabel('time [s]'); ylabel('spike__train1_binary');
plot(spike_train2, 0.5*ones(size(spike_train2,1),1),'k.'); hold all; % train 1
stairs(spike_train2, states2, 'Color','r','LineWidth',1.25,'LineStyle',':'); hold all;
xlim([t1 t2]); ylim([min(states2)-1 max(states2)+1]);
ylabel('Spike train 2','FontSize',12,'FontWeight','bold');
xlabel('Time [s]','FontSize',12,'FontWeight','bold');
title(strcat('CFI_{MI}=', num2str(CFI_MI)), 'FontSize',12,'FontWeight','bold');

