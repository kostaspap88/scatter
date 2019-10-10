% SCATTER attack based on "SCATTER: a new dimension in side-channel" by 
% Thiebeauld et al.
% author: Kostas Papagiannoulos -- kpcrypto.net


clear all;
close all;

% USER INPUT-------------

no_traces = 100;
oscilloscope_resolution = 8; % in bits
sigma = 0.1;

% END OF USER INPUT------


% DATASET SIMULATION---------------------------------

% PRESENT sbox
sbox = [12 6 5 11 9 0 10 13 3 14 15 8 4 7 1 2];

plaintext = randi(16,no_traces,1)-1;
k_true = 5;
y = bitxor(plaintext, k_true);
z = sbox(y+1)';
% lets add some unrelated samples r so that not all samples are useful to
% the SCATTER attack
r = (randi(16,no_traces,1)-1);
no_samples = 9; % e.g. we observe the leakage of z in 4 samples

% simulate noise with std sigma
noise = normrnd(0, sigma, no_traces, no_samples);

% add noise to compute the analog value of the leakage
% 4 sample points are related to key-dependent value z
% 5 sample points are related to unrelated value r
traces_analog = [ hw(r) hw(z) hw(r) hw(r) hw(z) hw(z) hw(r) hw(z) hw(r) ] + noise;

% assume that the following countermeasure is in place:
% - the device can shuffle every sample in the implementation
% - thus here we generate a <no_samples>-bit random sequence for every
% encryption 
% - then we permute every trace's samples using the random sequence
shuffle_pattern = zeros(no_traces, no_samples);
traces_analog_shuffled = zeros(no_traces, no_samples);
for i=1:no_traces
    shuffle_pattern = randperm(no_samples);
    current_trace = traces_analog(i,:);
    shuffled_trace = current_trace(shuffle_pattern);
    traces_analog_shuffled(i,:) = shuffled_trace;
end


% the analog value will be digitized with quantization levels that depend
% on the oscilloscope resolution
% we assume that our oscilloscope window is set appropriately to capture
% the max and min leakage values, i.e. there is no signal trimming
max_val = max(max(traces_analog));
min_val = min(min(traces_analog));

window_size = abs(max_val - min_val);
step = window_size/(2^oscilloscope_resolution);
levels = linspace(min_val, max_val, 2^oscilloscope_resolution + 1);
bin_edges = levels;

% SCATTER ATTACK---------------------------------------
% we perform an attack on both the normal and the unshuffled traces
% scatter ignores the time dimension, assisting us in the case of shuffled
% or simply misaligned samples



% Accumulator construction
accumulator = zeros(2^oscilloscope_resolution, 16, 5);
accumulator_shuffled = zeros(2^oscilloscope_resolution, 16, 5);
for i=1:no_traces
   d = distribute_trace_data(traces_analog(i,:), bin_edges);
   d_shuffled = distribute_trace_data(traces_analog_shuffled(i,:), bin_edges);
   for g=0:15
       h = hw(sbox(bitxor(plaintext(i), g)+1));
       accumulator(:, g+1, h+1) = accumulator(:, g+1, h+1) + d';
       accumulator_shuffled(:, g+1, h+1) = accumulator_shuffled(:, g+1, h+1) + d_shuffled';
   end
end

% Accumulator normalization
sum_acc = zeros(16, 5);
sum_acc_shuffled = zeros(16, 5);
for h=0:4
    for g = 0:15
        sum_acc(g+1, h+1) = sum(accumulator(:, g+1, h+1));
        sum_acc_shuffled(g+1, h+1) = sum(accumulator_shuffled(:, g+1, h+1));
    end
    accumulator(:,g+1, h+1) = accumulator(:, g+1, h+1)/sum_acc(g+1, h+1);
    accumulator_shuffled(:,g+1, h+1) = accumulator_shuffled(:, g+1, h+1)/sum_acc_shuffled(g+1, h+1);
end

% scatter chi-square distinguisher
chi2 = zeros(16, 5);
chi2_shuffled = zeros(16, 5);
for h=0:4
    for g = 0:15
        for u=0:(2^oscilloscope_resolution - 1)
            term0 = 1/16 * sum(accumulator(u+1, :, h+1));
            term0_shuffled = 1/16 * sum(accumulator_shuffled(u+1, :, h+1));
            
            term1 = (accumulator(u+1, g+1, h+1) - term0)^2;
            term1_shuffled = (accumulator_shuffled(u+1, g+1, h+1) - term0_shuffled)^2;
            
            if (term0==0)&&(term1==0)
                chi2(g+1, h+1) = chi2(g+1, h+1);
                chi2_shuffled(g+1, h+1) = chi2_shuffled(g+1, h+1);
            else
                chi2(g+1, h+1) = chi2(g+1, h+1) + term1/term0;
                chi2_shuffled(g+1, h+1) = chi2_shuffled(g+1, h+1) + term1_shuffled/term0_shuffled;
            end        
        end
    end
end

scatter_chi2 = zeros(16, 1);
scatter_chi2_shuffled = zeros(16, 1);
for g = 0:15
    scatter_chi2(g+1) = prod(chi2(g+1, :));
    scatter_chi2_shuffled(g+1) = prod(chi2_shuffled(g+1, :));
end

[val, index] = max(scatter_chi2);
top_guess_scatter_chi2 = index - 1

[val_shuffled, index_shuffled] = max(scatter_chi2_shuffled);
top_guess_scatter_chi2_shuffled = index_shuffled - 1







