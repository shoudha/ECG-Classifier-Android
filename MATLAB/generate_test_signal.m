function ecg = generate_test_signal(t_select, recs_vec)

%% Sampling frequency
fs = 360;

%% Load data
n_rec = numel(recs_vec);
i_rec = randi([1 n_rec]);

val = load("data/"+recs_vec(i_rec)+"m.mat");
val = val.val;

%% Select window of signal
n = size(val, 2);
if t_select ~= 0
    ns = round(fs*t_select);
else
    ns = n;
end

%% Random window position
i_start = 1;
i_end = numel(val(1,:))-ns+1;

i_select = randi([i_start i_end]);

ecg = val(1,i_select:i_select+ns-1);