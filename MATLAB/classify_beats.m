function [n_abnorm, n_beats_det] = classify_beats(ecg, SVMModel)

%% Sampling frequency
fs = 360;
ts = 1/fs;

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Extract QRS complex %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%% Pan-Tompkin algorithm %%%%%%%%%%%%%%%%%%%%%%%%%

%% DC drift cancellation 
ecg_m = ecg - mean(ecg); %cancel DC conponents

%% Normalization
ecg_m = ecg_m./max(abs(ecg_m));

%% Low-pass filtering by filter (1-z^-6)^2/(1-z^-1)^2
b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a = [1 -2 1];
h_LP = filter(b,a,[1 zeros(1,12)]); % Transfer function

ecg_LP = conv(ecg_m, h_LP, 'same');
ecg_LP = ecg_LP/max(abs(ecg_LP)); % normalize

%% High-pass filtering by filter (1-z^-6)^2/(1-z^-1)^2 z^-16-[(1-z^-32)/(1-z^-1)]
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];
h_HP = filter(b,a,[1 zeros(1,32)]); % Transfer function

ecg_HP = conv(ecg_LP, h_HP, 'same');
ecg_HP = ecg_HP/max(abs(ecg_HP)); % normalize

%% Derivative Filter
% Impulse response
h_d = [-1 -2 0 2 1]/8;

% Apply filter
ecg_D = conv(ecg_HP, h_d, 'same');
ecg_D = ecg_D/max(abs(ecg_D)); % normalize

%% Squaring
ecg_SQ = ecg_D.^2;
ecg_SQ = ecg_SQ./max(abs(ecg_SQ)); % normalize

%% Moving Window Integration
% Impulse response
h_MW = ones(1,31)/31;

% Apply filter
ecg_MW = conv(ecg_SQ, h_MW, 'same');
ecg_MW = ecg_MW/max(abs(ecg_MW));

%% Apply threshold
poss_reg =(ecg_MW > mean(ecg_MW))';

%% Find boundary of QRS complex
left_ind = find(diff([0 poss_reg.'])==1).';
right_ind = find(diff([poss_reg.' 0])==-1).';

%% Find QRS peaks
R_value = zeros(numel(left_ind), 1);
R_loc = zeros(numel(left_ind), 1);

Q_value = zeros(numel(left_ind), 1);
Q_loc = zeros(numel(left_ind), 1);

S_value = zeros(numel(left_ind), 1);
S_loc = zeros(numel(left_ind), 1);

for i = 1:numel(left_ind)
    [R_value(i), R_loc(i)] = max(ecg_m(left_ind(i):right_ind(i)));
    R_loc(i) = R_loc(i)-1+left_ind(i); % add offset

    [Q_value(i), Q_loc(i)] = min(ecg_m(left_ind(i):R_loc(i)));
    Q_loc(i) = Q_loc(i)-1+left_ind(i); % add offset

    [S_value(i), S_loc(i)] = min(ecg_m(R_loc(i):right_ind(i)));
    S_loc(i) = S_loc(i)-1+R_loc(i); % add offset
end

%% %%%%%%%%%%%%%%%%%%%%% Discard incomplete ECG beats %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if R_loc(1)*ts < .2
    Q_loc(1) = [];
    R_loc(1) = [];
    S_loc(1) = [];
    
    left_ind(1) = [];
    right_ind(1) = [];
end
if (numel(ecg)-R_loc(end))*ts < .2
    Q_loc(end) = [];
    R_loc(end) = [];
    S_loc(end) = [];
    
    left_ind(end) = [];
    right_ind(end) = [];
end
n_beats_det = numel(R_loc);

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Extract P peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select boundary for P peaks
t_PQ = 0.2;
P_right_ind = left_ind;
P_left_ind = left_ind - round(t_PQ*fs);
if P_left_ind(1) < 1
    P_left_ind(1) = 1;
end

%% Find P peaks
P_value = zeros(numel(R_loc), 1);
P_loc = zeros(numel(R_loc), 1);

for i = 1:numel(P_loc)
    [P_value(i), P_loc(i)] = max(ecg_m(P_left_ind(i):P_right_ind(i)));
    P_loc(i) = P_loc(i)-1+P_left_ind(i); % add offset
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Extract T peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select boundary for T peaks
t_ST = 0.4;
T_left_ind = right_ind;
T_right_ind = T_left_ind + t_ST*fs;
if T_right_ind(end) > numel(ecg_m)
    T_right_ind(end) = numel(ecg_m);
end

%% Find T peaks
T_value = zeros(numel(R_loc), 1);
T_loc = zeros(numel(R_loc), 1);

for i = 1:numel(T_loc)
    [T_value(i), T_loc(i)] = max(ecg_m(T_left_ind(i):T_right_ind(i)));
    T_loc(i) = T_loc(i)-1+T_left_ind(i); % add offset
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature extraction %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RS Width %%%%%%%%%%%%%%%%%%%%%%%%%%%%
RS_width = ts*(S_loc-R_loc);

%% %%%%%%%%%%%%%%%%%%%%%%% Mean Power Spectral density %%%%%%%%%%%%%%%%%%%%
MPSD = zeros(n_beats_det, 1);

for i = 1:n_beats_det
    ecg_trunc = ecg_m(P_loc(i):T_loc(i));
    MPSD(i) = mean(abs(fft(ecg_trunc)).^2);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Area Under QR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
area_QR = zeros(n_beats_det, 1);

for i = 1:n_beats_det
    ecg_trunc_QR = ecg_m(Q_loc(i):R_loc(i));
    area_QR(i) = trapz(ecg_trunc_QR);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Area Under RS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
area_RS = zeros(n_beats_det, 1);

for i = 1:n_beats_det
    ecg_trunc_RS = ecg_m(R_loc(i):S_loc(i));
    area_RS(i) = trapz(ecg_trunc_RS);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QS Width %%%%%%%%%%%%%%%%%%%%%%%%%%%%
QS_width = ts*(S_loc-Q_loc);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-RR Interval %%%%%%%%%%%%%%%%%%%%%%%%%
pre_RR_int = [0; ts*(R_loc(2:end))-ts*(R_loc(1:end-1))];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-RR Interval %%%%%%%%%%%%%%%%%%%%%%%%
post_RR_int = [pre_RR_int(2:end); 0];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QS Width %%%%%%%%%%%%%%%%%%%%%%%%%%%%
QR_width = ts*(R_loc-Q_loc);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature vector %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Current record feature vector
feature_vector = [  QS_width,...
                    pre_RR_int,...
                    post_RR_int,...
                    QR_width,...
                    RS_width,...
                    MPSD,...
                    area_QR,...
                    area_RS,...
                    ];

%% Test SVM model
[y_pred, ~] = predict(SVMModel, feature_vector);
n_abnorm = sum(y_pred==0);






