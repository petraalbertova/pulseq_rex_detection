% ----------------------------------------------------------------------------------------------------
%  Rotary Excitation of non-sinusoidal pulsed magnetic fields:
%  Toward noninvasive direct detection of cardiac conduction
%  Petra Albertova, Maximilian Gram, Martin Blaimer, Wolfgang R. Bauer, Peter M. Jakob, Peter Nordbeck
%  Magnetic Resonance in Medicine; doi: 10.1002/mrm.30190
% ----------------------------------------------------------------------------------------------------

%  Contact: 
%  Petra Albertova, University of Würzburg,  
%  Experimental Physics 5,  
%  Am Hubland, Würzburg,  
%  DE D-97074,  
%  petra.albertova@uni-wuerzburg.de 

% V1; June 10 2024

%% initialize basic pulseq objects
clear;

% our system specifications: Siemens MAGNETOM Skyra
% carefully check your MRI system specifications and adjust!
system = mr.opts( 'MaxGrad',         45, ...
                  'GradUnit',        'mT/m', ...
                  'MaxSlew',         200, ...
                  'SlewUnit',        'T/m/s', ...
                  'rfRingdownTime',  100e-6, ...
                  'rfDeadTime',      100e-6, ...
                  'adcDeadTime',     10e-6, ...
                  'gamma',           2.6752218744 *1e8 /2/pi, ...
                  'B0',              3.0 );

% create seq object
seq = mr.Sequence(system);

% set definitions
scan_id = datestr(datetime('now'), 'yymmddHHMM');
seq.setDefinition('Name',       'rex_detection');
seq.setDefinition('Scan_ID',    int64(str2num(scan_id(2:end)))); % 9 digits can be saved in the wip-mem-block
seq.setDefinition('FOV',        [0.24, 0.24, 0.05]);             % we used a FOV with 240x240mm and 5mm slice thickness
seq.setDefinition('Rot_Matrix', [1 0 0; 0 1 0; 0 0 1]);          % rotation can be done at the scanner

%% slice excitation

% e.g. with a 90° sinc pulse,
% 2ms exciatiot time,
% time bandwith product 6
% slice thickness 5mm

% note: we used SLR optimized excitation pulses.
% install the SIGPY package via 'pip install sigpy'

[rf_exc, gz, gz_reph] = mr.makeSincPulse( pi/2, ...
                                          system, ...
                                          'Duration', 2*1e-3, ...
                                          'timeBwProduct', 6, ...
                                          'SliceThickness', 5*1e-3, ...
                                          'use', 'excitation');

%% spiral readout

% note: this is an interleaved readout using 4 center-out spirals
% our FOV was 240x240mm with a 128x128 matrix
% each adc uses 5504 samples at 400kHz sampling rate
% if you want or need adaption of those parameters, install the variable
% densitiy toolbox provided by Brian Hargreaves;
% http://mrsrl.stanford.edu/~brian/vdspiral/vdspiral.tar.gz
load('spiral_readout.mat')
adc       = mr.makeAdc(5504, 'Dwell', 2.5*1e-6); 
adc.delay = system.adcDeadTime;

% gradient spoilig
% e.g. 2ms, 8x2Pi Twists inside the 5mm slice
% note: we performed gradient spoiling in z direction
% the spiral gradients were "looped-in" to balance the x and y moments
gz_spoil = mr.makeTrapezoid('z', 'Area', 8 / 0.005, 'maxGrad', system.maxGrad*0.75, 'maxSlew', system.maxSlew*0.75, 'system', system);

%% spin-lock preparation
tSL = 100 *1e-3; % [s]  spin-lock pulse duration
fSL = 25;        % [Hz] spin-lock pulse amplitude

% create block pulse waveform for spin-locking
% note: we observed problems for spin-lock pulses which were implemented
% with the mr.makeBlockPulse function of the Pulseq Git. We assume that
% pulseq block pulses are not suitable for longer duration continuous
% waves, since the waveforms only consits on two steps
Nringdown = round(system.rfRingdownTime/system.rfRasterTime);
N         = round(tSL/system.rfRasterTime) + Nringdown;
t         = (1:N)' * system.rfRasterTime;
signal    = ones(N,1) * fSL;
if Nringdown>0
    signal(end-Nringdown+1:end,1) = 0.0;    
end
SL.type         = 'rf';
SL.t            = t;
SL.signal       = signal;
SL.freqOffset   = 0;
SL.phaseOffset  = 0;
SL.deadTime     = system.rfDeadTime;
SL.ringdownTime = system.rfRingdownTime;
SL.delay        = system.rfDeadTime;
SL.shape_dur    = SL.t(end);
SL.use          = 'preparation';
clear Nringdown N t signal;

% load adiabatic excitation pulse
% we used an B0 and B1+ optimized AHP pulse for excitation prior to
% spin-locking; duration = 3ms; max f1 = 600Hz; 
load('ahp_pulse.mat');

% crusher gradiets after spin-locking
% e.g. 8x2Pi Twists in each voxel
gx_crush = mr.makeTrapezoid('x', 'Area', 8/(0.24/128), 'maxGrad', system.maxGrad * 0.75, 'maxSlew', system.maxSlew * 0.75, 'system', system);
gy_crush = mr.makeTrapezoid('y', 'Area', 8/(0.24/128), 'maxGrad', system.maxGrad * 0.75, 'maxSlew', system.maxSlew * 0.75, 'system', system);
gz_crush = mr.makeTrapezoid('z', 'Area', 8/(0.005),    'maxGrad', system.maxGrad * 0.75, 'maxSlew', system.maxSlew * 0.75, 'system', system);
gy_crush.delay = ceil(gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;                   % reduce effective slew rate
gz_crush.delay = gy_crush.delay + ceil(gy_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;  % reduce effective slew rate

%% transmission of Rotary EXcitation (tREX)
% in this section we will calculate gradient objects (e.g. in z direction)
% which are applied during spin-locking; thoses gradients are used to
% emulate the presence of a biomagnetic field interacting with the
% spin-locked magnetization; the gradients are set to extremely small values
% and by shifting the measurement slice close to the iso-center, it is
% possible to probe the rotary excitation effect for fields in the lower nT
% range. this concept was previously publisehd by our group in:
% Towards robust in vivo quantification of oscillating biomagnetic fields using Rotary Excitation based MRI
% Sci Rep. 2022 Sep 13;12(1):15375. doi: 10.1038/s41598-022-19275-5

% choose stimulus type
stim_type = 'sinc'; % 'sine', 'gauss', 'sinc', 'qrs'

% choose stimulus peak magnitude
stim_mag = 100 *1e-9; % [T]

% choose slice offset
dz = 10 *1e-3; % [m]

% choose number of interactions
% for a detection experiment, we need to repeat the interaction with
% different relative time shifts between spin-lock and stimulus
Ninter = 10;

% calculate time shifts
t_shift = linspace(-1/2, 1/2, Ninter) / fSL; % this is a good choice for most experiments

% calculate tREX gradients
figure(); hold on
title(['Magnetic Stimuli: ' stim_type]);
for j=1:Ninter

    % calculate waveform
    if strcmp(stim_type, 'sine')
        fstim = 25; % choose oscillation frequency
        y     = sin(2*pi*100*((1:round(tSL/system.gradRasterTime))*system.gradRasterTime-t_shift(j)));
    else
        if strcmp(stim_type, 'gauss')
            tau_stim = 50 *1e-3; % choose stimulus duration
            BW_stim  = 6;        % choose time bandwitdh product
            y        = get_gauss_stimulus(tau_stim, BW_stim, system);
        elseif strcmp(stim_type, 'sinc')
            tau_stim = 50 *1e-3; % choose stimulus duration
            BW_stim  = 6;        % choose time bandwitdh product
            y        = get_sinc_stimulus(tau_stim, BW_stim, system);
        elseif strcmp(stim_type, 'qrs')
            QRS_stim = 100 *1e-3; % choose QRS duration
            y        = get_qrs_stimulus(QRS_stim, system);    
        end
        n_tSL   = round(tSL/system.gradRasterTime);
        n_tau   = numel(y);
        n_zero  = n_tSL - n_tau;
        y       = [zeros(floor(n_zero/2),1); y; zeros(ceil(n_zero/2),1)]; % zero padding
        n_shift = round(t_shift(j)/system.gradRasterTime);
        if n_shift>0
            y(end-n_shift+1:end) = [];  % right: zero deleting
            y = [zeros(n_shift,1); y];  % left:  zero padding
        end
        if n_shift<0
            y(1:abs(n_shift)) = [];          % left:  zero deleting
            y = [y; zeros(abs(n_shift),1)];  % right: zero padding
        end
        if numel(y) > n_tSL
            n_diff = numel(y)-n_tSL;
            y      = y(ceil(n_diff/2): end - ceil(n_diff/2));
        end        
        clear n_tSL n_tau n_zero n_shift n_diff;
        y = y.'; % use row vector
    end

    % adjust magnitude
    y = y / max(abs(y(:)));
    y = y * stim_mag / dz;
    plot((1:numel(y))*system.gradRasterTime*1e3, y*1e9/1e3, 'LineWidth', 2); % [ms], [nT/mm]
    xlabel('spin-lock duration tSL [ms]', 'FontSize', 12);
    ylabel('tREX gradient [nT/mm]', 'FontSize', 12);
    set(gca, 'FontWeight', 'bold');

    % tREX gradient objects
    g_trex(j)       = mr.makeArbitraryGrad('z', y);
    g_trex(j).first = 0;
    g_trex(j).last  = 0;
    clear y;

end

% perform a simple fourier analysis of the magnetic stimulus
get_fourier_analysis(g_trex(1).tt, g_trex(1).waveform, fSL);

%% shift slice offset
% note: don't change the slice position at your scanner
% this would change the effective stimulus field strength
rf_exc.freqOffset = gz.amplitude * dz;

%% add recovery delay between readout and spin-lock
Trec = mr.makeDelay(1);

%% create pulseq sequence

% start with some magnetization resets for constant longitudinal magnetization
for j=1:5
    seq.addBlock(AHP);
    seq.addBlock(gx_crush, gy_crush, gz_crush);
    seq.addBlock(Trec);
end

% start of tREX sequence
for j=1:Ninter
    for k=1:numel(gx)
        seq.addBlock(AHP);
        seq.addBlock(SL, g_trex(j));
        seq.addBlock(gx_crush, gy_crush, gz_crush);
        seq.addBlock(rf_exc, gz);
        seq.addBlock(gz_reph);
        seq.addBlock(gx(k), gy(k), adc);
        seq.addBlock(gz_spoil);
        seq.addBlock(Trec);
    end
    % calculate kspace trajectory
    if j==1
        warning('OFF', 'mr:restoreShape');
        [ktraj_adc, ~, ktraj_full] = seq.calculateKspacePP();
    end
end

%% perform timing and gradient checks
[timings_ok, error_report] = seq.checkTiming;
test_report                = seq.testReport;

disp(' '); disp(' ');
disp('---------- check timings and gradient limits!!! ----------');
disp(' ');
if (timings_ok)
    disp('Timing check passed successfully');
else
    disp(   'Timing check failed! Error listing follows:');
    fprintf([error_report{:}]);
    fprintf('\n');
end
disp(' '); disp(' ');
for j=1:numel(test_report)
    fprintf(test_report{1,j});
end
disp(' ');
disp('----------------------------------------------------------');
disp(' ');

%% write pulseq file
if timings_ok
    seq.write(['external_' scan_id '_' stim_type '_' num2str(round(tSL*1e3)) 'ms_' num2str(round(fSL)) 'Hz.seq']);
end

%% plot sequence diagram
seq.plot();
TotalDuration = sum(seq.blockDurations);

%% plot kspace trajectory
figure();
hold on
plot( ktraj_full(1,:), ktraj_full(2,:), 'b-');
plot( ktraj_adc(1,:),  ktraj_adc(2,:), 'r.');
axis('equal');
xlabel('kx [1/m]');
ylabel('ky [1/m]');
xlim([-1 1]*max(abs(ktraj_adc(:)*1.05)))
ylim([-1 1]*max(abs(ktraj_adc(:)*1.05)))
title('kspace trajectory');
set(gca, 'FontWeight', 'bold');

%% additional functions for waveform generation

function y = get_gauss_stimulus(tau_stim, BW_stim, system)
    system.rfRasterTime = system.gradRasterTime;
    y = mr.makeGaussPulse(pi, system, 'Duration', tau_stim, 'timeBwProduct', BW_stim, 'use', 'excitation');
    y = real(y.signal)';
end

function y = get_sinc_stimulus(tau_stim, BW_stim, system)
    system.rfRasterTime = system.gradRasterTime;
    y = mr.makeSincPulse(pi, system, 'Duration', tau_stim, 'timeBwProduct', BW_stim, 'apodization', 0.5, 'use', 'excitation');
    y = real(y.signal)';
end

function yq = get_qrs_stimulus(QRS_stim, system)
    raw    = load('qrs.mat');
    y      = raw.qrs;
    t_full = QRS_stim / 0.4382; % since QS duration is 43.82 % of full duration
    t      = linspace(0, 1, numel(y))';
    tq     = linspace(0, 1, round(t_full/system.gradRasterTime))';
    yq     = squeeze(interp1(t, y, tq));
end

function [P, f] = get_fourier_analysis(t, y, fSL)

dt = mean(diff(t));
Fs = 1 / dt;
n  = length(y);
y  = [zeros(1,10*n), y, zeros(1,10*n)];
n  = length(y);
Y  = fft(y);
P  = abs(Y).^2/n;
f  = (0:n-1)*(Fs/n);

figure;
plot(f, P, 'k-', 'LineWidth', 2);
xline(fSL, 'r--', 'spin-lock frequency fSL', 'LineWidth', 2, 'LabelOrientation', 'horizontal')
xlim([0 200]);
title('Power Spectral Density of Stimulus');
xlabel('Frequency [Hz]');
ylabel('Power/Frequency [dB/Hz]');
grid on;
set(gca, 'FontWeight', 'bold');

end