clear; clc;
%% sistem parametreleri
pd = 0.9;            
pfa = 1e-6;         
range_res = 50;    
target_rcs = 1;         
prop_speed = physconst('LightSpeed');   
fc = 8e9;
lambda = prop_speed/fc;

%% TxRx parametreleri
Gt = 20;
Gr = 20;
L=1;
Rt= 15000; %güç fazlaysa max menzili değiştir (karesiyle değişir)
Rr= 15000;
Rmax=ceil(sqrt(Rt*Rr));

Urt = [0;0;0];
Urx = [16000;0;0]; % vericiye göre alıcının konumu 
[Lbase,Langle]= rangeangle(Urt,Urx);

%% Pulse tanımları
prf = prop_speed/(2*Rmax);
pulse_width = 2*range_res/prop_speed;
pulse_bw = 1/pulse_width;
fs = 2*pulse_bw; %nyquist
noise_bw=pulse_bw;
num_pulse = 50;

waveform = phased.RectangularWaveform(...
    'PulseWidth',pulse_width,...
    'PRF',prf,...
    'SampleRate',fs);

%% Verici gücü hesabı
snr_min = albersheim(pd, pfa, num_pulse) % N sample'da verilen pd pfa elde edilmesi için gereken min SNR

Pt = (4*pi)^3*Rmax^4*db2pow(snr_min)*...
    L*noisepow(pulse_bw)/(Gt*Gr*lambda^2*target_rcs)


%% Tx kurulumu
transmitter = phased.Transmitter(...
    'Gain',Gt,...
    'PeakPower',Pt);  
txAntenna = phased.ShortDipoleAntennaElement('AxisDirection','Z');
txarray = phased.ULA(4,lambda/2,'Element',txAntenna);

radmotion = phased.Platform(...
   'OrientationAxesOutputPort',true,...
   'OrientationAxes',azelaxes(0,0));

radiator = phased.Radiator(...
    'Sensor',txarray,'PropagationSpeed',prop_speed,...
    'OperatingFrequency',fc, 'Polarization','Combined');
%% Rx kurulumu
receiver = phased.ReceiverPreamp(...
    'Gain',Gr,...
    'NoiseFigure',0,...
    'SampleRate',fs,...
    'EnableInputPort',false);

rxAntenna = phased.ShortDipoleAntennaElement('AxisDirection','Z');
rxarray = phased.ULA(4,lambda/2,'Element',rxAntenna);

recvmotion = phased.Platform(...
    'InitialPosition',Urx,...
    'OrientationAxesOutputPort',true,...
    'OrientationAxes',azelaxes(180,0));%azimuth 180

collector = phased.Collector(...
     'Sensor',rxarray,'PropagationSpeed',prop_speed,...
    'OperatingFrequency',fc, 'Polarization','Combined');

dopplerResponse=phased.RangeDopplerResponse('SampleRate',pulse_bw,...
    'PropagationSpeed',prop_speed,...
    'DopplerOutput','Speed',...
    'OperatingFrequency',fc,...
    'DopplerFFTLengthSource','Property',...
    'DopplerFFTLength',512);

beamformer = phased.PhaseShiftBeamformer(...
    'SensorArray',collector.Sensor,...
    'PropagationSpeed',prop_speed,'OperatingFrequency',fc,...
    'DirectionSource','Input port',...
    'WeightsNormalization','Preserve Power','WeightsOutputPort',true);
%%
%pattern(rxarray,fc,[-180:180],0)
%% Hedef sistem tanımı
TN=3; % hedef sayısı
tgtpos =[[8000;6000;700],[8000;0;100],[9000;15000;8000]];
tgtvel = [[0;100;0],[0;100;0],[100;100;100]];
ScatteringMatrices = {[1 0;0 1];[1 0;0 1];[1 0;0 1]};

for i=1:TN
    target{i} = phased.RadarTarget('EnablePolarization',true,...
        'Mode','Bistatic','ScatteringMatrix',ScatteringMatrices{i},...
        'PropagationSpeed',prop_speed,'OperatingFrequency',fc);
end
tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel,...
        'OrientationAxesOutputPort',true,'OrientationAxes',azelaxes(0,0));

txchannel = phased.FreeSpace('SampleRate',fs,...
    'OperatingFrequency',fc,'PropagationSpeed',prop_speed);
rxchannel = phased.FreeSpace('SampleRate',fs,...
    'OperatingFrequency',fc,'PropagationSpeed',prop_speed);

%% Simülasyon ve noise tanımları
fast_time_grid = unigrid(0,1/fs,1/prf,'[)');
slow_time_grid = (0:num_pulse-1)/prf;
receiver.SeedSource = 'Property';
receiver.Seed = 2007;

rxpulses = zeros(numel(fast_time_grid),num_pulse);

npower = noisepow(noise_bw,receiver.NoiseFigure,receiver.ReferenceTemperature);
threshold = npower * db2pow(npwgnthresh(pfa,num_pulse,'noncoherent'));

% 3d plot
Plot3d(Urt,Urx,tgtpos,Rmax/1000,Langle);

%% simulasyon
for m=1:num_pulse
    % sistemdeki platformların dinamiklerinin hesabı
    [radpos,radvel,radax] = radmotion(1/prf);
    [recvpos,recvel,recvax] = recvmotion(1/prf);
    [tgtpos,tgtvel,tgtax] = tgtmotion(1/prf);
    
    % hedef açısının ve uzaklığının hesabı
    [radrange,radang] = rangeangle(tgtpos,radpos,radax); % gidiş rotası
    
    % pulse propagasyonu - verici yönü
    pulse = waveform();
    txsig = transmitter(pulse); 
    txsig = radiator(txsig,radang,radax); % pulse radiated
    txsig = txchannel(txsig,radpos,tgtpos,radvel,tgtvel);
    
    %reflection
    for i=1:TN
        [~,fwrang] = rangeangle(radpos,tgtpos(:,i),tgtax(:,:,i)); % geliş açısı
        [recvr(i),backang] = rangeangle(recvpos,tgtpos(:,i),tgtax(:,:,i)); %dönüş açısı
        tgtsig(i) = target{i}(txsig(i),fwrang,backang,tgtax(:,:,i));% hedef yansıması
    end
    % alıcı yönü
    rxsig = rxchannel(tgtsig,tgtpos,recvpos,tgtvel,recvel);
    [~,inang] = rangeangle(tgtpos,recvpos,recvax);
    
   rxsig = collector(rxsig,inang,recvax); % collector combined all incident pulses!!!
   [yc,w]= beamformer(rxsig,inang); % beamformer recalculate pulse amp.'s based on inang 
   rxpulses(:,m) = receiver(sum(yc,2)); % beamformed signal is doubled based on inang --> must add them together
end

clear yc rxsig tgtsig txsig pulse;
%% işlem öncesi plot 
figure('name','range output (before filter)');
num_pulse_plot = 2;
helperRadarPlot(rxpulses,threshold,fast_time_grid,slow_time_grid,num_pulse_plot);

%% matched filter
matchingcoeff = getMatchedFilter(waveform);
matchedfilter = phased.MatchedFilter(...
    'Coefficients',matchingcoeff,...
    'GainOutputPort',true);

[rxpulses_m, mfgain] = matchedfilter(rxpulses);
matchingdelay = size(matchingcoeff,1)-1;
rxpulses_m = buffer(rxpulses_m(matchingdelay+1:end),size(rxpulses_m,1));

threshold = threshold * db2pow(mfgain);

%helperRadarPlot(rxpulses_m,threshold,fast_time_grid,slow_time_grid,num_pulse_plot);
%% pulse integration
range_gates = prop_speed*fast_time_grid;%range 2x

tvg = phased.TimeVaryingGain(...
    'RangeLoss',fspl(range_gates,lambda),...
    'ReferenceLoss',fspl(2.5*Rmax,lambda));

rxpulses_t = tvg(rxpulses_m);
clear rxpulses_m;
rxpulses_t = pulsint(rxpulses_t,'noncoherent'); % pulse integration

figure('name','range output');
helperRadarPlot(rxpulses_t,threshold,fast_time_grid,slow_time_grid,1);
%% Range hesabı
[~,range_detect] = findpeaks(rxpulses_t,'MinPeakHeight',sqrt(threshold));
range_estimates = range_gates(range_detect);
estimates = round(Bistatik_range_calculation(flip(range_estimates),inang,Lbase))
true_range = round(recvr)
error= estimates-true_range;
%% doppler response
[resp,rng_grid,dop_grid]=dopplerResponse(rxpulses,matchingcoeff);
hPlots.himg= phased.scopes.MatrixViewer( ...
    'XStart',  dop_grid(1)*17/5, ...
    'XScale',  (dop_grid(2)-dop_grid(1))*17/5, ...
    'YStart',  rng_grid(1)/2000, ...
    'YScale',  (rng_grid(2)-rng_grid(1))/2000, ...
    'YInvert', true, ...
    'XLabel',  'Speed (km/h)', ...
    'YLabel',  'Range (km)', ...
    'Title',   'Range-Doppler Map');
hPlots.himg(mag2db(abs(resp)));
drawnow;

%% doppler 3
% angR = pi*inang/180;
% angT = pi*radang/180;
% 
% for i=1:2
%     tetaR = acos(cos(angR(1,i))*cos(angR(2,i)))-pi/2;
%     tetaT = pi/2-acos(cos(angT(1,i))*cos(angT(2,i)));
%     beta = (tetaT-tetaR)*180/pi;
% 
%     vec_delta=((Urt+Urx)/2)-tgtpos(:,i);
%     delta = atan2d(norm(cross(vec_delta,tgtvel(:,i))),dot(vec_delta,tgtvel(:,i)));
% 
%     [p, f] = periodogram(rxpulses(range_detect(2),:),[],256,prf, ...
%         'power','centered');
%     speed_vec = dop2speed(f,lambda)/2;
%     spectrum_data = p/max(p);
% 
%     [~,dop_detect] = findpeaks(pow2db(spectrum_data),'MinPeakHeight',-5);
%     speed_vec(dop_detect)
%     sp(i) = speed_vec(dop_detect)/(2*cos(delta)*cos(beta/2));
% end
% sp
