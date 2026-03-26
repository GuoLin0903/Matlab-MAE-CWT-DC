%------------------------------------------------------------------------------------------%
% ExtractFeati_OPT.m
%
% OPTIMIZED FEATURES EXTRACTION
%
% Synopsis : Feat = ExtractFeati_OPT(T2,V2,E2,heads_info, [pa])
% Extracts time, spectral and wavelets features from a waveform
%
% INPUTS :
% - T2 : Truncated Time vector of the waveform // double
% - V2 : Truncated Amplitude vector of the waveform (in volts) // double
% - E2 : Cumulative Energy of the trigger-truncated signal // double 
% - heads_info : Some headers information // cell
%       > fs : Sampling Frequency (in Hz) // double
%       > channel : Channel Number // double
%       > hit : Hit Number // double
%       > hittime : Time Of Test (in sec) // double
%
% OUTPUTS :
% - Feat : Feature vector // double
%           PF, frequency value and magnitude, signals
%
% Author : Nicolas Morizet, Postdoc Researcher
% (Adapted from Emmanuel Maillet, PhD)
% IREINE Project
% INSAVALOR, MATEIS LAB
%
% CREATION DATE : 27/06/2013
% LAST UPDATE : 18/07/2013
%------------------------------------------------------------------------------------------%

function Feat = ExtractFeati_OPT(T2,V2,E2,heads_info,varargin)

% Optional debug output folder (legacy support)
pa = '';
if ~isempty(varargin)
    pa = varargin{1};
end
% NOTE: In this repo version, debug dump lines are kept commented by default.

% Debug dump (optional): dump the exact segment used for feature extraction
%   - c_signals_EnergyCriterion.txt : one row per call, tab-separated V2
%   - c_time_EnergyCriterion.txt    : one row per call, tab-separated time (us)
% WARNING: This can become very large for big datasets. Keep cfg.features.debug_dump = false by default.
doDump = ~isempty(pa) && (ischar(pa) || isstring(pa)) && isfolder(pa);
if doDump
    try
        p1 = fullfile(char(pa), 'c_signals_EnergyCriterion.txt');
        fid = fopen(p1,'at');
        if fid ~= -1
            fprintf(fid,'%g	', V2(:));
            fprintf(fid,'
');
            fclose(fid);
        end
        p2 = fullfile(char(pa), 'c_time_EnergyCriterion.txt');
        fid = fopen(p2,'at');
        if fid ~= -1
            fprintf(fid,'%g	', T2(:)*1e6); % us
            fprintf(fid,'
');
            fclose(fid);
        end
    catch
        % ignore dump errors
    end
end
% % Output all recalculated signals
% path = fullfile(pa,'\c_signals_EnergyCriterion.txt');
% fid = fopen(path,'at');
% fprintf(fid,'%f\t',V2(:));
% fprintf(fid,'\n');
% fclose(fid);
% path = fullfile(pa,'\c_time_EnergyCriterion.txt');
% fid = fopen(path,'at');
% fprintf(fid,'%f\t',T2(:)*1e6); % unit us
% fprintf(fid,'\n');
% fclose(fid);

if (nargin < 1)
    error('At least one parameter is required !');
end

%==========================================================%
% PARAMETERS
%==========================================================%
% NEW !! : Max Frequency set to 1MHz
fMAX = 1e6; % in Hz
% NEW !! : Partial Powers Frequency Bands
% intf = 1e3*[0 100 225 400 800]'; 
  intf = 1e3*[0 100 250 500 1000]';
% intf = 1e3*[0 100 225 300 500]'; % in kHz
% Roll-off frequency factor
roll_off_factor = 0.95;
% Roll-on frequency factor
roll_on_factor = 0.05;
%==========================================================%

%==========================================================%
% RETRIEVING FS, CHANNEL, HIT & HITTIME
%==========================================================%
fs = heads_info{1,1};
channel = heads_info{1,2};
hit = heads_info{1,3};
hittime = heads_info{1,4};
%==========================================================%

%==========================================================%
% A. TIME FEATURES
%==========================================================%
% Amplitude (A, V)
[A,b] = max(abs(V2));
% Duration (D, s)
D = T2(end);
% Energy (E, V^2)
LV2 = numel(V2);

if ~isempty(E2)
    % -------------------------------------------------------------
    % 情况 1：外部传进来了累积能量向量 E2（AEwin 风格）
    %         原版代码就是取最后一点的值，这里加一下防越界保护
    % -------------------------------------------------------------
    nE2 = numel(E2);
    idx = min(max(LV2,1), nE2);   % 1 <= idx <= nE2
    E   = E2(idx);

else
    % -------------------------------------------------------------
    % 情况 2：E2 为空（这是你 AE_Main 现在的情况）
    %         直接从波形计算离散能量：
    %             E ∝ Σ V(t)^2
    %     对同一文件来说只是差一个常数因子，
    %     对聚类 / Weasel 来说不影响判别。
    % -------------------------------------------------------------
    if LV2 == 0 || all(~isfinite(V2))
        E = NaN;   % 极端异常情况，避免 NaN 扩散到相关系数里
    else
        % 版本 A：完全离散的能量（单位“V^2·sample”，更接近 AEwin 的列名）
        E = sum(V2.^2);


    end
end

% Zero-Crossings (ZC, -) & Zero-Crossing Rate (ZCR, %)
c = (V2(1:LV2-1)<=0) & (V2(2:LV2)>0);
ZC = sum(c);      
ZCR = 100*ZC/LV2;
% Rise time (RT, s)
RT = T2(b);
% Temporal centroid (TC, s)
V2_RMS = abs(V2) / sqrt(LV2); % RMS of the signal
TC = sum(T2.*V2_RMS)./sum(V2_RMS);
% Temporal decrease (alpha, -)
[~,b] = max(V2_RMS); % Computed from the RMS MAX !
p = polyfit(T2(b+1:end),V2_RMS(b+1:end),1);
alpha = -p(1);
%==========================================================%

%==========================================================%
% B. FREQUENCY FEATURES
%==========================================================%
% Regular FFT (preferable to find a good peak frequency)
NFFT = 2^nextpow2(LV2);
Y = fft(V2,NFFT)/LV2;
f = fs/2*linspace(0,1,NFFT/2+1); % watch out when SPI resampling is used here !!..

% OUTPUT PF in low frequency and high frequency
% Modified by Xi CHEN
fMAX1 = 0.1e6; % 100 kHz in Hz
FcL = f((f<=fMAX1)&(0<f))';
YL = abs(Y((f<=fMAX1)&(0<f)))';
PFL = FcL(YL == max(YL))/1e3; % peak frequency in kHz
fMAX2 = 1e6; % in Hz
FcH = f((f<=fMAX2)&(fMAX1<f))';
YH = abs(Y((f<=fMAX2)&(fMAX1<f)))';
PFH = FcH(YH == max(YH))/1e3; % peak frequency in kHz
% freqout = [PFL PFH max(YL)/max(YH)];
% path = fullfile(pa,'\PF_pics.txt');
% fid = fopen(path,'at');
% fprintf(fid,'%f\t',freqout(:));
% fprintf(fid,'\n');
% fclose(fid);

% Focusing on useful frequencies
Fc = f(f<=fMAX);
Y = abs(Y(f<=fMAX))';

% % OUTPUT frequency map amplitude and frequency Hz
% % Modified by Xi CHEN
% path = fullfile(pa,'\frequency_value.txt');
% fid = fopen(path,'at');
% fprintf(fid,'%f\t',Fc(:));
% fprintf(fid,'\n');
% fclose(fid);
% path = fullfile(pa,'\frequency_magnitude.txt');
% fid = fopen(path,'at');
% fprintf(fid,'%f\t',Y(:));
% fprintf(fid,'\n');
% fclose(fid);

% Partial Powers (PP1 to PP4, %)
Ptot = sum(Y);
PP2 = zeros(length(intf)-1,1);
for i=1:length(intf)-1
    PP2(i) = 100*sum(Y(Fc>=intf(i) & Fc<intf(i+1))) / Ptot;
end

% Frequency centroid (FC2, Hz)
FC2 = sum(Fc.*Y)/sum(Y);
% Peak frequency (PF2, Hz)
PF2 = Fc(Y == max(Y));
% Spectral spread (SSpread, Hz)
SSpread = sqrt(sum(((Fc-FC2).^2).*Y)/sum(Y));
% Spectral skewness (SSkew, -)
SSkew = (1/SSpread^(3))*(sum(((Fc-FC2).^3).*Y)/sum(Y));
% Spectral kurtosis (SKurt, -)
SKurt = (1/SSpread^(4))*(sum(((Fc-FC2).^4).*Y)/sum(Y));
% Spectral slope (SSlope, -) [linear regression of the normalized spectral amplitude]
p = polyfit(Fc/fMAX,Y/max(Y),1);
SSlope = p(1);
% Roll-off frequency (SRoff, Hz)
cum_energy = cumsum(Y);
R = roll_off_factor * max(cum_energy);
indroff = find(cum_energy < R, 1, 'last');
SRoff = Fc(indroff);
% Spectral spread to peak (SSpreadP, Hz)
SSpreadP = sqrt(sum(((Fc-PF2).^2).*Y)/sum(Y));
% Spectral skewness to peak (SSkewP)
SSkewP = (1/SSpreadP^(3))*(sum(((Fc-PF2).^3).*Y)/sum(Y));
% Spectral kurtosis to peak (SKurtP)
SKurtP = (1/SSpreadP^(4))*(sum(((Fc-PF2).^4).*Y)/sum(Y));
% Roll-on frequency (SRon, Hz)
R = roll_on_factor * max(cum_energy);
indron = find(cum_energy < R, 1, 'last');
if isempty(indron)
    indron = 1;
end
SRon = Fc(indron);
%==========================================================%
%==========================================================%
%% Convert Unit Data
% Amplitude in dB
A = 20.*(log10(A.*1e6))-42;
% Absolute Energy in aJ
% E = 4*E./(D.*10);
% E = 1024*E.*D;   % Prise en compte du resampling 
% % Duration, Rise Time and Temporal Centroid in Ä¾s
D  = (D.*10^(6));
RT = (RT.*10^(6));
TC = (TC.*10^(6));
% Frequencies in kHz
FC2      = FC2./1000;
PF2      = PF2./1000;
SRoff    = SRoff./1000;
SRon     = SRon./1000;
SSpread  = SSpread./1000;
SSpreadP = SSpreadP./1000;
% % Slope in kHz-1
% SSlope = abs(SSlope.*1000);

%==========================================================%
% C. WAVELET PACKETS FEATURES
%==========================================================%
% Wavelet Packets Energy (WPE1 to WPE8, %)
Twp = wpdec(V2,3,'sym8','shannon');
Ewp = wenergy(Twp);
%Shannon entropy criteria
[N,~] = histcounts(abs(V2),100);
L = N(N~=0);
L = L./size(V2,1);
Entropy = sum(L.*log2(L))*(-1);
%==========================================================%

%==========================================================%
%% Ratio features (from M. Moevus)
% Tr  = RT./D;      % Temps de montÃ©e relatif
% D_A = D./A;       % DurÃ©e / Amplitude
% T_E = D-RT;       % Temps d'extinction
% A_M = A./RT;      % Angle de montÃ©e
% A_D = A./(D-RT);  % Temps de montÃ©e sur descente
% E_R = E./A;       % Energie relative
% A_F = A./FC2;     % Amplitude / frÃ©quence centroÄ�de
% A_E = A./Entropy; % Amplitude / Entropie

%==========================================================%
% D. BUILDING FINAL FEATURES VECTOR
%==========================================================%
% Feat = [hit hittime channel A D E ZCR RT TC alpha,...
%         PP2' FC2 PF2 SSpread SSkew SKurt SSlope SRoff,...
%         SSpreadP SSkewP SKurtP SRon Tr D_A T_E A_M A_D,...
%         E_R A_F A_E Entropy]';
    
     
Feat = [hit hittime channel A D E ZCR RT TC alpha,...
      PP2' FC2 PF2 SSpread SSkew SKurt SSlope SRoff,...
        sqrt(SSpreadP) SSkewP SKurtP SRon Ewp Entropy]';
    

%==========================================================%


