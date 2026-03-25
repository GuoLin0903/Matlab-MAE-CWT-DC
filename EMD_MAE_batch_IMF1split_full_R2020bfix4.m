%% ================== EMD_MAE_batch_IMF1split_full.m ==================
% Batch EMD-MAE with dispersion-guided IMF1 splitting (exclusive masking)
% For your MATLAB R2020b workflow (waveReader + CWT + EMD)
%
% What this script does:
%  1) Multi-select *.wave/*.tradb (Hexagon/DigitalWave) via uigetfile
%  2) For each file: choose channel; choose events (multi-select if 3D)
%  3) For each event: crop a fixed window; run EMD; compute CWT; overlay dispersion curves
%  4) A0/S0 reconstruction:
%     - IMF2..K: grouped by low-band CWT energy ratio (old rule)
%     - IMF1: if dispersion curves available => split into IMF1(A0-part) + IMF1(S0-part)
%            using EXCLUSIVE time-frequency masks around A0/S0 curves
%            (with time gate to avoid reflections)
%     - If IMF1 split fails => FALLBACK to the old rule (no forced IMF1->S0)
%  5) Save figures + results.csv + results.mat
%
% Notes:
%  - This is a signal-processing heuristic (dispersion-guided masking), not a proof of physical purity.
%  - EMD mode-mixing is expected for transient multi-packet guided waves (Huang et al., 1998).
%
% ---------------------------------------------------------------------
% Dependencies:
%  - waveReader.m on path (yours)
%  - Signal Processing Toolbox: emd
%  - Wavelet Toolbox: cwt, icwt
%
% Author: generated for GuoLin0903/Matlab-MAE-CWT-DC workflow
% =====================================================================

clear; clc; close all;

%% ================= 0) USER OPTIONS =================
% -------- Output --------
outOpt.save_dir = fullfile(pwd, ['EMD_MAE_batch_IMF1split_out_', datestr(now,'yyyymmdd_HHMMSS')]);
outOpt.save_png = true;
outOpt.save_pdf = false;
outOpt.dpi      = 300;
outOpt.show_fig = false;          % false for batch

% -------- Window crop --------
dataOpt.start_us  = 0;
dataOpt.sigLen_us = 500;          % 500 us @ 5 MHz => 2500 samples
dataOpt.demean    = true;

% -------- EMD params (compatible wrapper) --------
emdOpt.max_imf = 5;
emdOpt.siftRelTol = 0.1;
emdOpt.interp = 'PCHIP';
emdOpt.maxEnergyRatio = 20;
emdOpt.siftMaxIter = 100;
emdOpt.maxNumExtrema = 1;
emdOpt.display = 0;
emdOpt.use_default_stop = false;

% -------- CWT / display --------
dispOpt.fmax_kHz = 400;
dispOpt.nFreqBins = 400;
dispOpt.voicesPerOctave = 32;
dispOpt.threshold = 0.01;
dispOpt.gamma = 1;
dispOpt.force_same_rawrow = true;

% -------- A0/S0 grouping (IMF2..K) --------
sepOpt.split_freq_Hz = 150e3;     % low/high boundary for lowRatio
sepOpt.ratio_th = 0.50;

% -------- IMF1 split (dispersion-guided exclusive masking) --------
sepOpt.split_imf1_enable   = true;
sepOpt.imf1_sigma_us       = 10;      % curve tolerance
sepOpt.imf1_mag_th         = 0.20;    % coefficient magnitude gate
sepOpt.imf1_delta_us       = 12;      % exclusivity margin (avoid ambiguous region)
sepOpt.imf1_time_gate_us   = [80 260];% only split in this time window (avoid reflections)
sepOpt.imf1_show_parts_in_fig1 = true;

% -------- Dispersion curve overlay --------
curveOpt.enable = true;
curveOpt.dist_mm_default = 100;
curveOpt.shift_us = 66;             % global shift for curve alignment
curveOpt.lineWidth = 1.6;
curveOpt.A0_style = '-';
curveOpt.S0_style = '--';
curveOpt.color = [1 1 1];           % white

% -------- waveReader headerlength --------
ioOpt.try_headerlength = true;
ioOpt.headerlength = 502;

%% ================= 1) Select wave files =================
fprintf('请选择数据文件（*.wave / *.tradb）... 可多选\n');
[fname, pname] = uigetfile({'*.wave;*.tradb;*.*'}, '选择 wave 文件（可多选）', 'MultiSelect','on');
if isequal(fname, 0), error('未选择 wave 文件'); end
if ischar(fname), fname = {fname}; end
nFiles = numel(fname);

if ~exist(outOpt.save_dir,'dir'); mkdir(outOpt.save_dir); end
figDir = fullfile(outOpt.save_dir, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

%% ================= 2) Select dispersion curve file (optional) =================
curve_raw = [];
curve_file = '';
if curveOpt.enable
    fprintf('请选择色散曲线文件（Vallen 导出 Excel/CSV）... 可取消\n');
    [cname, cpname] = uigetfile({'*.xlsx;*.xls;*.csv;*.txt;*.*'}, '选择色散曲线文件（可取消）');
    if isequal(cname, 0)
        warning('未选择曲线文件：将不叠加色散曲线，也不会启用 IMF1 split。');
        curveOpt.enable = false;
    else
        curve_file = fullfile(cpname, cname);
        try
            curve_raw = readcell(curve_file);
        catch ME
            warning('曲线文件 readcell 失败：%s\n将不叠加曲线。', ME.message);
            curveOpt.enable = false;
        end
    end
end

%% ================= 3) Batch loop =================
rows = {};
head = {'file','distance_mm','fs_Hz','ch','event','Nwin','K_imf', ...
        'idx_A0_base','idx_S0_base','imf1_split_ok', ...
        'imf1_A0_energy_frac','imf1_S0_energy_frac', ...
        'lowRatio_all','A0_lowRatio_mean','S0_lowRatio_mean'};

for fi = 1:nFiles
    wavefile = fullfile(pname, fname{fi});
    fprintf('\n==================== [%d/%d] %s ====================\n', fi, nFiles, fname{fi});

    [wave, fs] = readWaveRobust_plus(wavefile, ioOpt);
    wave = double(wave); fs = double(fs);

    [chSel, evList, sigList] = pickMultiEvents_from_wave(wave);
    nEv = numel(evList);

    dist_mm = parse_distance_mm(fname{fi});
    if isnan(dist_mm), dist_mm = curveOpt.dist_mm_default; end

    curve_data = struct('ok',false,'msg','','A0',struct(),'S0',struct());
    if curveOpt.enable && ~isempty(curve_raw)
        [ok, A0, S0, msg] = parse_vallen_like_cell(curve_raw, dist_mm, curveOpt.shift_us);
        if ok
            curve_data.ok = true; curve_data.A0 = A0; curve_data.S0 = S0; curve_data.msg = msg;
        else
            curve_data = load_dispersion_curve(curve_file, dist_mm, curveOpt.shift_us);
        end
        if ~curve_data.ok
            warning('曲线解析失败（dist=%gmm）：%s；本文件将不叠加曲线/不 split IMF1。', dist_mm, curve_data.msg);
        end
    end

    for i = 1:nEv
        evSel = evList(i);
        raw_full = sigList{i};
        raw_full(~isfinite(raw_full)) = 0;
        raw_full = raw_full(:);
        if dataOpt.demean, raw_full = raw_full - mean(raw_full); end

        start_idx = max(1, round(dataOpt.start_us*1e-6*fs) + 1);
        if isinf(dataOpt.sigLen_us)
            end_idx = numel(raw_full);
        else
            end_idx = min(numel(raw_full), start_idx + round(dataOpt.sigLen_us*1e-6*fs) - 1);
        end
        if end_idx <= start_idx
            warning('file=%s event=%d: 截取范围无效，跳过。', fname{fi}, evSel);
            continue;
        end
        raw_sig = raw_full(start_idx:end_idx);
        N = numel(raw_sig);
        t_us = (0:N-1)/fs*1e6;

        one = run_one_event_emd_mae_IMF1split(raw_sig, t_us, fs, curve_data, curveOpt, emdOpt, dispOpt, sepOpt);

        base = sprintf('%s__d%gmm__ch%d__ev%04d', strip_ext(fname{fi}), dist_mm, chSel, evSel);
        if outOpt.show_fig, set(0,'DefaultFigureVisible','on'); else, set(0,'DefaultFigureVisible','off'); end

        % Fig1: decomposition (optionally includes IMF1 parts)
        plot_dual_view_AGU_curve(one.list_decomp, t_us, fs, dispOpt, ...
            sprintf('Figure 1: EMD Decomposition | %s', base), ...
            one.ref_max_raw, one.raw_ylim, curve_data, curveOpt);
        fig1 = gcf; set(fig1,'Renderer','opengl');

        % Fig2: recon
        plot_dual_view_AGU_curve(one.list_recon, t_us, fs, dispOpt, ...
            sprintf('Figure 2: A0/S0 Reconstruction | %s', base), ...
            one.ref_max_raw, one.raw_ylim, curve_data, curveOpt);
        fig2 = gcf; set(fig2,'Renderer','opengl');

        if outOpt.save_png
            save_fig(fig1, fullfile(figDir, [base,'__Fig1_IMFs.png']), outOpt.dpi);
            save_fig(fig2, fullfile(figDir, [base,'__Fig2_Recon.png']), outOpt.dpi);
        end
        if outOpt.save_pdf
            save_fig_pdf(fig1, fullfile(figDir, [base,'__Fig1_IMFs.pdf']));
            save_fig_pdf(fig2, fullfile(figDir, [base,'__Fig2_Recon.pdf']));
        end

        if ~outOpt.show_fig, close(fig1); close(fig2); end

        rows(end+1,:) = {fname{fi}, dist_mm, fs, chSel, evSel, N, ...
            one.K, mat2str(one.idx_A0_base), mat2str(one.idx_S0_base), one.imf1_split_ok, ...
            one.imf1_A0_energy_frac, one.imf1_S0_energy_frac, ...
            mat2str(one.lowRatio,3), one.A0_lowRatio_mean, one.S0_lowRatio_mean}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', head);
writetable(T, fullfile(outOpt.save_dir, 'results.csv'));
save(fullfile(outOpt.save_dir, 'results.mat'), 'T', 'outOpt','dataOpt','emdOpt','dispOpt','sepOpt','curveOpt');

fprintf('\n==================== DONE ====================\n');
fprintf('Output: %s\n', outOpt.save_dir);
disp(T);

%% =====================================================================
%% ============================ MAIN EVENT ==============================
function one = run_one_event_emd_mae_IMF1split(raw_sig, t_us, fs, curve_data, curveOpt, emdOpt, dispOpt, sepOpt)

    [imfs, ~] = emd_compat(raw_sig, emdOpt);
    if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
    K = size(imfs,2);

    [pkHz, cHz] = estimate_mode_freqs(imfs, fs);

    % lowRatio for all IMFs (used for IMF2..K grouping)
    lowRatio = zeros(1,K);
    for k = 1:K
        lowRatio(k) = lowband_energy_ratio_cwt(imfs(:,k), fs, sepOpt.split_freq_Hz, dispOpt.voicesPerOctave);
    end
    idx_A0 = find(lowRatio >= sepOpt.ratio_th);
    idx_S0 = setdiff(1:K, idx_A0);

    % protect empty groups
    if isempty(idx_A0)
        [~, imax] = max(lowRatio); idx_A0 = imax; idx_S0 = setdiff(1:K, idx_A0);
    end
    if isempty(idx_S0)
        [~, imin] = min(lowRatio); idx_S0 = imin; idx_A0 = setdiff(1:K, idx_S0);
    end

    % -------- IMF1 split (exclusive masks). If fails => fallback to old logic ----------
    imf1_split_ok = false;
    imf1_A0 = zeros(size(raw_sig));
    imf1_S0 = zeros(size(raw_sig));
    imf1_A0_energy_frac = NaN;
    imf1_S0_energy_frac = NaN;

    if sepOpt.split_imf1_enable && curveOpt.enable && isfield(curve_data,'ok') && curve_data.ok && K>=1
        [imf1_A0, imf1_S0, imf1_split_ok] = split_by_dispersion_mask_imf1( ...
            imfs(:,1), fs, curve_data, ...
            sepOpt.imf1_sigma_us, sepOpt.imf1_mag_th, sepOpt.imf1_delta_us, sepOpt.imf1_time_gate_us, ...
            dispOpt.voicesPerOctave);

        if imf1_split_ok
            e0 = sum(imfs(:,1).^2) + eps;
            imf1_A0_energy_frac = sum(imf1_A0.^2) / e0;
            imf1_S0_energy_frac = sum(imf1_S0.^2) / e0;
        end
    end

    if imf1_split_ok
        idx_A0_base = idx_A0; idx_S0_base = idx_S0;
        idx_A0_base(idx_A0_base==1) = [];
        idx_S0_base(idx_S0_base==1) = [];

        A0_recon = sum(imfs(:, idx_A0_base), 2) + imf1_A0;
        S0_recon = sum(imfs(:, idx_S0_base), 2) + imf1_S0;
    else
        idx_A0_base = idx_A0;
        idx_S0_base = idx_S0;
        A0_recon = sum(imfs(:, idx_A0), 2);
        S0_recon = sum(imfs(:, idx_S0), 2);
    end

    % -------- plot lists --------
    list_decomp = {};
    list_decomp{1}.sig = raw_sig;
    list_decomp{1}.name = 'Original signal (mix)';

    for k = 1:K
        list_decomp{end+1}.sig  = imfs(:,k);
        list_decomp{end}.name   = sprintf('IMF %d (centroid %.1f kHz, peak %.1f kHz)', k, cHz(k)/1000, pkHz(k)/1000);

        if k==1 && imf1_split_ok && sepOpt.imf1_show_parts_in_fig1
            list_decomp{end+1}.sig = imf1_S0;
            list_decomp{end}.name  = sprintf('IMF1 -> S0-part (mask, σ=%.1fµs, δ=%.1fµs)', sepOpt.imf1_sigma_us, sepOpt.imf1_delta_us);
            list_decomp{end+1}.sig = imf1_A0;
            list_decomp{end}.name  = sprintf('IMF1 -> A0-part (mask, σ=%.1fµs, δ=%.1fµs)', sepOpt.imf1_sigma_us, sepOpt.imf1_delta_us);
        end
    end

    list_recon = cell(1,3);
    list_recon{1}.sig = raw_sig;
    list_recon{1}.name = 'Original signal (mix)';
    if imf1_split_ok
        list_recon{2}.sig = S0_recon;
        list_recon{2}.name = sprintf('S0 reconstruction (base IMFs %s + IMF1(S0-part))', mat2str(idx_S0_base));
        list_recon{3}.sig = A0_recon;
        list_recon{3}.name = sprintf('A0 reconstruction (base IMFs %s + IMF1(A0-part))', mat2str(idx_A0_base));
    else
        list_recon{2}.sig = S0_recon;
        list_recon{2}.name = sprintf('S0 reconstruction (IMFs %s)', mat2str(idx_S0));
        list_recon{3}.sig = A0_recon;
        list_recon{3}.name = sprintf('A0 reconstruction (IMFs %s)', mat2str(idx_A0));
    end

    [ref_max_raw, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt);

    one = struct();
    one.K = K;
    one.lowRatio = lowRatio;
    one.idx_A0_base = idx_A0_base;
    one.idx_S0_base = idx_S0_base;
    one.imf1_split_ok = imf1_split_ok;
    one.imf1_A0_energy_frac = imf1_A0_energy_frac;
    one.imf1_S0_energy_frac = imf1_S0_energy_frac;
    one.A0_lowRatio_mean = mean(lowRatio(idx_A0),'omitnan');
    one.S0_lowRatio_mean = mean(lowRatio(idx_S0),'omitnan');
    one.list_decomp = list_decomp;
    one.list_recon  = list_recon;
    one.ref_max_raw = ref_max_raw;
    one.raw_ylim    = raw_ylim;
end

%% =====================================================================
%% ============================ I/O & GUI ===============================
function [wave, fs] = readWaveRobust_plus(wavefile, ioOpt)
    wave = []; fs = [];
    % 13 outputs
    try
        out = cell(1,13);
        [out{:}] = waveReader(wavefile);
        wave = out{2}; fs = out{7};
        if ~isempty(fs), return; end
    catch
    end
    % 7 outputs
    try
        out = cell(1,7);
        [out{:}] = waveReader(wavefile);
        wave = out{2}; fs = out{7};
        if ~isempty(fs), return; end
    catch
    end
    if ioOpt.try_headerlength
        try
            out = cell(1,13);
            [out{:}] = waveReader(wavefile, ioOpt.headerlength);
            wave = out{2}; fs = out{7};
            if ~isempty(fs), return; end
        catch
        end
        out = cell(1,7);
        [out{:}] = waveReader(wavefile, ioOpt.headerlength);
        wave = out{2}; fs = out{7};
        if ~isempty(fs), return; end
    end
    error('waveReader failed for file: %s', wavefile);
end

function [chSel, evList, sigList] = pickMultiEvents_from_wave(wave)
    chSel = 1; sigList = {};
    if isvector(wave)
        evList = 1; sigList{1} = wave(:); return;
    end
    nd = ndims(wave);
    if nd == 2
        [~,C] = size(wave);
        if C > 1
            listC = compose('ch %d', 1:C);
            [chSel, ok] = listdlg('ListString',listC,'SelectionMode','single', ...
                'PromptString','选择通道 channel（单选）', 'InitialValue',1,'ListSize',[240 280]);
            if isempty(ok) || ok==0, chSel = 1; end
        end
        evList = 1;
        sigList{1} = wave(:,chSel);
        return;
    end
    C  = size(wave,2);
    Ev = size(wave,3);
    if C > 1
        listC = compose('ch %d', 1:C);
        [chSel, ok] = listdlg('ListString',listC,'SelectionMode','single', ...
            'PromptString','选择通道 channel（单选）', 'InitialValue',1,'ListSize',[240 280]);
        if isempty(ok) || ok==0, chSel = 1; end
    end
    if Ev > 1
        listE = compose('event %d', 1:Ev);
        [evList, ok] = listdlg('ListString',listE,'SelectionMode','multiple', ...
            'PromptString','选择事件 event（可多选：Ctrl/Shift）', ...
            'InitialValue',1:min(Ev,10),'ListSize',[260 360]);
        if isempty(ok) || ok==0, evList = 1; end
    else
        evList = 1;
    end
    sigList = cell(1,numel(evList));
    for i = 1:numel(evList)
        sigList{i} = wave(:, chSel, evList(i));
    end
end

function dist_mm = parse_distance_mm(fname)
    dist_mm = NaN;
    tok = regexp(fname,'(\d+)\s*mm','tokens','once');
    if ~isempty(tok)
        dist_mm = str2double(tok{1});
    end
end

function s = strip_ext(fname)
    [~,s,~] = fileparts(fname);
end

%% =====================================================================
%% ============================ SAVE FIG ================================
function save_fig(figH, outPng, dpi)
    try
        exportgraphics(figH, outPng, 'Resolution', dpi);
    catch
        print(figH, outPng, '-dpng', sprintf('-r%d', dpi));
    end
end

function save_fig_pdf(figH, outPdf)
    try
        exportgraphics(figH, outPdf, 'ContentType','vector');
    catch
        print(figH, outPdf, '-dpdf', '-painters');
    end
end

%% =====================================================================
%% ============================ EMD & FEATURES ==========================
function [imfs, res] = emd_compat(x, emdOpt)
    x = x(:);
    if exist('emd','file') ~= 2
        error('emd() not found. Install Signal Processing Toolbox.');
    end
    if emdOpt.use_default_stop
        try
            [imfs, res] = emd(x, 'MaxNumIMF', emdOpt.max_imf);
        catch
            [imfs, res] = emd(x);
            if size(imfs,2) > emdOpt.max_imf, imfs = imfs(:,1:emdOpt.max_imf); end
        end
        return;
    end
    try
        [imfs, res] = emd(x, ...
            'MaxNumIMF', emdOpt.max_imf, ...
            'SiftRelativeTolerance', emdOpt.siftRelTol, ...
            'Interpolation', emdOpt.interp, ...
            'MaxEnergyRatio', emdOpt.maxEnergyRatio, ...
            'SiftMaxIterations', emdOpt.siftMaxIter, ...
            'MaxNumExtrema', emdOpt.maxNumExtrema, ...
            'Display', emdOpt.display);
    catch
        try
            [imfs, res] = emd(x, ...
                'MaxNumIMF', emdOpt.max_imf, ...
                'SiftRelativeTolerance', emdOpt.siftRelTol, ...
                'Interpolation', emdOpt.interp, ...
                'MaxEnergyRatio', emdOpt.maxEnergyRatio, ...
                'Display', emdOpt.display);
        catch
            [imfs, res] = emd(x, 'MaxNumIMF', emdOpt.max_imf);
        end
    end
end

function [pkHz, cHz] = estimate_mode_freqs(modes, fs)
    if size(modes,1) < size(modes,2), modes = modes.'; end
    K = size(modes,2);
    pkHz = zeros(1,K);
    cHz  = zeros(1,K);
    N = size(modes,1);
    nfft = 2^nextpow2(N);
    f = (0:nfft-1)/nfft*fs;
    half = 1:floor(nfft/2);
    for k = 1:K
        x = modes(:,k);
        X = abs(fft(x, nfft)).^2;
        P = X(half);
        ff = f(half);
        if all(P==0), pkHz(k)=0; cHz(k)=0; continue; end
        [~, imax] = max(P);
        pkHz(k) = ff(imax);
        cHz(k)  = sum(ff(:).*P(:)) / (sum(P(:))+eps);
    end
end

function rLow = lowband_energy_ratio_cwt(x, fs, split_freq_Hz, voicesPerOctave)
    [cfs, f] = cwt(x, fs, 'amor', 'VoicesPerOctave', voicesPerOctave);
    mag2 = abs(cfs).^2;
    if f(1) > f(end)
        f = flipud(f); mag2 = flipud(mag2);
    end
    low = (f <= split_freq_Hz);
    Elow = sum(mag2(low,:), 'all');
    Etot = sum(mag2, 'all');
    rLow = Elow / (Etot + eps);
end

function [ref_max, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt)
    [cfs, ~] = cwt(raw_sig, fs, 'amor', 'VoicesPerOctave', dispOpt.voicesPerOctave);
    ref_max = max(abs(cfs(:)));
    if ref_max == 0, ref_max = 1; end
    ymax = max(abs(raw_sig(:)));
    raw_ylim = [-1.05*ymax, 1.05*ymax];
end

%% =====================================================================
%% ============================ DISPERSION CURVES =======================
function curve = load_dispersion_curve(curvefile, dist_mm, shift_us)
    curve = struct('ok',false,'msg','', 'A0',struct(),'S0',struct());
    if ~isfile(curvefile), curve.msg = 'curve file not found'; return; end
    try
        raw = readcell(curvefile);
        [ok, A0, S0, msg] = parse_vallen_like_cell(raw, dist_mm, shift_us);
        if ok
            curve.ok = true; curve.A0 = A0; curve.S0 = S0; curve.msg = msg;
            return;
        end
    catch
    end
    curve.msg = 'failed to parse dispersion curve';
end

function [ok, A0, S0, msg] = parse_vallen_like_cell(raw, dist_mm, shift_us)
    ok = false; msg = '';
    A0 = struct(); S0 = struct();
    if isempty(raw) || size(raw,1) < 10, msg = 'raw too small'; return; end
    [rS0, cS0] = find_cell(raw, 'Group S0');
    [rA0, cA0] = find_cell(raw, 'Group A0');
    if isempty(rS0) || isempty(rA0), msg = 'Group S0/A0 not found'; return; end
    r_data = rS0(1) + 2;
    fS = cell2num(raw(r_data:end, cS0(1)));
    vS = cell2num(raw(r_data:end, cS0(1)+1));
    fA = cell2num(raw(r_data:end, cA0(1)));
    vA = cell2num(raw(r_data:end, cA0(1)+1));
    [fS, vS] = trim_nan_pairs(fS, vS);
    [fA, vA] = trim_nan_pairs(fA, vA);
    if numel(fS) < 5 || numel(fA) < 5, msg = 'too few points'; return; end
    S0.f_kHz = fS * 1000;
    A0.f_kHz = fA * 1000;
    vS_mm_per_us = vS; % m/ms == mm/us
    vA_mm_per_us = vA;
    S0.t_us = dist_mm ./ (vS_mm_per_us + eps) + shift_us;
    A0.t_us = dist_mm ./ (vA_mm_per_us + eps) + shift_us;
    S0.dist_mm = dist_mm; A0.dist_mm = dist_mm;
    S0.shift_us = shift_us; A0.shift_us = shift_us;
    ok = true;
    msg = 'parsed Vallen-like dispersion export';
end

function [r, c] = find_cell(raw, key)
    r = []; c = [];
    key = lower(string(key));
    for i = 1:size(raw,1)
        for j = 1:size(raw,2)
            v = raw{i,j};
            if ischar(v) || isstring(v)
                s = lower(string(v));
                if contains(s, key)
                    r(end+1) = i; %#ok<AGROW>
                    c(end+1) = j; %#ok<AGROW>
                end
            end
        end
    end
end

function x = cell2num(col)
    x = nan(numel(col),1);
    for i = 1:numel(col)
        v = col{i};
        if isnumeric(v) && isfinite(v)
            x(i) = double(v);
        elseif ischar(v) || isstring(v)
            vv = str2double(string(v));
            if isfinite(vv), x(i) = vv; end
        end
    end
end

function [a,b] = trim_nan_pairs(a,b)
    a = a(:); b = b(:);
    ok = isfinite(a) & isfinite(b);
    a = a(ok); b = b(ok);
end

%% =====================================================================
%% ============================ IMF1 SPLIT (EXCLUSIVE MASK) =============
function [xA0, xS0, ok] = split_by_dispersion_mask_imf1(x, fs, curve_data, sigma_us, mag_th, delta_us, time_gate_us, voicesPerOctave)
    ok = false;
    x = x(:);
    N = numel(x);
    t_us = (0:N-1)/fs*1e6;

    [fA, tA] = prep_curve(curve_data.A0.f_kHz(:), curve_data.A0.t_us(:));
    [fS, tS] = prep_curve(curve_data.S0.f_kHz(:), curve_data.S0.t_us(:));

    try
        [wt, f] = cwt(x, fs, 'amor', 'VoicesPerOctave', voicesPerOctave);
    catch
        [wt, f] = cwt(x, fs, 'amor');
    end
    % R2020b icwt() expects f to be strictly decreasing.
    % cwt() may return increasing or decreasing; enforce decreasing here.
    if f(1) < f(end)
        f = flipud(f);
        wt = flipud(wt);
    end

    % ---- Enforce strictly decreasing frequency vector for R2020b icwt() ----
    % cwt() may return non-strict monotonic f due to rounding; sort + reorder rows.
    [f, idxSort] = sort(f(:), 'descend');
    wt = wt(idxSort, :);
    df = diff(f);
    if any(df >= 0)
        % Make strictly decreasing by subtracting a tiny ramp (preserves order).
        tiny = eps(max(f));
        f = f - (0:numel(f)-1)' * tiny;
    end
    f_kHz = f/1e3;

    tAq = interp1(fA, tA, f_kHz, 'linear', 'extrap');
    tSq = interp1(fS, tS, f_kHz, 'linear', 'extrap');

    Tg = repmat(t_us, numel(f_kHz), 1);
    dA = abs(Tg - tAq);
    dS = abs(Tg - tSq);

    A = abs(wt);
    A = A ./ (max(A(:)) + eps);
    gate = (A >= mag_th);

    if ~isempty(time_gate_us) && numel(time_gate_us)==2
        tg = (t_us >= time_gate_us(1)) & (t_us <= time_gate_us(2));
        gate = gate & repmat(tg, numel(f_kHz), 1);
    end

    maskA = gate & (dA <= sigma_us) & ((dA + delta_us) <= dS);
    maskS = gate & (dS <= sigma_us) & ((dS + delta_us) <= dA);

    if nnz(maskA) < 200 || nnz(maskS) < 200
        xA0 = zeros(N,1);
        xS0 = zeros(N,1);
        ok = false;
        return;
    end

freqrange = [min(f) max(f)];
    try
        % Some MATLAB versions expect icwt(wt, freqrange, f, ...)
        xA0 = icwt(wt .* maskA, freqrange, f, 'SignalMean', 0);
        xS0 = icwt(wt .* maskS, freqrange, f, 'SignalMean', 0);
    catch
        % R2020b often expects icwt(wt, f, freqrange, ...) and requires f strictly decreasing
        xA0 = icwt(wt .* maskA, f, freqrange, 'SignalMean', 0);
        xS0 = icwt(wt .* maskS, f, freqrange, 'SignalMean', 0);
    end
xA0 = xA0(:); xS0 = xS0(:);
    ok = true;
end

function [f_sorted, t_sorted] = prep_curve(f_kHz, t_us)
    f_kHz = f_kHz(:); t_us = t_us(:);
    good = isfinite(f_kHz) & isfinite(t_us) & (f_kHz>0);
    f_kHz = f_kHz(good); t_us = t_us(good);
    [f_kHz, ord] = sort(f_kHz);
    t_us = t_us(ord);
    [f_sorted, ia] = unique(f_kHz, 'stable');
    t_sorted = t_us(ia);
end

%% =====================================================================
%% ============================ PLOTTING ================================
function plot_dual_view_AGU_curve(plotList, t_us, fs, dispOpt, figTitle, ref_max_raw, raw_ylim, curve_data, curveOpt)
    nRows = numel(plotList);
    f_lin = linspace(0, dispOpt.fmax_kHz*1000, dispOpt.nFreqBins);

    stored_cwt = cell(1,nRows);
    all_max = 0;
    for i = 1:nRows
        try
            [cfs, f] = cwt(plotList{i}.sig, fs, 'amor', 'VoicesPerOctave', dispOpt.voicesPerOctave);
        catch
            [cfs, f] = cwt(plotList{i}.sig, fs, 'amor');
        end
        mag = abs(cfs);
        if f(1) > f(end), f = flipud(f); mag = flipud(mag); end
        mag_lin = interp1(f, mag, f_lin, 'linear', 0);
        stored_cwt{i} = mag_lin;
        all_max = max(all_max, max(mag_lin(:)));
    end
    if all_max == 0, all_max = 1; end

    fig_h = min(1100, 170*nRows);
    figure('Name', figTitle, 'Color','w', 'Position', [60, 60, 1250, fig_h]);
    tl = tiledlayout(nRows, 2, 'TileSpacing','compact', 'Padding','compact');
    f_kHz_axis = f_lin/1000;

    for i = 1:nRows
        nexttile;
        plot(t_us, plotList{i}.sig, 'b', 'LineWidth', 1);
        grid on; xlim([0, t_us(end)]);
        ylabel('Amp');
        title(plotList{i}.name, 'FontWeight','bold', 'FontSize', 10, 'Interpreter','none');
        if i == 1 && dispOpt.force_same_rawrow && ~isempty(raw_ylim)
            ylim(raw_ylim);
        end
        if i < nRows, xticklabels([]); else, xlabel('Time [\mus]'); end

        nexttile;
        img = stored_cwt{i};
        local_max = max(img(:));
        is_S0_like = contains(lower(string(plotList{i}.name)), 's0');
        is_weak = (local_max < 0.2*all_max);

        if (is_S0_like || is_weak) && local_max > 0 && i > 1
            norm_base = local_max;
        else
            norm_base = all_max;
        end
        if i == 1 && dispOpt.force_same_rawrow
            norm_base = ref_max_raw;
        end

        img_show = img ./ (norm_base + eps);
        img_show(img_show < dispOpt.threshold) = 0;
        img_show(img_show > 1) = 1;
        img_show = img_show .^ dispOpt.gamma;

        imagesc(t_us, f_kHz_axis, img_show);
        axis xy; colormap(jet(256)); caxis([0 1]);
        ylim([0, dispOpt.fmax_kHz]); xlim([0, t_us(end)]);
        ylabel('Freq [kHz]');
        if i < nRows, xticklabels([]); else, xlabel('Time [\mus]'); end

        if curveOpt.enable && isfield(curve_data,'ok') && curve_data.ok
            hold on;
            plot(curve_data.A0.t_us, curve_data.A0.f_kHz, curveOpt.A0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
            plot(curve_data.S0.t_us, curve_data.S0.f_kHz, curveOpt.S0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
            if i == 1, add_curve_legend(gca, curveOpt); end
            hold off;
        end
    end

    cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'Norm |CWT|';
    title(tl, figTitle, 'FontSize', 12, 'FontWeight','bold', 'Interpreter','none');
end

function add_curve_legend(ax, curveOpt)
    axes(ax); %#ok<LAXES>
    xl = xlim(ax); yl = ylim(ax);
    x0 = xl(2) - 0.22*(xl(2)-xl(1));
    y0 = yl(2) - 0.08*(yl(2)-yl(1));
    hold on;
    plot([x0, x0+0.05*(xl(2)-xl(1))], [y0, y0], curveOpt.A0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
    text(x0+0.06*(xl(2)-xl(1)), y0, 'A0', 'Color', curveOpt.color, 'FontWeight','bold', 'VerticalAlignment','middle');
    y1 = y0 - 0.05*(yl(2)-yl(1));
    plot([x0, x0+0.05*(xl(2)-xl(1))], [y1, y1], curveOpt.S0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
    text(x0+0.06*(xl(2)-xl(1)), y1, 'S0', 'Color', curveOpt.color, 'FontWeight','bold', 'VerticalAlignment','middle');
    hold off;
end
