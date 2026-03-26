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
curveOpt.shift_us = 66;             % default/global shift
curveOpt.use_manual_shift_map = true;
curveOpt.manual_shift_dist_mm = [100 150 200];
curveOpt.manual_shift_us_map  = [66 40 20];
curveOpt.manual_shift_tol_mm  = 5;
curveOpt.allow_default_shift  = true;
curveOpt.lineWidth = 1.6;
curveOpt.A0_style = '-';
curveOpt.S0_style = '--';
curveOpt.A0_color = [1 1 1];        % white
curveOpt.S0_color = [1 0.85 0.15];  % yellow/orange for clear separation
curveOpt.legend_location = 'northwest';

% -------- waveReader headerlength --------
ioOpt.try_headerlength = true;
ioOpt.headerlength = 502;

% -------- ExtractFeati + unified modal identification --------
outOpt.save_diag_fig = true;
outOpt.save_feature_rank = true;

srcOpt.enable = true;
srcOpt.parse_from_filename = true;
srcOpt.unknown_label = 'unknown';
srcOpt.use_default_when_unknown = true;   % current surface-only stage
srcOpt.default_source_label = 'surface';
srcOpt.surface_keywords = {'surface','surf','top','face','flat','center','centre','su'};
srcOpt.edge_keywords    = {'edge','side','lateral','bord','sd'};
srcOpt.feature_names = {'logRE','logRP','fc_ratio','fp_ratio','bw_ratio','dur_ratio','dt_ratio','rhoRecon','split_conf'};
srcOpt.min_ref_per_class = 3;
srcOpt.do_leave_one_distance_out = true;
srcOpt.min_finite_frac = 0.80;
srcOpt.max_eta2 = 0.35;
srcOpt.max_selected_features = 12;

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
rowStruct = struct([]);

for fi = 1:nFiles
    wavefile = fullfile(pname, fname{fi});
    fprintf('\n==================== [%d/%d] %s ====================\n', fi, nFiles, fname{fi});

    [wave, fs] = readWaveRobust_plus(wavefile, ioOpt);
    wave = double(wave); fs = double(fs);

    [chSel, evList, sigList] = pickMultiEvents_from_wave(wave);
    nEv = numel(evList);

    dist_mm = parse_distance_mm(fname{fi});
    if isnan(dist_mm), dist_mm = curveOpt.dist_mm_default; end
    srcLabel = parse_source_label(fname{fi}, srcOpt);
    shift_us_this = get_shift_us_for_file(fname{fi}, dist_mm, curveOpt);

    curve_data = struct('ok',false,'msg','','A0',struct(),'S0',struct());
    if curveOpt.enable && ~isempty(curve_raw)
        [ok, A0, S0, msg] = parse_vallen_like_cell(curve_raw, dist_mm, shift_us_this);
        if ok
            curve_data.ok = true; curve_data.A0 = A0; curve_data.S0 = S0; curve_data.msg = msg;
        else
            curve_data = load_dispersion_curve(curve_file, dist_mm, shift_us_this);
        end
        if ~curve_data.ok
            warning('曲线解析失败（dist=%gmm，shift=%.1fus）：%s；本文件将不叠加曲线/不 split IMF1。', dist_mm, shift_us_this, curve_data.msg);
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
        one.shift_us = shift_us_this;

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

        featRaw = extractfeati_struct(t_us, raw_sig,       fs, chSel, evSel);
        featS0  = extractfeati_struct(t_us, one.S0_recon,  fs, chSel, evSel);
        featA0  = extractfeati_struct(t_us, one.A0_recon,  fs, chSel, evSel);
        modalFeat = build_modal_extract_features(featS0, featA0);

        rowS = build_event_row_struct(fname{fi}, dist_mm, srcLabel, fs, chSel, evSel, N, one, ...
                                      featRaw, featS0, featA0, modalFeat);
        if isempty(rowStruct)
            rowStruct = rowS;
        else
            rowStruct(end+1,1) = rowS; %#ok<AGROW>
        end
    end
end

if isempty(rowStruct)
    T = table();
else
    T = struct2table(rowStruct, 'AsArray', true);
end

model = struct();
cvres = struct();
if srcOpt.enable && ~isempty(T)
    [T, model, cvres] = build_unified_modal_model(T, srcOpt);

    if outOpt.save_feature_rank && isfield(model,'feature_stats') && istable(model.feature_stats)
        try
            writetable(model.feature_stats, fullfile(outOpt.save_dir, 'Modal_feature_ranking.csv'));
        catch ME
            warning('保存特征排序失败：%s', ME.message);
        end
    end
    if isfield(model,'selected_features') && ~isempty(model.selected_features)
        try
            writecell(model.selected_features(:), fullfile(outOpt.save_dir, 'Selected_modal_features.csv'));
        catch ME
            warning('保存已选特征失败：%s', ME.message);
        end
    end

    if outOpt.save_diag_fig
        try
            plot_unified_diagnostics(T, outOpt.save_dir);
        catch ME
            warning('诊断图生成失败：%s', ME.message);
        end
    end
end

writetable(T, fullfile(outOpt.save_dir, 'results.csv'));
save(fullfile(outOpt.save_dir, 'results.mat'), 'T', 'model', 'cvres', ...
    'outOpt','dataOpt','emdOpt','dispOpt','sepOpt','curveOpt','srcOpt');

fprintf('\n==================== DONE ====================\n');
fprintf('Output: %s\n', outOpt.save_dir);
if isfield(cvres,'overall_accuracy') && isfinite(cvres.overall_accuracy)
    fprintf('LODO accuracy = %.3f\n', cvres.overall_accuracy);
end
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
    one.raw_sig = raw_sig;
    one.S0_recon = S0_recon;
    one.A0_recon = A0_recon;
    one.imfs = imfs;

    % -------- unified modal features for source identification --------
    [fcS_kHz, fpS_kHz, bwS_kHz] = signal_fft_features(S0_recon, fs);
    [fcA_kHz, fpA_kHz, bwA_kHz] = signal_fft_features(A0_recon, fs);

    ES = sum(S0_recon(:).^2);
    EA = sum(A0_recon(:).^2);
    PS = max(abs(S0_recon(:)));
    PA = max(abs(A0_recon(:)));

    [tS_pk_us, durS_us] = envelope_peak_and_duration(S0_recon, fs, 0.10);
    [tA_pk_us, durA_us] = envelope_peak_and_duration(A0_recon, fs, 0.10);
    dt_pk_us = tA_pk_us - tS_pk_us;

    dt_pred_us = predict_dt_from_curve(curve_data, fcS_kHz, fcA_kHz);
    if isfinite(dt_pred_us) && abs(dt_pred_us) > eps
        dt_ratio = dt_pk_us ./ dt_pred_us;
    else
        dt_ratio = NaN;
    end

    ET = ES + EA + eps;
    PT = PS + PA + eps;
    rhoRecon = ET / (sum(raw_sig(:).^2) + eps);
    rhoResidual = max(0, 1 - rhoRecon);

    rhoE = ES / ET;
    rhoP = PS / PT;
    logRE = log10((ES + eps) / (EA + eps));
    logRP = log10((PS + eps) / (PA + eps));
    fc_ratio  = safe_sym_ratio(fcS_kHz, fcA_kHz);
    fp_ratio  = safe_sym_ratio(fpS_kHz, fpA_kHz);
    bw_ratio  = safe_sym_ratio(bwS_kHz, bwA_kHz);
    dur_ratio = safe_sym_ratio(durS_us, durA_us);

    imf1_support = 0;
    if isfinite(imf1_A0_energy_frac), imf1_support = imf1_support + imf1_A0_energy_frac; end
    if isfinite(imf1_S0_energy_frac), imf1_support = imf1_support + imf1_S0_energy_frac; end
    split_conf = min(1, max(0, 0.55*double(imf1_split_ok) + 0.25*min(rhoRecon,1) + 0.20*min(imf1_support,1)));

    one.ES = ES;
    one.EA = EA;
    one.rhoE = rhoE;
    one.logRE = logRE;

    one.PS = PS;
    one.PA = PA;
    one.rhoP = rhoP;
    one.logRP = logRP;

    one.fcS_kHz = fcS_kHz;
    one.fcA_kHz = fcA_kHz;
    one.fc_ratio = fc_ratio;

    one.fpS_kHz = fpS_kHz;
    one.fpA_kHz = fpA_kHz;
    one.fp_ratio = fp_ratio;

    one.bwS_kHz = bwS_kHz;
    one.bwA_kHz = bwA_kHz;
    one.bw_ratio = bw_ratio;

    one.durS_us = durS_us;
    one.durA_us = durA_us;
    one.dur_ratio = dur_ratio;

    one.tS_pk_us = tS_pk_us;
    one.tA_pk_us = tA_pk_us;
    one.dt_pk_us = dt_pk_us;
    one.dt_pred_us = dt_pred_us;
    one.dt_ratio = dt_ratio;

    one.rhoRecon = rhoRecon;
    one.rhoResidual = rhoResidual;
    one.split_conf = split_conf;
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
    toks = regexp(char(fname), '(\d+)\s*mm', 'tokens');
    if isempty(toks), return; end
    vals = nan(numel(toks),1);
    for ii = 1:numel(toks)
        vals(ii) = str2double(toks{ii}{1});
    end
    vals = vals(isfinite(vals));
    if isempty(vals), return; end
    % Prefer the largest mm value so filenames like "3mmPMMA_100mm_edge-1.wave"
    % are interpreted as distance=100 mm rather than thickness=3 mm.
    dist_mm = max(vals);
end

function s = strip_ext(fname)
    [~,s,~] = fileparts(fname);
end

function label = parse_source_label(fname, srcOpt)
    label = srcOpt.unknown_label;
    if ~srcOpt.parse_from_filename
        if isfield(srcOpt,'use_default_when_unknown') && srcOpt.use_default_when_unknown
            label = srcOpt.default_source_label;
        end
        return;
    end
    s = lower(char(fname));
    if isfield(srcOpt,'edge_keywords') && ~isempty(srcOpt.edge_keywords)
        for ii = 1:numel(srcOpt.edge_keywords)
            if contains(s, lower(char(srcOpt.edge_keywords{ii})))
                label = 'edge';
                return;
            end
        end
    end
    % If the filename does not contain any explicit edge keyword, treat it as surface.
    label = 'surface';
end

function shift_us = get_shift_us_for_file(fname, dist_mm, curveOpt)
    shift_us = curveOpt.shift_us;
    if isfield(curveOpt,'use_manual_shift_map') && curveOpt.use_manual_shift_map
        dlist = curveOpt.manual_shift_dist_mm(:);
        slist = curveOpt.manual_shift_us_map(:);
        tolmm = curveOpt.manual_shift_tol_mm;
        if isempty(tolmm) || ~isfinite(tolmm), tolmm = 5; end

        idx = find(abs(dlist - dist_mm) <= tolmm, 1, 'first');
        if isempty(idx)
            d2 = parse_distance_mm(fname);
            if isfinite(d2)
                idx = find(abs(dlist - d2) <= tolmm, 1, 'first');
            end
        end

        if ~isempty(idx)
            shift_us = slist(idx);
            return;
        end

        if ~(isfield(curveOpt,'allow_default_shift') && curveOpt.allow_default_shift)
            error('未找到距离 %.1f mm 对应的手动 shift 设置。', dist_mm);
        end
    end
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

function [fc_kHz, fp_kHz, bw_kHz] = signal_fft_features(x, fs)
    x = x(:);
    x(~isfinite(x)) = 0;
    x = x - mean(x);
    if all(abs(x) < eps)
        fc_kHz = NaN; fp_kHz = NaN; bw_kHz = NaN;
        return;
    end
    nfft = 2^nextpow2(max(numel(x), 1024));
    X = abs(fft(x, nfft)).^2;
    f = (0:nfft-1)/nfft*fs;
    half = 2:floor(nfft/2);
    P = X(half);
    ff = f(half);
    if isempty(P) || sum(P) <= 0
        fc_kHz = NaN; fp_kHz = NaN; bw_kHz = NaN;
        return;
    end
    [~, imax] = max(P);
    fp_Hz = ff(imax);
    fc_Hz = sum(ff(:).*P(:)) / (sum(P(:)) + eps);
    bw_Hz = sqrt(sum(((ff(:)-fc_Hz).^2).*P(:)) / (sum(P(:)) + eps));
    fc_kHz = fc_Hz / 1e3;
    fp_kHz = fp_Hz / 1e3;
    bw_kHz = bw_Hz / 1e3;
end

function [tpk_us, dur_us] = envelope_peak_and_duration(x, fs, thr_frac)
    x = x(:);
    if nargin < 3 || isempty(thr_frac), thr_frac = 0.10; end
    if all(abs(x) < eps)
        tpk_us = NaN; dur_us = NaN;
        return;
    end
    env = abs(hilbert(x));
    [~, ipk] = max(env);
    tpk_us = (ipk-1)/fs*1e6;
    thr = thr_frac * max(env);
    idx = find(env >= thr);
    if isempty(idx)
        dur_us = NaN;
    else
        dur_us = (idx(end) - idx(1))/fs*1e6;
    end
end

function r = safe_sym_ratio(a, b)
    if ~isfinite(a) || ~isfinite(b)
        r = NaN;
        return;
    end
    r = (a - b) / (a + b + eps);
end

function dt_pred_us = predict_dt_from_curve(curve_data, fcS_kHz, fcA_kHz)
    dt_pred_us = NaN;
    if ~isfield(curve_data,'ok') || ~curve_data.ok
        return;
    end
    if ~isfinite(fcS_kHz) || ~isfinite(fcA_kHz)
        return;
    end
    try
        tS = interp1(curve_data.S0.f_kHz(:), curve_data.S0.t_us(:), fcS_kHz, 'linear', 'extrap');
        tA = interp1(curve_data.A0.f_kHz(:), curve_data.A0.t_us(:), fcA_kHz, 'linear', 'extrap');
        if isfinite(tS) && isfinite(tA)
            dt_pred_us = tA - tS;
        end
    catch
        dt_pred_us = NaN;
    end
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

    % Use the raw-signal CWT maximum as the shared normalization reference
    % so every TF subplot can be compared directly in one common color scale.
    global_norm_base = ref_max_raw;
    if isempty(global_norm_base) || ~isfinite(global_norm_base) || global_norm_base <= 0
        global_norm_base = all_max;
    end
    if global_norm_base <= 0, global_norm_base = 1; end

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

        ax_tf = nexttile;
        img = stored_cwt{i};
        img_show = img ./ (global_norm_base + eps);
        img_show(img_show < dispOpt.threshold) = 0;
        img_show(img_show > 1) = 1;
        img_show = img_show .^ dispOpt.gamma;

        imagesc(t_us, f_kHz_axis, img_show, 'Parent', ax_tf);
        axis(ax_tf, 'xy'); colormap(ax_tf, jet(256)); caxis(ax_tf, [0 1]);
        ylim(ax_tf, [0, dispOpt.fmax_kHz]); xlim(ax_tf, [0, t_us(end)]);
        ylabel(ax_tf, 'Freq [kHz]');
        if i < nRows, xticklabels(ax_tf, []); else, xlabel(ax_tf, 'Time [\mus]'); end

        if curveOpt.enable && isfield(curve_data,'ok') && curve_data.ok
            hold(ax_tf, 'on');
            hA = plot(ax_tf, curve_data.A0.t_us, curve_data.A0.f_kHz, curveOpt.A0_style, ...
                'Color', curveOpt.A0_color, 'LineWidth', curveOpt.lineWidth);
            hS = plot(ax_tf, curve_data.S0.t_us, curve_data.S0.f_kHz, curveOpt.S0_style, ...
                'Color', curveOpt.S0_color, 'LineWidth', curveOpt.lineWidth);
            if i == 1
                add_curve_legend(ax_tf, hA, hS, curveOpt);
            end
            hold(ax_tf, 'off');
        end
    end

    cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'Norm |CWT| (shared vs raw)';
    title(tl, figTitle, 'FontSize', 12, 'FontWeight','bold', 'Interpreter','none');
end

function add_curve_legend(ax, hA, hS, curveOpt)
    lgd = legend(ax, [hA hS], {'A0','S0'}, 'Location', curveOpt.legend_location, ...
        'Box', 'on', 'AutoUpdate', 'off');
    try
        lgd.TextColor = 'k';
        lgd.Color = [1 1 1];
    catch
    end
end


%% =====================================================================
%% ======================= UNIFIED SOURCE IDENTIFICATION ================

function [T, model, cvres] = build_unified_modal_model(T, srcOpt)
    model = struct();
    cvres = struct('overall_accuracy', NaN, 'n_eval', 0, 'note', '');

    if isempty(T)
        return;
    end

    vnames = T.Properties.VariableNames;
    baseFeat = srcOpt.feature_names(ismember(srcOpt.feature_names, vnames));
    modalFeat = vnames(startsWith(vnames, 'Mlr_') | startsWith(vnames, 'Msd_'));
    cand = unique([baseFeat, modalFeat], 'stable');

    if isempty(cand)
        warning('未找到可用于统一模态识别的特征列。');
        return;
    end

    statsRows = {};
    D = double(T.distance_mm);
    for j = 1:numel(cand)
        x = to_numeric_column(T.(cand{j}));
        finite_frac = mean(isfinite(x));
        eta2 = eta2_by_distance(x, D);
        spread = robust_spread(x);
        statsRows(end+1,:) = {cand{j}, finite_frac, eta2, spread}; %#ok<AGROW>
    end
    featStats = cell2table(statsRows, 'VariableNames', {'feature','finite_frac','eta2_distance','robust_spread'});
    featStats = sortrows(featStats, {'eta2_distance','finite_frac','robust_spread'}, {'ascend','descend','descend'});
    model.feature_stats = featStats;
    model.candidate_features = cand;

    keep = featStats.finite_frac >= srcOpt.min_finite_frac & isfinite(featStats.eta2_distance);
    if isfield(srcOpt,'max_eta2') && isfinite(srcOpt.max_eta2)
        keep = keep & (featStats.eta2_distance <= srcOpt.max_eta2);
    end
    stable = featStats(keep,:);
    if isempty(stable)
        stable = featStats(1:min(height(featStats), max(1, numel(baseFeat))), :);
    end
    nSel = min(height(stable), srcOpt.max_selected_features);
    selected = stable.feature(1:nSel);
    model.selected_features = selected;

    X = table_to_numeric(T, selected);
    [Xz, scaler] = robust_standardize(X);
    model.scaler = scaler;

    [score, coeff, latent, mu] = local_pca(Xz, min(3, size(Xz,2)));
    model.pca_coeff = coeff;
    model.pca_latent = latent;
    model.pca_mu = mu;

    if size(score,2) >= 1, T.modal_PC1 = score(:,1); else, T.modal_PC1 = nan(height(T),1); end
    if size(score,2) >= 2, T.modal_PC2 = score(:,2); else, T.modal_PC2 = nan(height(T),1); end
    if size(score,2) >= 3, T.modal_PC3 = score(:,3); else, T.modal_PC3 = nan(height(T),1); end

    [T, srcModel, cvres] = maybe_build_source_classifier(T, selected, srcOpt);
    model.source_model = srcModel;
end

function [T, model, cvres] = maybe_build_source_classifier(T, featNames, srcOpt)
    model = struct();
    cvres = struct('overall_accuracy', NaN, 'n_eval', 0, 'note', '');

    if ~ismember('source_label', T.Properties.VariableNames)
        cvres.note = 'source_label 缺失，跳过源分类。';
        return;
    end

    labels = lower(string(T.source_label));
    known = labels ~= "" & labels ~= string(srcOpt.unknown_label) & labels ~= "unknown";
    ucls = unique(labels(known));
    if numel(ucls) < 2
        cvres.note = '当前仅有一个 source_label 类别；已完成距离鲁棒特征筛选，但跳过源分类器。';
        return;
    end

    Xall = table_to_numeric(T, featNames);
    [Xz, scaler] = robust_standardize(Xall);
    model.feature_names = featNames;
    model.scaler = scaler;
    classes = cellstr(ucls(:)).';
    model.classes = classes;

    dmat = nan(height(T), numel(classes));
    mu = cell(1, numel(classes));
    invS = cell(1, numel(classes));
    nClass = zeros(1, numel(classes));

    for c = 1:numel(classes)
        idx = known & strcmp(labels, classes{c});
        nClass(c) = nnz(idx);
        [mu{c}, invS{c}] = fit_mahal_model(Xz(idx,:));
        dmat(:,c) = mahal_predict(Xz, mu{c}, invS{c});
    end

    model.mu = mu;
    model.invS = invS;
    model.nClass = nClass;

    pred = repmat(string(srcOpt.unknown_label), height(T), 1);
    pred_margin = nan(height(T),1);
    pred_conf = nan(height(T),1);

    for i = 1:height(T)
        di = dmat(i,:);
        if ~all(isfinite(di)), continue; end
        [ds, ord] = sort(di, 'ascend');
        pred(i) = string(classes{ord(1)});
        if numel(ds) >= 2
            pred_margin(i) = ds(2) - ds(1);
            pred_conf(i) = 1 ./ (1 + exp(-(ds(2) - ds(1))));
        end
    end

    for c = 1:numel(classes)
        T.(['score_', classes{c}]) = dmat(:,c);
    end
    T.pred_source_label = cellstr(pred);
    T.pred_margin = pred_margin;
    T.pred_conf = pred_conf;

    if srcOpt.do_leave_one_distance_out
        cvres = leave_one_distance_out_cv_generic(T, featNames, srcOpt, classes);
    end
end

function val = robust_spread(x)
    x = x(isfinite(x));
    if isempty(x)
        val = NaN;
        return;
    end
    medx = median(x);
    val = median(abs(x - medx));
    if ~isfinite(val) || val < eps
        val = std(x);
    end
end

function eta2 = eta2_by_distance(x, d)
    x = x(:); d = d(:);
    good = isfinite(x) & isfinite(d);
    x = x(good); d = d(good);
    if numel(x) < 3 || numel(unique(d)) < 2
        eta2 = NaN;
        return;
    end
    mu = mean(x);
    SSb = 0;
    SSt = sum((x - mu).^2);
    ud = unique(d);
    for i = 1:numel(ud)
        idx = d == ud(i);
        if ~any(idx), continue; end
        ni = nnz(idx);
        mi = mean(x(idx));
        SSb = SSb + ni * (mi - mu).^2;
    end
    eta2 = SSb / max(SSt, eps);
end

function [score, coeff, latent, mu] = local_pca(X, nComp)
    if isempty(X)
        score = []; coeff = []; latent = []; mu = [];
        return;
    end
    X = double(X);
    mu = mean(X,1,'omitnan');
    X(~isfinite(X)) = 0;
    Xc = X - mu;
    [~, S, V] = svd(Xc, 'econ');
    coeff = V(:,1:min(nComp, size(V,2)));
    score = Xc * coeff;
    s = diag(S);
    latent = s.^2 / max(size(X,1)-1,1);
    latent = latent(1:min(nComp, numel(latent)));
end

function cvres = leave_one_distance_out_cv_generic(T, featNames, srcOpt, classes)
    cvres = struct('overall_accuracy', NaN, 'n_eval', 0, 'by_distance', table(), 'note', '');
    labels = lower(string(T.source_label));
    known = labels ~= "" & labels ~= string(srcOpt.unknown_label) & labels ~= "unknown";
    D = T.distance_mm;
    uD = unique(D(known));
    pred_all = strings(0,1);
    true_all = strings(0,1);
    dist_log = [];
    rows = {};

    for ii = 1:numel(uD)
        d0 = uD(ii);
        test = known & (D == d0);
        train = known & (D ~= d0);
        if nnz(test) < 1 || nnz(train) < numel(classes)
            continue;
        end

        Xtr = table_to_numeric(T(train,:), featNames);
        Xte = table_to_numeric(T(test,:),  featNames);

        [Xtrz, scaler] = robust_standardize(Xtr);
        Xtez = apply_scaler(Xte, scaler);

        ytr = lower(string(T.source_label(train)));
        yte = lower(string(T.source_label(test)));

        dmat = nan(size(Xtez,1), numel(classes));
        for c = 1:numel(classes)
            idxc = strcmp(ytr, classes{c});
            if nnz(idxc) < 1, continue; end
            [mu, invS] = fit_mahal_model(Xtrz(idxc,:));
            dmat(:,c) = mahal_predict(Xtez, mu, invS);
        end

        pred = repmat("unknown", size(Xtez,1), 1);
        for k = 1:size(Xtez,1)
            if all(isfinite(dmat(k,:)))
                [~, ic] = min(dmat(k,:));
                pred(k) = classes{ic};
            end
        end

        acc = mean(pred == yte);
        rows(end+1,:) = {d0, nnz(test), acc}; %#ok<AGROW>
        pred_all = [pred_all; pred(:)]; %#ok<AGROW>
        true_all = [true_all; yte(:)]; %#ok<AGROW>
        dist_log = [dist_log; repmat(d0, numel(yte), 1)]; %#ok<AGROW>
    end

    if ~isempty(rows)
        cvres.by_distance = cell2table(rows, 'VariableNames', {'distance_mm','n_test','accuracy'});
    end
    good = (pred_all ~= "unknown") & (true_all ~= "unknown");
    if any(good)
        cvres.overall_accuracy = mean(pred_all(good) == true_all(good));
        cvres.n_eval = nnz(good);
        cvres.dist_log = dist_log(good);
        cvres.pred = cellstr(pred_all(good));
        cvres.truth = cellstr(true_all(good));
    end
end

function [Xz, scaler] = robust_standardize(X)
    Xz = X;
    p = size(X,2);
    scaler.center = nan(1,p);
    scaler.scale  = nan(1,p);
    for j = 1:p
        col = X(:,j);
        medj = median(col(isfinite(col)));
        if isempty(medj), medj = 0; end
        madj = median(abs(col(isfinite(col)) - medj));
        if isempty(madj) || madj < eps
            sdj = std(col(isfinite(col)));
            if isempty(sdj) || sdj < eps, sdj = 1; end
            madj = sdj;
        end
        scaler.center(j) = medj;
        scaler.scale(j) = madj;
        good = isfinite(col);
        Xz(good,j) = (col(good) - medj) / (madj + eps);
        Xz(~good,j) = 0;
    end
end

function X = table_to_numeric(T, featNames)
    X = nan(height(T), numel(featNames));
    for j = 1:numel(featNames)
        x = T.(featNames{j});
        X(:,j) = to_numeric_column(x);
    end
end

function x = to_numeric_column(x)
    if iscell(x)
        try
            x = cellfun(@double, x);
        catch
            x = nan(numel(x),1);
        end
    end
    x = double(x);
end

function [mu, invS] = fit_mahal_model(X)
    if isempty(X)
        mu = [];
        invS = [];
        return;
    end
    X = double(X);
    mu = mean(X,1,'omitnan');
    X(~isfinite(X)) = 0;
    Xc = X - mu;
    if size(X,1) <= 1
        S = eye(size(X,2));
    else
        S = (Xc' * Xc) / max(size(X,1)-1,1);
    end
    if ~all(isfinite(S(:))) || isempty(S)
        S = eye(size(X,2));
    end
    lam = 0.15 * trace(S) / max(size(S,1),1);
    if ~isfinite(lam) || lam <= 0, lam = 1e-3; end
    S = S + lam * eye(size(S));
    invS = pinv(S);
end

function d = mahal_predict(X, mu, invS)
    if isempty(mu) || isempty(invS)
        d = nan(size(X,1),1);
        return;
    end
    X = double(X);
    X(~isfinite(X)) = 0;
    Xm = X - mu;
    d = sum((Xm * invS) .* Xm, 2);
end

function Xz = apply_scaler(X, scaler)
    Xz = X;
    for j = 1:size(X,2)
        col = X(:,j);
        good = isfinite(col);
        Xz(~good,j) = 0;
        Xz(good,j) = (col(good) - scaler.center(j)) / (scaler.scale(j) + eps);
    end
end

function plot_unified_diagnostics(T, outDir)
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    lab = lower(string(T.source_label));
    c = T.distance_mm;
    valid = isfinite(T.logRE) & isfinite(T.logRP);

    fig1 = figure('Color','w', 'Position', [60 60 900 680]);
    scatter(T.logRE(valid), T.logRP(valid), 70, c(valid), 'filled');
    hold on;
    idx1 = valid & strcmp(lab,'surface');
    idx2 = valid & strcmp(lab,'edge');
    plot(T.logRE(idx1), T.logRP(idx1), 'ko', 'MarkerSize', 9, 'LineWidth', 1.1);
    plot(T.logRE(idx2), T.logRP(idx2), 'k^', 'MarkerSize', 9, 'LineWidth', 1.1);
    grid on;
    xlabel('log_{10}(E_S / E_A)');
    ylabel('log_{10}(P_S / P_A)');
    title('Unified modal signature: logRE vs logRP');
    cb = colorbar; cb.Label.String = 'Distance [mm]';
    saveas(fig1, fullfile(outDir, 'Diag_logRE_logRP.png'));
    close(fig1);

    if ismember('modal_PC1', T.Properties.VariableNames) && ismember('modal_PC2', T.Properties.VariableNames)
        valid2 = isfinite(T.modal_PC1) & isfinite(T.modal_PC2);
        fig2 = figure('Color','w', 'Position', [80 80 900 680]);
        scatter(T.modal_PC1(valid2), T.modal_PC2(valid2), 70, T.distance_mm(valid2), 'filled');
        hold on;
        plot(T.modal_PC1(valid2 & strcmp(lab,'surface')), T.modal_PC2(valid2 & strcmp(lab,'surface')), ...
            'ko', 'MarkerSize', 9, 'LineWidth', 1.1);
        plot(T.modal_PC1(valid2 & strcmp(lab,'edge')), T.modal_PC2(valid2 & strcmp(lab,'edge')), ...
            'k^', 'MarkerSize', 9, 'LineWidth', 1.1);
        grid on;
        xlabel('Modal PC1');
        ylabel('Modal PC2');
        title('Distance-robust modal feature space');
        cb = colorbar; cb.Label.String = 'Distance [mm]';
        saveas(fig2, fullfile(outDir, 'Diag_modal_PC1_PC2.png'));
        close(fig2);
    end

    scoreVars = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'score_'));
    if numel(scoreVars) >= 2
        valid3 = isfinite(T.(scoreVars{1})) & isfinite(T.(scoreVars{2}));
        fig3 = figure('Color','w', 'Position', [100 100 900 680]);
        scatter(T.(scoreVars{1})(valid3), T.(scoreVars{2})(valid3), 70, T.distance_mm(valid3), 'filled');
        hold on;
        plot(xlim, xlim, 'k--', 'LineWidth', 1);
        grid on;
        xlabel(strrep(scoreVars{1}, '_', '\_'));
        ylabel(strrep(scoreVars{2}, '_', '\_'));
        title('Unified source-ID score space');
        cb = colorbar; cb.Label.String = 'Distance [mm]';
        saveas(fig3, fullfile(outDir, 'Diag_source_score_space.png'));
        close(fig3);
    end

    sumVars = {'file','distance_mm','source_label','shift_us','split_conf','logRE','logRP'};
    sumVars = sumVars(ismember(sumVars, T.Properties.VariableNames));
    if ismember('pred_source_label', T.Properties.VariableNames)
        sumVars = [sumVars, {'pred_source_label','pred_margin','pred_conf'}];
    end
    if ~isempty(sumVars)
        Tout = T(:, sumVars);
        writetable(Tout, fullfile(outDir, 'Predictions_summary.csv'));
    end
end

%% ===================== ExtractFeati integration ======================
function feat = extractfeati_struct(t_us, sig, fs, ch, ev)
    baseNames = extractfeati_feature_names();
    feat = nan_feature_struct(baseNames);
    if isempty(sig)
        return;
    end
    sig = double(sig(:));
    sig(~isfinite(sig)) = 0;
    T2 = double(t_us(:)) * 1e-6;
    if numel(T2) ~= numel(sig)
        T2 = (0:numel(sig)-1)'/fs;
    end
    heads_info = {double(fs), double(ch), double(ev), 0};
    persistent warned;
    if exist('ExtractFeati_OPT','file') ~= 2
        if isempty(warned)
            warning('未找到 ExtractFeati_OPT.m，相关特征将填 NaN。');
            warned = true;
        end
        return;
    end
    try
        F = ExtractFeati_OPT(T2, sig, [], heads_info);
        F = double(F(:));
        nUse = min(numel(F)-3, numel(baseNames));
        if nUse > 0
            for k = 1:nUse
                feat.(baseNames{k}) = F(k+3);
            end
        end
    catch ME
        if isempty(warned)
            warning('ExtractFeati_OPT 调用失败：%s', ME.message);
            warned = true;
        end
    end
end

function names = extractfeati_feature_names()
    names = {'A_dB','D_us','E','ZCR','RT_us','TC_us','alpha', ...
             'PP1','PP2','PP3','PP4', ...
             'FC2_kHz','PF2_kHz','SSpread_kHz','SSkew','SKurt','SSlope','SRoff_kHz', ...
             'sqrtSSpreadP','SSkewP','SKurtP','SRon_kHz', ...
             'WPE1','WPE2','WPE3','WPE4','WPE5','WPE6','WPE7','WPE8','Entropy'};
end

function feat = nan_feature_struct(names)
    feat = struct();
    for k = 1:numel(names)
        feat.(names{k}) = NaN;
    end
end

function tf = is_positive_extract_feature(name)
    signedNames = {'SSkew','SSkewP','SSlope'};
    tf = ~ismember(name, signedNames);
end

function out = build_modal_extract_features(featS0, featA0)
    out = struct();
    fns = fieldnames(featS0);
    for k = 1:numel(fns)
        nm = fns{k};
        s = featS0.(nm);
        a = featA0.(nm);
        if is_positive_extract_feature(nm)
            out.(sprintf('Mlr_%s', nm)) = safe_log_ratio(s, a);
        else
            out.(sprintf('Msd_%s', nm)) = safe_signed_diff(s, a);
        end
    end
end

function r = safe_log_ratio(a, b)
    if ~(isfinite(a) && isfinite(b))
        r = NaN;
        return;
    end
    aa = abs(a) + eps;
    bb = abs(b) + eps;
    r = log10(aa / bb);
end

function r = safe_signed_diff(a, b)
    if ~(isfinite(a) && isfinite(b))
        r = NaN;
        return;
    end
    r = (a - b) / (abs(a) + abs(b) + eps);
end

function out = prefix_struct_fields(s, prefix)
    out = struct();
    fns = fieldnames(s);
    for k = 1:numel(fns)
        out.([prefix, fns{k}]) = s.(fns{k});
    end
end

function row = build_event_row_struct(fileName, dist_mm, srcLabel, fs, chSel, evSel, N, one, featRaw, featS0, featA0, modalFeat)
    row = struct();
    row.file = fileName;
    row.distance_mm = dist_mm;
    row.source_label = srcLabel;
    row.shift_us = one.shift_us;
    row.fs_Hz = fs;
    row.ch = chSel;
    row.event = evSel;
    row.Nwin = N;
    row.K_imf = one.K;
    row.idx_A0_base = mat2str(one.idx_A0_base);
    row.idx_S0_base = mat2str(one.idx_S0_base);
    row.imf1_split_ok = one.imf1_split_ok;
    row.imf1_A0_energy_frac = one.imf1_A0_energy_frac;
    row.imf1_S0_energy_frac = one.imf1_S0_energy_frac;
    row.lowRatio_all = mat2str(one.lowRatio, 3);
    row.A0_lowRatio_mean = one.A0_lowRatio_mean;
    row.S0_lowRatio_mean = one.S0_lowRatio_mean;
    row.ES = one.ES;
    row.EA = one.EA;
    row.rhoE = one.rhoE;
    row.logRE = one.logRE;
    row.PS = one.PS;
    row.PA = one.PA;
    row.rhoP = one.rhoP;
    row.logRP = one.logRP;
    row.fcS_kHz = one.fcS_kHz;
    row.fcA_kHz = one.fcA_kHz;
    row.fc_ratio = one.fc_ratio;
    row.fpS_kHz = one.fpS_kHz;
    row.fpA_kHz = one.fpA_kHz;
    row.fp_ratio = one.fp_ratio;
    row.bwS_kHz = one.bwS_kHz;
    row.bwA_kHz = one.bwA_kHz;
    row.bw_ratio = one.bw_ratio;
    row.durS_us = one.durS_us;
    row.durA_us = one.durA_us;
    row.dur_ratio = one.dur_ratio;
    row.tS_pk_us = one.tS_pk_us;
    row.tA_pk_us = one.tA_pk_us;
    row.dt_pk_us = one.dt_pk_us;
    row.dt_pred_us = one.dt_pred_us;
    row.dt_ratio = one.dt_ratio;
    row.rhoRecon = one.rhoRecon;
    row.rhoResidual = one.rhoResidual;
    row.split_conf = one.split_conf;

    r1 = prefix_struct_fields(featRaw, 'raw_');
    r2 = prefix_struct_fields(featS0,  'S0_');
    r3 = prefix_struct_fields(featA0,  'A0_');
    rr = [fieldnames(r1); fieldnames(r2); fieldnames(r3); fieldnames(modalFeat)];
    % merge
    fns = fieldnames(r1);
    for k = 1:numel(fns), row.(fns{k}) = r1.(fns{k}); end
    fns = fieldnames(r2);
    for k = 1:numel(fns), row.(fns{k}) = r2.(fns{k}); end
    fns = fieldnames(r3);
    for k = 1:numel(fns), row.(fns{k}) = r3.(fns{k}); end
    fns = fieldnames(modalFeat);
    for k = 1:numel(fns), row.(fns{k}) = modalFeat.(fns{k}); end
end
