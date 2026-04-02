
%% ================== EMD_MAE_batch_IMF1split_full_v9_onV8.m ==================
% Batch EMD-MAE with dispersion-guided IMF1 splitting (exclusive masking)
% V9 (based on your V8 code, keeping the original file/event selection logic)
%
% Main updates vs V8:
%  1) Keep the original V8 file selection + per-file channel/event selection.
%  2) Keep per-event figure saving (mandatory for later checking).
%  3) Keep raw / S0 / A0 full feature export (no residual ExtractFeati export).
%  4) Replace summary diagnostic colorbar-style grouping with discrete legend groups
%     based on (source_label, distance_mm), e.g. surface 100 mm / surface 150 mm /
%     surface 200 mm / edge 100 mm.
%  5) Output additional grouped diagnostic figures for core feature pairs and for
%     all numeric features versus distance, all saved under the same output folder.
%
% Notes:
%  - The signal reading, event picking, EMD decomposition, IMF1 split, and per-event
%    reconstruction logic are intentionally preserved from V8 to avoid breaking a
%    workflow that already runs on your machine.
% ==============================================================================

clear; clc; close all;

%% ================= 0) USER OPTIONS =================
% -------- Output --------
outOpt.save_dir = fullfile(pwd, ['EMD_MAE_batch_IMF1split_out_', datestr(now,'yyyymmdd_HHMMSS')]);
outOpt.save_png = true;
outOpt.save_pdf = false;
outOpt.dpi      = 300;
outOpt.show_fig = false;          % false for batch
outOpt.save_diag_fig = true;
outOpt.save_event_mat = false;    % optional per-event MAT dump

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
curveOpt.use_manual_shift_map = true;
curveOpt.shift_map_dist_mm = [100 150 200];
curveOpt.shift_map_us      = [66 40 20];   % manual shift
curveOpt.shift_us_fallback = 66;
curveOpt.lineWidth = 1.6;
curveOpt.A0_style = '-';
curveOpt.S0_style = '--';
curveOpt.color = [1 1 1];           % white

% -------- waveReader headerlength --------
ioOpt.try_headerlength = true;
ioOpt.headerlength = 502;

% -------- Source label parsing --------
srcOpt.parse_from_filename = true;
srcOpt.edge_keywords    = {'edge','side','lateral','bord','sd'};
srcOpt.surface_default  = 'surface';

% -------- Basic feature exploration / plotting --------
basicOpt.enable = true;
basicOpt.topN = 8;
basicOpt.feature_names = {};   % auto-filled after results table is built

diagOpt.enable = true;
diagOpt.tile_per_fig = 6;
diagOpt.point_size = 70;
diagOpt.marker_edge_color = [0 0 0];
diagOpt.legend_location = 'eastoutside';
diagOpt.max_legend_columns = 1;
diagOpt.jitter_frac = 0.035;         % horizontal jitter fraction for x=distance plots
diagOpt.core_pairs = { ...
    {'logRE','logRP','Unified modal signature: logRE vs logRP','log_{10}(E_S / E_A)','log_{10}(P_S / P_A)'}, ...
    {'rhoE','rhoP','Unified modal signature: rhoE vs rhoP','\rho_E = E_S/(E_S+E_A)','\rho_P = P_S/(P_S+P_A)'}, ...
    {'fcS_kHz','fcA_kHz','Modal centroid frequency: f_c(S0) vs f_c(A0)','f_c(S0) [kHz]','f_c(A0) [kHz]'}, ...
    {'fpS_kHz','fpA_kHz','Modal peak frequency: f_p(S0) vs f_p(A0)','f_p(S0) [kHz]','f_p(A0) [kHz]'}, ...
    {'bwS_kHz','bwA_kHz','Modal bandwidth: BW(S0) vs BW(A0)','BW(S0) [kHz]','BW(A0) [kHz]'}, ...
    {'durS_us','durA_us','Envelope duration: Dur(S0) vs Dur(A0)','Dur(S0) [\mus]','Dur(A0) [\mus]'}, ...
    {'ES','EA','Modal energy: E_S vs E_A','E_S','E_A'}, ...
    {'PS','PA','Modal peak amplitude: P_S vs P_A','P_S','P_A'} ...
    };
diagOpt.modal_topN = 12;
diagOpt.feature_pages_only_valid = true;

%% ================= 1) Select wave files =================
fprintf('请选择数据文件（*.wave / *.tradb）... 可多选\n');
[fname, pname] = uigetfile({'*.wave;*.tradb;*.*'}, '选择 wave 文件（可多选）', 'MultiSelect','on');
if isequal(fname, 0), error('未选择 wave 文件'); end
if ischar(fname), fname = {fname}; end
nFiles = numel(fname);

if ~exist(outOpt.save_dir,'dir'); mkdir(outOpt.save_dir); end
figDir = fullfile(outOpt.save_dir, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end
diagDir = fullfile(outOpt.save_dir, 'diagnostics');
if ~exist(diagDir,'dir'); mkdir(diagDir); end
if outOpt.save_event_mat
    eventMatDir = fullfile(outOpt.save_dir, 'event_mat');
    if ~exist(eventMatDir,'dir'); mkdir(eventMatDir); end
else
    eventMatDir = '';
end

if outOpt.show_fig
    set(0,'DefaultFigureVisible','on');
else
    set(0,'DefaultFigureVisible','off');
end

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
extractNames = extractfeati_feature_names();
extractNamesNoMeta = setdiff(extractNames, {'hit','hittime_s','channel'}, 'stable');
head = [{'file','distance_mm','source_label','shift_us','fs_Hz','ch','event','Nwin','K_imf', ...
        'idx_A0_base','idx_S0_base','imf1_split_ok', ...
        'imf1_A0_energy_frac','imf1_S0_energy_frac', ...
        'lowRatio_all','A0_lowRatio_mean','S0_lowRatio_mean', ...
        'ES','EA','rhoE','logRE', ...
        'PS','PA','rhoP','logRP', ...
        'fcS_kHz','fcA_kHz','fc_ratio', ...
        'fpS_kHz','fpA_kHz','fp_ratio', ...
        'bwS_kHz','bwA_kHz','bw_ratio', ...
        'durS_us','durA_us','dur_ratio', ...
        'tS_pk_us','tA_pk_us','dt_pk_us', ...
        'dt_pred_us','dt_ratio', ...
        'rhoRecon','rhoResidual','split_conf'}, ...
        strcat('raw_', extractNames), strcat('S0_', extractNames), strcat('A0_', extractNames)];

% candidate features for ranking + plotting
basicOpt.feature_names = [{'logRE','logRP','rhoE','rhoP', ...
    'fc_ratio','fp_ratio','bw_ratio','dur_ratio','dt_ratio', ...
    'rhoRecon','split_conf', ...
    'fcS_kHz','fcA_kHz','fpS_kHz','fpA_kHz','bwS_kHz','bwA_kHz', ...
    'durS_us','durA_us','ES','EA','PS','PA'}, ...
    strcat('raw_',extractNamesNoMeta), strcat('S0_',extractNamesNoMeta), strcat('A0_',extractNamesNoMeta)];

for fi = 1:nFiles
    wavefile = fullfile(pname, fname{fi});
    fprintf('\n==================== [%d/%d] %s ====================\n', fi, nFiles, fname{fi});

    [wave, fs] = readWaveRobust_plus(wavefile, ioOpt);
    wave = double(wave); fs = double(fs);

    % ===== Keep V8 selection logic unchanged =====
    [chSel, evList, sigList] = pickMultiEvents_from_wave(wave);
    nEv = numel(evList);

    dist_mm = parse_distance_mm(fname{fi});
    if isnan(dist_mm), dist_mm = curveOpt.dist_mm_default; end
    srcLabel = parse_source_label(fname{fi}, srcOpt);
    shift_us = manual_shift_us_from_distance(dist_mm, curveOpt);

    curve_data = struct('ok',false,'msg','','A0',struct(),'S0',struct());
    if curveOpt.enable && ~isempty(curve_raw)
        [ok, A0, S0, msg] = parse_vallen_like_cell(curve_raw, dist_mm, shift_us);
        if ok
            curve_data.ok = true; curve_data.A0 = A0; curve_data.S0 = S0; curve_data.msg = msg;
        else
            curve_data = load_dispersion_curve(curve_file, dist_mm, shift_us);
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

        featRaw = extractfeati_struct(t_us, raw_sig,       fs, chSel, evSel);
        featS0  = extractfeati_struct(t_us, one.S0_recon,  fs, chSel, evSel);
        featA0  = extractfeati_struct(t_us, one.A0_recon,  fs, chSel, evSel);

        base = sprintf('%s__d%gmm__ch%d__ev%04d', strip_ext(fname{fi}), dist_mm, chSel, evSel);

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

        if outOpt.save_event_mat
            save(fullfile(eventMatDir, [base,'.mat']), ...
                'raw_sig','t_us','fs','one','featRaw','featS0','featA0','curve_data','curveOpt','emdOpt','dispOpt','sepOpt');
        end

        if ~outOpt.show_fig, close(fig1); close(fig2); end

        rows(end+1,:) = [{fname{fi}, dist_mm, srcLabel, shift_us, fs, chSel, evSel, N, ...
            one.K, mat2str(one.idx_A0_base), mat2str(one.idx_S0_base), one.imf1_split_ok, ...
            one.imf1_A0_energy_frac, one.imf1_S0_energy_frac, ...
            mat2str(one.lowRatio,3), one.A0_lowRatio_mean, one.S0_lowRatio_mean, ...
            one.ES, one.EA, one.rhoE, one.logRE, ...
            one.PS, one.PA, one.rhoP, one.logRP, ...
            one.fcS_kHz, one.fcA_kHz, one.fc_ratio, ...
            one.fpS_kHz, one.fpA_kHz, one.fp_ratio, ...
            one.bwS_kHz, one.bwA_kHz, one.bw_ratio, ...
            one.durS_us, one.durA_us, one.dur_ratio, ...
            one.tS_pk_us, one.tA_pk_us, one.dt_pk_us, ...
            one.dt_pred_us, one.dt_ratio, ...
            one.rhoRecon, one.rhoResidual, one.split_conf}, ...
            feature_struct_to_cell(featRaw, extractNames), ...
            feature_struct_to_cell(featS0,  extractNames), ...
            feature_struct_to_cell(featA0,  extractNames)]; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', head);

rankTab = table();
selectedTab = table();
modalRankTab = table();
modalSelectedTab = table();
if basicOpt.enable && ~isempty(T)
    try
        [rankTab, selectedTab] = rank_basic_features(T, basicOpt);
        writetable(rankTab, fullfile(outOpt.save_dir, 'Feature_ranking_basic.csv'));
        writetable(selectedTab, fullfile(outOpt.save_dir, 'Selected_basic_features.csv'));

        [modalRankTab, modalSelectedTab] = rank_modal_features(T, extractNamesNoMeta, diagOpt.modal_topN);
        writetable(modalRankTab, fullfile(outOpt.save_dir, 'Feature_ranking_modal_S0_vs_A0.csv'));
        writetable(modalSelectedTab, fullfile(outOpt.save_dir, 'Selected_modal_features_S0_vs_A0.csv'));

        if outOpt.save_diag_fig && diagOpt.enable
            plot_basic_diagnostics_v9(T, rankTab, selectedTab, outOpt.save_dir, diagDir, diagOpt, basicOpt);
            plot_mode_specific_feature_pages(T, diagDir, diagOpt, extractNamesNoMeta, 'S0');
            plot_mode_specific_feature_pages(T, diagDir, diagOpt, extractNamesNoMeta, 'A0');
            plot_modal_selected_feature_pages(T, diagDir, diagOpt, modalSelectedTab);
        end
    catch ME
        warning('基础特征诊断图生成失败：%s', ME.message);
    end
end

writetable(T, fullfile(outOpt.save_dir, 'results.csv'));
save(fullfile(outOpt.save_dir, 'results.mat'), 'T', 'rankTab', 'selectedTab', 'modalRankTab', 'modalSelectedTab', 'extractNames', ...
    'outOpt','dataOpt','emdOpt','dispOpt','sepOpt','curveOpt','srcOpt','basicOpt','diagOpt');

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
            list_decomp{end}.name  = sprintf('IMF1 -> S0-part (mask, \\sigma=%.1f\\mus, \\delta=%.1f\\mus)', sepOpt.imf1_sigma_us, sepOpt.imf1_delta_us);
            list_decomp{end+1}.sig = imf1_A0;
            list_decomp{end}.name  = sprintf('IMF1 -> A0-part (mask, \\sigma=%.1f\\mus, \\delta=%.1f\\mus)', sepOpt.imf1_sigma_us, sepOpt.imf1_delta_us);
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
    one.raw_sig     = raw_sig(:);
    one.S0_recon    = S0_recon(:);
    one.A0_recon    = A0_recon(:);

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
    toks = regexp(lower(fname),'(\d+)\s*mm','tokens');
    if isempty(toks), return; end
    vals = nan(numel(toks),1);
    for ii = 1:numel(toks)
        vals(ii) = str2double(toks{ii}{1});
    end
    vals = vals(isfinite(vals));
    if ~isempty(vals)
        dist_mm = max(vals);   % e.g. 3mmPMMA_100mm_edge -> 100 mm, not 3 mm
    end
end

function s = strip_ext(fname)
    [~,s,~] = fileparts(fname);
end

function label = parse_source_label(fname, srcOpt)
    label = srcOpt.surface_default;  % default: no 'edge' => surface
    if ~srcOpt.parse_from_filename, return; end
    s = lower(string(fname));
    if any(cellfun(@(k) contains(s, lower(string(k))), srcOpt.edge_keywords))
        label = 'edge';
    end
end

function shift_us = manual_shift_us_from_distance(dist_mm, curveOpt)
    shift_us = curveOpt.shift_us_fallback;
    if ~isfinite(dist_mm)
        return;
    end
    if isfield(curveOpt,'use_manual_shift_map') && curveOpt.use_manual_shift_map && ...
            isfield(curveOpt,'shift_map_dist_mm') && isfield(curveOpt,'shift_map_us') && ...
            numel(curveOpt.shift_map_dist_mm) == numel(curveOpt.shift_map_us) && numel(curveOpt.shift_map_us) >= 2
        try
            shift_us = interp1(curveOpt.shift_map_dist_mm(:), curveOpt.shift_map_us(:), dist_mm, 'linear', 'extrap');
        catch
            shift_us = curveOpt.shift_us_fallback;
        end
    end
end

function names = extractfeati_feature_names()
    names = {'hit','hittime_s','channel', ...
        'A_dB','D_us','E','ZCR','RT_us','TC_us','alpha', ...
        'PP2_1','PP2_2','PP2_3','PP2_4', ...
        'FC2_kHz','PF2_kHz','SSpread_kHz','SSkew','SKurt','SSlope','SRoff_kHz', ...
        'SSpreadP_sqrtkHz','SSkewP','SKurtP','SRon_kHz', ...
        'WPE1','WPE2','WPE3','WPE4','WPE5','WPE6','WPE7','WPE8','Entropy'};
end

function feat = extractfeati_struct(t_us, sig, fs, chSel, evSel)
    names = extractfeati_feature_names();
    vals = nan(1, numel(names));
    feat = cell2struct(num2cell(vals), names, 2);

    if exist('ExtractFeati_OPT','file') ~= 2
        warning('ExtractFeati_OPT 不在当前 MATLAB 路径中。');
        return;
    end

    x = double(sig(:));
    t_sec = double(t_us(:)) * 1e-6;
    [t_sec, x] = sanitize_signal_for_extractfeati(t_sec, x);
    if numel(x) < 8
        return;
    end

    % heads_info format required by your ExtractFeati_OPT:
    % {fs, channel, hit, hittime}
    heads_info = {double(fs), double(chSel), double(evSel), 0};

    try
        Feat = ExtractFeati_OPT(t_sec, x, [], heads_info);
    catch ME1
        try
            % Fallback: same feature definition, but with guard rails for short / edge cases.
            Feat = ExtractFeati_OPT_safe_local(t_sec, x, [], heads_info);
        catch ME2
            warning('ExtractFeati failed | ch=%d ev=%d | repo: %s | fallback: %s', ...
                chSel, evSel, ME1.message, ME2.message);
            return;
        end
    end

    Feat = double(Feat(:))';
    n = min(numel(Feat), numel(names));
    for ii = 1:n
        feat.(names{ii}) = Feat(ii);
    end
end

function [T2, V2] = sanitize_signal_for_extractfeati(T2, V2)
    T2 = double(T2(:));
    V2 = double(V2(:));
    n = min(numel(T2), numel(V2));
    T2 = T2(1:n);
    V2 = V2(1:n);

    good = isfinite(T2) & isfinite(V2);
    T2 = T2(good);
    V2 = V2(good);
    if isempty(T2)
        T2 = 0;
        V2 = 0;
        return;
    end

    % ensure time starts at 0 and is strictly nondecreasing
    T2 = T2 - T2(1);
    [T2, ia] = unique(T2, 'stable');
    V2 = V2(ia);

    if numel(T2) >= 2
        dt = median(diff(T2));
        if ~isfinite(dt) || dt <= 0
            dt = 1;
            T2 = (0:numel(V2)-1)' * dt;
        end
    else
        T2 = [0; 1e-6];
        V2 = [V2(1); V2(1)];
    end

    V2(~isfinite(V2)) = 0;
    if all(abs(V2) < eps)
        % keep tiny nonzero to avoid degenerate divisions inside the original code
        V2(1) = eps;
    end
end

function Feat = ExtractFeati_OPT_safe_local(T2,V2,E2,heads_info,varargin)
    %#ok<INUSD>
    if nargin < 4
        error('ExtractFeati_OPT_safe_local needs T2, V2, E2, heads_info');
    end
    [T2, V2] = sanitize_signal_for_extractfeati(T2, V2);

    fMAX = 1e6;
    intf = 1e3*[0 100 250 500 1000]';
    roll_off_factor = 0.95;
    roll_on_factor  = 0.05;

    fs = double(heads_info{1,1});
    channel = double(heads_info{1,2});
    hit = double(heads_info{1,3});
    hittime = double(heads_info{1,4});

    % -------- A. TIME FEATURES --------
    [A, bA] = max(abs(V2));
    if isempty(A) || ~isfinite(A), A = 0; end
    D = T2(end);
    LV2 = numel(V2);

    if ~isempty(E2)
        nE2 = numel(E2);
        idx = min(max(LV2,1), nE2);
        E = double(E2(idx));
    else
        E = sum(V2.^2, 'omitnan');
    end

    if LV2 >= 2
        c = (V2(1:LV2-1)<=0) & (V2(2:LV2)>0);
        ZC = sum(c);
    else
        ZC = 0;
    end
    ZCR = 100*ZC/max(LV2,1);
    RT = T2(min(max(bA,1),numel(T2)));

    V2_RMS = abs(V2) / sqrt(max(LV2,1));
    denTC = sum(V2_RMS);
    if denTC > 0
        TC = sum(T2.*V2_RMS)./denTC;
    else
        TC = NaN;
    end

    [~, bR] = max(V2_RMS);
    if numel(T2) - bR >= 2
        p = polyfit(T2(bR+1:end), V2_RMS(bR+1:end), 1);
        alpha = -p(1);
    else
        alpha = NaN;
    end

    % -------- B. FREQUENCY FEATURES --------
    NFFT = 2^nextpow2(max(LV2, 8));
    Yc = fft(V2, NFFT) / max(LV2,1);
    f = fs/2 * linspace(0,1,NFFT/2+1);

    Fc = f(f<=fMAX);
    Y  = abs(Yc(f<=fMAX)).';
    if isempty(Fc)
        Fc = 0; Y = 0;
    end

    Ptot = sum(Y);
    PP2 = zeros(length(intf)-1,1);
    if Ptot > 0
        for i=1:length(intf)-1
            PP2(i) = 100*sum(Y(Fc>=intf(i) & Fc<intf(i+1))) / Ptot;
        end
        FC2 = sum(Fc.*Y)/sum(Y);
        [~, iPF] = max(Y);
        PF2 = Fc(iPF);
        SSpread = sqrt(sum(((Fc-FC2).^2).*Y)/sum(Y));
        if isfinite(SSpread) && SSpread > eps
            SSkew = (1/SSpread^(3))*(sum(((Fc-FC2).^3).*Y)/sum(Y));
            SKurt = (1/SSpread^(4))*(sum(((Fc-FC2).^4).*Y)/sum(Y));
        else
            SSkew = NaN; SKurt = NaN;
        end
        if max(Y) > 0 && numel(Fc) >= 2
            p = polyfit(Fc/fMAX, Y/max(Y), 1);
            SSlope = p(1);
        else
            SSlope = NaN;
        end
        cum_energy = cumsum(Y);
        Roff = roll_off_factor * max(cum_energy);
        indroff = find(cum_energy < Roff, 1, 'last');
        if isempty(indroff), indroff = numel(Fc); end
        SRoff = Fc(indroff);
        SSpreadP = sqrt(sum(((Fc-PF2).^2).*Y)/sum(Y));
        if isfinite(SSpreadP) && SSpreadP > eps
            SSkewP = (1/SSpreadP^(3))*(sum(((Fc-PF2).^3).*Y)/sum(Y));
            SKurtP = (1/SSpreadP^(4))*(sum(((Fc-PF2).^4).*Y)/sum(Y));
        else
            SSkewP = NaN; SKurtP = NaN;
        end
        Ron = roll_on_factor * max(cum_energy);
        indron = find(cum_energy < Ron, 1, 'last');
        if isempty(indron), indron = 1; end
        SRon = Fc(indron);
    else
        FC2 = NaN; PF2 = NaN; SSpread = NaN; SSkew = NaN; SKurt = NaN; SSlope = NaN;
        SRoff = NaN; SSpreadP = NaN; SSkewP = NaN; SKurtP = NaN; SRon = NaN;
    end

    % Unit conversion consistent with your ExtractFeati_OPT
    A = 20.*(log10(max(A,eps).*1e6))-42;
    D  = D.*1e6;
    RT = RT.*1e6;
    TC = TC.*1e6;
    FC2 = FC2./1000;
    PF2 = PF2./1000;
    SRoff = SRoff./1000;
    SRon = SRon./1000;
    SSpread = SSpread./1000;
    SSpreadP = SSpreadP./1000;

    % -------- C. WAVELET PACKET FEATURES --------
    Ewp = nan(8,1);
    try
        Twp = wpdec(V2,3,'sym8','shannon');
        Ewp0 = wenergy(Twp);
        Ewp0 = Ewp0(:);
        m = min(numel(Ewp0), 8);
        Ewp(1:m) = Ewp0(1:m);
    catch
        % keep NaN if wavelet packet part fails
    end

    [N,~] = histcounts(abs(V2),100);
    L = N(N~=0);
    if isempty(L)
        Entropy = NaN;
    else
        L = L./sum(L);
        Entropy = -sum(L.*log2(L));
    end

    Feat = [hit hittime channel A D E ZCR RT TC alpha, ...
        PP2' FC2 PF2 SSpread SSkew SKurt SSlope SRoff, ...
        sqrt(SSpreadP) SSkewP SKurtP SRon Ewp' Entropy]';
end

function c = feature_struct_to_cell(feat, names)
    c = cell(1, numel(names));
    for ii = 1:numel(names)
        if isfield(feat, names{ii})
            c{ii} = feat.(names{ii});
        else
            c{ii} = NaN;
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

    if isempty(imfs)
        error('EMD 返回空。');
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
    if f(1) < f(end)
        f = flipud(f);
        wt = flipud(wt);
    end

    % Enforce strictly decreasing frequency vector for R2020b icwt()
    [f, idxSort] = sort(f(:), 'descend');
    wt = wt(idxSort, :);
    df = diff(f);
    if any(df >= 0)
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

    % icwt compatibility across MATLAB versions:
    % Use explicit wavelet name first, because some versions interpret a
    % numeric second positional input as wname and error out.
    freqrange = sort([min(f) max(f)], 'ascend');
    wtA = wt .* maskA;
    wtS = wt .* maskS;
    try
        xA0 = icwt(wtA, 'amor', f, freqrange, 'SignalMean', 0);
        xS0 = icwt(wtS, 'amor', f, freqrange, 'SignalMean', 0);
    catch ME1
        try
            xA0 = icwt(wtA, [], f, freqrange, 'SignalMean', 0);
            xS0 = icwt(wtS, [], f, freqrange, 'SignalMean', 0);
        catch ME2
            try
                xA0 = icwt(wtA, f, freqrange, 'SignalMean', 0);
                xS0 = icwt(wtS, f, freqrange, 'SignalMean', 0);
            catch ME3
                error('icwt reconstruction failed.\n  try1: %s\n  try2: %s\n  try3: %s', ...
                    ME1.message, ME2.message, ME3.message);
            end
        end
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
        if dispOpt.force_same_rawrow && isfinite(ref_max_raw) && ref_max_raw > 0
            norm_base = ref_max_raw;
        else
            norm_base = all_max;
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
    x0 = xl(2) - 0.18*(xl(2)-xl(1));
    y0 = yl(2) - 0.08*(yl(2)-yl(1));
    dy = 0.18*(yl(2)-yl(1));
    hold on;
    plot([x0, x0+0.04*(xl(2)-xl(1))], [y0, y0], curveOpt.A0_style, ...
        'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
    text(x0+0.05*(xl(2)-xl(1)), y0, 'A0', 'Color', curveOpt.color, ...
        'FontWeight','bold', 'VerticalAlignment','middle', 'HorizontalAlignment','left');
    y1 = y0 - dy;
    plot([x0, x0+0.04*(xl(2)-xl(1))], [y1, y1], curveOpt.S0_style, ...
        'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
    text(x0+0.05*(xl(2)-xl(1)), y1, 'S0', 'Color', curveOpt.color, ...
        'FontWeight','bold', 'VerticalAlignment','middle', 'HorizontalAlignment','left');
    hold off;
end

%% =====================================================================
%% ======================= BASIC FEATURE DIAGNOSTICS V9 ================
function plot_basic_diagnostics_v9(T, rankTab, selectedTab, outDir, diagDir, diagOpt, basicOpt)
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if ~exist(diagDir, 'dir'), mkdir(diagDir); end

    G = build_group_style(T);
    writetable(G.group_table, fullfile(outDir, 'Legend_groups.csv'));

    % 1) exact replacement of the old logRE-logRP overview:
    if ismember('logRE', T.Properties.VariableNames) && ismember('logRP', T.Properties.VariableNames)
        fig1 = figure('Color','w', 'Position', [60 60 980 760]);
        ax = axes(fig1); hold(ax, 'on');
        valid = isfinite(double(T.logRE)) & isfinite(double(T.logRP));
        draw_grouped_scatter(ax, double(T.logRE), double(T.logRP), G, valid, diagOpt.point_size);
        grid(ax, 'on');
        xlabel(ax, 'log_{10}(E_S / E_A)');
        ylabel(ax, 'log_{10}(P_S / P_A)');
        title(ax, 'Unified modal signature: logRE vs logRP', 'FontWeight','bold');
        add_group_legend(ax, G, diagOpt.legend_location, diagOpt.max_legend_columns);
        save_fig(fig1, fullfile(outDir, 'Diag_logRE_logRP_discreteLegend.png'), 300);
        close(fig1);
    end

    % 2) core pair scatter figures (same style as the sample figure)
    plot_core_pair_scatter_pages(T, G, diagDir, diagOpt);

    % 3) all numeric features vs distance (same discrete legend style, paged)
    plot_all_features_vs_distance_pages(T, G, diagDir, diagOpt, basicOpt);

    % 4) keep ranking tables already saved; additionally save selected-only quick view
    if ~isempty(selectedTab)
        try
            plot_selected_feature_pages(T, G, selectedTab, diagDir, diagOpt);
        catch ME
            warning('selectedTab quick view failed: %s', ME.message);
        end
    end
end

function G = build_group_style(T)
    n = height(T);
    src = repmat("unknown", n, 1);
    dist = nan(n,1);

    if ismember('source_label', T.Properties.VariableNames)
        src = lower(string(T.source_label));
    end
    if ismember('distance_mm', T.Properties.VariableNames)
        dist = double(T.distance_mm);
    end

    labels = strings(n,1);
    for i = 1:n
        if isfinite(dist(i))
            labels(i) = sprintf('%s %g mm', char(src(i)), dist(i));
        else
            labels(i) = char(src(i));
        end
    end
    labels = strrep(labels, 'unknown ', '');
    [groupNames, ~, gid] = unique(labels, 'stable');

    baseColors = [ ...
        0.1216 0.4667 0.7059; ...
        1.0000 0.4980 0.0549; ...
        0.1725 0.6275 0.1725; ...
        0.8392 0.1529 0.1569; ...
        0.5804 0.4039 0.7412; ...
        0.5490 0.3373 0.2941; ...
        0.8902 0.4667 0.7608; ...
        0.4980 0.4980 0.4980; ...
        0.7373 0.7412 0.1333; ...
        0.0902 0.7451 0.8118];
    markers = {'o','s','^','d','v','>','<','p','h','x','+','*'};

    nG = numel(groupNames);
    styles = struct('color', cell(1,nG), 'marker', cell(1,nG), 'label', cell(1,nG));
    rows = cell(nG, 4);
    for g = 1:nG
        styles(g).color = baseColors(mod(g-1, size(baseColors,1))+1, :);
        styles(g).marker = markers{mod(g-1, numel(markers))+1};
        styles(g).label = char(groupNames(g));
        rows(g,:) = {g, char(groupNames(g)), mat2str(styles(g).color,4), styles(g).marker};
    end

    G = struct();
    G.labels = groupNames;
    G.gid = gid;
    G.styles = styles;
    G.group_table = cell2table(rows, 'VariableNames', {'group_id','group_label','rgb','marker'});
end

function draw_grouped_scatter(ax, x, y, G, validMask, pointSize)
    if nargin < 5 || isempty(validMask), validMask = true(size(x)); end
    if nargin < 6 || isempty(pointSize), pointSize = 60; end

    x = double(x(:)); y = double(y(:));
    validMask = validMask(:) & isfinite(x) & isfinite(y);
    hold(ax, 'on');
    for g = 1:numel(G.labels)
        idx = validMask & (G.gid(:) == g);
        if ~any(idx), continue; end
        scatter(ax, x(idx), y(idx), pointSize, ...
            'Marker', G.styles(g).marker, ...
            'MarkerFaceColor', G.styles(g).color, ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 0.9);
    end
end

function hLeg = add_group_legend(ax, G, legendLocation, ncol)
    if nargin < 3 || isempty(legendLocation), legendLocation = 'eastoutside'; end
    if nargin < 4 || isempty(ncol), ncol = 1; end
    hold(ax, 'on');
    h = gobjects(numel(G.labels),1);
    names = cell(numel(G.labels),1);
    for g = 1:numel(G.labels)
        h(g) = plot(ax, nan, nan, ...
            'LineStyle','none', ...
            'Marker', G.styles(g).marker, ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', G.styles(g).color, ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 0.9);
        names{g} = G.styles(g).label;
    end
    hLeg = legend(ax, h, names, 'Location', legendLocation, 'Box','off');
    try
        hLeg.NumColumns = ncol;
    catch
    end
end

function plot_core_pair_scatter_pages(T, G, diagDir, diagOpt)
    pairs = diagOpt.core_pairs;
    nPairs = numel(pairs);
    if nPairs == 0, return; end

    nPerFig = diagOpt.tile_per_fig;
    nFig = ceil(nPairs / nPerFig);

    for fi = 1:nFig
        idx1 = (fi-1)*nPerFig + 1;
        idx2 = min(fi*nPerFig, nPairs);
        thisIdx = idx1:idx2;
        nThis = numel(thisIdx);

        fig = figure('Color','w', 'Position', [60 60 1500 max(650, 340*ceil(nThis/2))]);
        tl = tiledlayout(ceil(nThis/2), 2, 'TileSpacing','compact', 'Padding','compact');

        firstAx = [];
        for ii = 1:nThis
            p = pairs{thisIdx(ii)};
            xName = p{1}; yName = p{2}; ttl = p{3}; xLab = p{4}; yLab = p{5};
            if ~ismember(xName, T.Properties.VariableNames) || ~ismember(yName, T.Properties.VariableNames)
                continue;
            end
            ax = nexttile; hold(ax, 'on');
            x = double(T.(xName));
            y = double(T.(yName));
            valid = isfinite(x) & isfinite(y);
            draw_grouped_scatter(ax, x, y, G, valid, diagOpt.point_size);
            grid(ax, 'on');
            xlabel(ax, xLab, 'Interpreter','tex');
            ylabel(ax, yLab, 'Interpreter','tex');
            title(ax, ttl, 'Interpreter','none', 'FontWeight','bold');
            if isempty(firstAx), firstAx = ax; end
        end

        if ~isempty(firstAx)
            add_group_legend(firstAx, G, diagOpt.legend_location, diagOpt.max_legend_columns);
        end
        title(tl, 'Core modal feature pairs (discrete legend by source-distance group)', 'FontWeight','bold');
        save_fig(fig, fullfile(diagDir, sprintf('Diag_core_pair_scatter_page%02d.png', fi)), 300);
        close(fig);
    end
end

function plot_all_features_vs_distance_pages(T, G, diagDir, diagOpt, basicOpt)
    featNames = basicOpt.feature_names;
    featNames = featNames(ismember(featNames, T.Properties.VariableNames));
    featNames = featNames(:)';

    % keep numeric scalar columns only
    keep = false(size(featNames));
    for i = 1:numel(featNames)
        x = T.(featNames{i});
        keep(i) = isnumeric(x) || islogical(x);
    end
    featNames = featNames(keep);
    if diagOpt.feature_pages_only_valid
        keep2 = false(size(featNames));
        for i = 1:numel(featNames)
            x = double(T.(featNames{i}));
            keep2(i) = any(isfinite(x));
        end
        featNames = featNames(keep2);
    end
    if isempty(featNames), return; end

    d = double(T.distance_mm);
    if ~any(isfinite(d)), return; end

    nPerFig = diagOpt.tile_per_fig;
    nFig = ceil(numel(featNames) / nPerFig);
    dspan = max(d(isfinite(d))) - min(d(isfinite(d)));
    if ~isfinite(dspan) || dspan <= 0, dspan = 1; end

    for fi = 1:nFig
        idx1 = (fi-1)*nPerFig + 1;
        idx2 = min(fi*nPerFig, numel(featNames));
        thisFeat = featNames(idx1:idx2);
        nThis = numel(thisFeat);

        fig = figure('Color','w', 'Position', [60 60 1550 max(700, 330*ceil(nThis/2))]);
        tl = tiledlayout(ceil(nThis/2), 2, 'TileSpacing','compact', 'Padding','compact');

        firstAx = [];
        for ii = 1:nThis
            fname = thisFeat{ii};
            ax = nexttile; hold(ax, 'on');
            y = double(T.(fname));
            valid = isfinite(d) & isfinite(y);

            for g = 1:numel(G.labels)
                idx = valid & (G.gid(:) == g);
                if ~any(idx), continue; end
                xg = d(idx) + local_group_jitter(g, numel(G.labels), dspan, diagOpt.jitter_frac);
                scatter(ax, xg, y(idx), diagOpt.point_size, ...
                    'Marker', G.styles(g).marker, ...
                    'MarkerFaceColor', G.styles(g).color, ...
                    'MarkerEdgeColor', [0 0 0], ...
                    'LineWidth', 0.9);
            end

            grid(ax, 'on');
            xlabel(ax, 'Distance [mm]');
            ylabel(ax, fname, 'Interpreter','none');
            title(ax, sprintf('%s vs distance', fname), 'Interpreter','none', 'FontWeight','bold');
            if isempty(firstAx), firstAx = ax; end
        end

        if ~isempty(firstAx)
            add_group_legend(firstAx, G, diagOpt.legend_location, diagOpt.max_legend_columns);
        end
        title(tl, sprintf('All numeric features vs distance (page %d/%d)', fi, nFig), 'FontWeight','bold');
        save_fig(fig, fullfile(diagDir, sprintf('Diag_all_features_vs_distance_page%02d.png', fi)), 300);
        close(fig);
    end
end

function plot_selected_feature_pages(T, G, selectedTab, diagDir, diagOpt)
    featNames = selectedTab.feature(:)';
    if isstring(featNames), featNames = cellstr(featNames); end
    featNames = featNames(ismember(featNames, T.Properties.VariableNames));
    keep = false(size(featNames));
    for i = 1:numel(featNames)
        x = double(T.(featNames{i}));
        keep(i) = any(isfinite(x));
    end
    featNames = featNames(keep);
    if isempty(featNames), return; end

    d = double(T.distance_mm);
    if ~any(isfinite(d)), return; end
    dspan = max(d(isfinite(d))) - min(d(isfinite(d)));
    if ~isfinite(dspan) || dspan <= 0, dspan = 1; end

    nPerFig = diagOpt.tile_per_fig;
    nFig = ceil(numel(featNames) / nPerFig);

    for fi = 1:nFig
        idx1 = (fi-1)*nPerFig + 1;
        idx2 = min(fi*nPerFig, numel(featNames));
        thisFeat = featNames(idx1:idx2);
        nThis = numel(thisFeat);

        fig = figure('Color','w', 'Position', [60 60 1550 max(700, 330*ceil(nThis/2))]);
        tl = tiledlayout(ceil(nThis/2), 2, 'TileSpacing','compact', 'Padding','compact');
        firstAx = [];

        for ii = 1:nThis
            fname = thisFeat{ii};
            ax = nexttile; hold(ax, 'on');
            y = double(T.(fname));
            valid = isfinite(d) & isfinite(y);

            for g = 1:numel(G.labels)
                idx = valid & (G.gid(:) == g);
                if ~any(idx), continue; end
                xg = d(idx) + local_group_jitter(g, numel(G.labels), dspan, diagOpt.jitter_frac);
                scatter(ax, xg, y(idx), diagOpt.point_size, ...
                    'Marker', G.styles(g).marker, ...
                    'MarkerFaceColor', G.styles(g).color, ...
                    'MarkerEdgeColor', [0 0 0], ...
                    'LineWidth', 0.9);
            end

            grid(ax, 'on');
            xlabel(ax, 'Distance [mm]');
            ylabel(ax, fname, 'Interpreter','none');
            title(ax, sprintf('Selected feature: %s', fname), 'Interpreter','none', 'FontWeight','bold');
            if isempty(firstAx), firstAx = ax; end
        end

        if ~isempty(firstAx)
            add_group_legend(firstAx, G, diagOpt.legend_location, diagOpt.max_legend_columns);
        end
        title(tl, sprintf('Selected high-score features (page %d/%d)', fi, nFig), 'FontWeight','bold');
        save_fig(fig, fullfile(diagDir, sprintf('Diag_selected_features_page%02d.png', fi)), 300);
        close(fig);
    end
end

function dx = local_group_jitter(g, nG, dspan, frac)
    if nargin < 4 || isempty(frac), frac = 0.03; end
    if nG <= 1
        dx = 0;
    else
        offsets = linspace(-0.5, 0.5, nG);
        dx = offsets(g) * frac * dspan;
    end
end

function [rankTab, selectedTab] = rank_modal_features(T, extractNamesNoMeta, topN)
    rows = {};
    for j = 1:numel(extractNamesNoMeta)
        base = extractNamesNoMeta{j};
        sName = ['S0_' base];
        aName = ['A0_' base];
        if ~ismember(sName, T.Properties.VariableNames) || ~ismember(aName, T.Properties.VariableNames)
            continue;
        end
        xS = double(T.(sName));
        xA = double(T.(aName));
        dmm = double(T.distance_mm);
        good = isfinite(xS) & isfinite(xA) & isfinite(dmm);
        if nnz(good) < 4
            continue;
        end
        xS = xS(good); xA = xA(good); dg = dmm(good);
        mdiff = mean(xS - xA, 'omitnan');
        pooled = sqrt((var(xS,0,'omitnan') + var(xA,0,'omitnan'))/2 + eps);
        modal_effect_d = abs(mdiff) / pooled;
        sign_consistency = max(mean((xS - xA) > 0), mean((xS - xA) < 0));
        dist_eta = 0.5 * (oneway_eta2(xS, dg) + oneway_eta2(xA, dg));
        score = modal_effect_d * sign_consistency / (1 + 1.5*dist_eta);
        rows(end+1,:) = {base, modal_effect_d, sign_consistency, dist_eta, mdiff, score}; %#ok<AGROW>
    end
    if isempty(rows)
        rankTab = table(); selectedTab = table(); return;
    end
    rankTab = cell2table(rows, 'VariableNames', ...
        {'feature_base','modal_effect_d','sign_consistency','distance_eta2_mean','mean_S0_minus_A0','modal_score'});
    rankTab = sortrows(rankTab, {'modal_score','modal_effect_d'}, {'descend','descend'});
    nSel = min(topN, height(rankTab));
    selectedTab = rankTab(1:nSel,:);
end

function plot_mode_specific_feature_pages(T, diagDir, diagOpt, extractNamesNoMeta, modeTag)
    if nargin < 5, modeTag = 'S0'; end
    prefix = [modeTag '_'];
    featNames = strcat(prefix, extractNamesNoMeta(:)');
    featNames = featNames(ismember(featNames, T.Properties.VariableNames));
    if isempty(featNames), return; end

    d = double(T.distance_mm);
    src = lower(string(T.source_label));
    validSrc = unique(src, 'stable');
    validSrc(validSrc=="") = [];
    if isempty(validSrc), validSrc = "unknown"; end

    keep = false(size(featNames));
    for i = 1:numel(featNames)
        x = double(T.(featNames{i}));
        keep(i) = any(isfinite(x));
    end
    featNames = featNames(keep);
    if isempty(featNames), return; end

    nPerFig = diagOpt.tile_per_fig;
    nFig = ceil(numel(featNames)/nPerFig);
    for fi = 1:nFig
        idx1 = (fi-1)*nPerFig + 1;
        idx2 = min(fi*nPerFig, numel(featNames));
        thisFeat = featNames(idx1:idx2);
        nThis = numel(thisFeat);
        fig = figure('Color','w', 'Position', [60 60 1550 max(700, 330*ceil(nThis/2))]);
        tl = tiledlayout(ceil(nThis/2), 2, 'TileSpacing','compact', 'Padding','compact');

        for ii = 1:nThis
            fname = thisFeat{ii};
            ax = nexttile; hold(ax,'on');
            y = double(T.(fname));
            for s = 1:numel(validSrc)
                idx = isfinite(d) & isfinite(y) & (src == validSrc(s));
                if ~any(idx), continue; end
                if validSrc(s) == "edge"
                    mk = '^';
                else
                    mk = 'o';
                end
                scatter(ax, d(idx), y(idx), diagOpt.point_size, ...
                    'Marker', mk, 'MarkerFaceColor',[0.2 0.6 0.85], 'MarkerEdgeColor','k', 'LineWidth',0.9);
            end
            grid(ax,'on');
            xlabel(ax, 'Distance [mm]');
            ylabel(ax, fname, 'Interpreter','none');
            title(ax, sprintf('%s feature vs distance', fname), 'Interpreter','none', 'FontWeight','bold');
        end
        title(tl, sprintf('%s-only ExtractFeati features vs distance (page %d/%d)', modeTag, fi, nFig), 'FontWeight','bold');
        save_fig(fig, fullfile(diagDir, sprintf('Diag_%s_features_vs_distance_page%02d.png', modeTag, fi)), 300);
        close(fig);
    end
end

function plot_modal_selected_feature_pages(T, diagDir, diagOpt, modalSelectedTab)
    if isempty(modalSelectedTab), return; end
    featBase = modalSelectedTab.feature_base(:)';
    if isstring(featBase), featBase = cellstr(featBase); end
    nPerFig = max(2, floor(diagOpt.tile_per_fig/2));
    nFig = ceil(numel(featBase) / nPerFig);
    d = double(T.distance_mm);
    G = build_group_style(T);

    for fi = 1:nFig
        idx1 = (fi-1)*nPerFig + 1;
        idx2 = min(fi*nPerFig, numel(featBase));
        thisFeat = featBase(idx1:idx2);
        nThis = numel(thisFeat);
        fig = figure('Color','w', 'Position', [60 60 1550 max(720, 300*nThis)]);
        tl = tiledlayout(nThis, 2, 'TileSpacing','compact', 'Padding','compact');
        firstAx = [];
        for ii = 1:nThis
            base = thisFeat{ii};
            sName = ['S0_' base];
            aName = ['A0_' base];
            if ~ismember(sName, T.Properties.VariableNames) || ~ismember(aName, T.Properties.VariableNames)
                continue;
            end
            % left: S0 vs distance
            ax1 = nexttile; hold(ax1,'on');
            y1 = double(T.(sName));
            valid1 = isfinite(d) & isfinite(y1);
            for g = 1:numel(G.labels)
                idx = valid1 & (G.gid(:)==g);
                if ~any(idx), continue; end
                scatter(ax1, d(idx), y1(idx), diagOpt.point_size, 'Marker', G.styles(g).marker, ...
                    'MarkerFaceColor', G.styles(g).color, 'MarkerEdgeColor','k', 'LineWidth',0.9);
            end
            grid(ax1,'on'); xlabel(ax1,'Distance [mm]'); ylabel(ax1, sName, 'Interpreter','none');
            title(ax1, sprintf('S0: %s', base), 'Interpreter','none', 'FontWeight','bold');
            if isempty(firstAx), firstAx = ax1; end

            % right: A0 vs distance
            ax2 = nexttile; hold(ax2,'on');
            y2 = double(T.(aName));
            valid2 = isfinite(d) & isfinite(y2);
            for g = 1:numel(G.labels)
                idx = valid2 & (G.gid(:)==g);
                if ~any(idx), continue; end
                scatter(ax2, d(idx), y2(idx), diagOpt.point_size, 'Marker', G.styles(g).marker, ...
                    'MarkerFaceColor', G.styles(g).color, 'MarkerEdgeColor','k', 'LineWidth',0.9);
            end
            grid(ax2,'on'); xlabel(ax2,'Distance [mm]'); ylabel(ax2, aName, 'Interpreter','none');
            title(ax2, sprintf('A0: %s', base), 'Interpreter','none', 'FontWeight','bold');
        end
        if ~isempty(firstAx)
            add_group_legend(firstAx, G, diagOpt.legend_location, diagOpt.max_legend_columns);
        end
        title(tl, sprintf('Top modal-discriminative ExtractFeati features (page %d/%d)', fi, nFig), 'FontWeight','bold');
        save_fig(fig, fullfile(diagDir, sprintf('Diag_modal_selected_features_page%02d.png', fi)), 300);
        close(fig);
    end
end

%% =====================================================================
%% ======================= UNIFIED SOURCE IDENTIFICATION ================
% Reserved from previous versions if later needed. Not called in this script.

function [T, model, cvres] = build_unified_source_model(T, srcOpt)
    model = struct();
    cvres = struct('overall_accuracy', NaN, 'n_eval', 0);

    if ~ismember('source_label', T.Properties.VariableNames)
        return;
    end

    known = ismember(lower(string(T.source_label)), {'surface','edge'});
    if nnz(known) < 2
        warning('已知 source_label 太少，跳过统一识别。');
        return;
    end

    featNames = srcOpt.feature_names;
    featNames = featNames(ismember(featNames, T.Properties.VariableNames));
    if isempty(featNames)
        warning('未找到统一识别所需特征列，跳过统一识别。');
        return;
    end

    Xall = table_to_numeric(T, featNames);
    [Xz, scaler] = robust_standardize(Xall);

    labels = lower(string(T.source_label));
    model.feature_names = featNames;
    model.scaler = scaler;

    classes = {'surface','edge'};
    dmat = nan(height(T), numel(classes));
    mu = cell(1, numel(classes));
    invS = cell(1, numel(classes));
    nClass = zeros(1, numel(classes));

    for c = 1:numel(classes)
        idx = known & strcmp(labels, classes{c});
        nClass(c) = nnz(idx);
        if nClass(c) < max(2, srcOpt.min_ref_per_class)
            warning('类别 %s 样本数不足（n=%d），但仍尝试建模。', classes{c}, nClass(c));
        end
        [mu{c}, invS{c}] = fit_mahal_model(Xz(idx,:));
        dmat(:,c) = mahal_predict(Xz, mu{c}, invS{c});
    end

    model.classes = classes;
    model.mu = mu;
    model.invS = invS;
    model.nClass = nClass;

    T.score_surface = dmat(:,1);
    T.score_edge    = dmat(:,2);

    pred = repmat("unknown", height(T), 1);
    pred_margin = nan(height(T),1);
    pred_conf = nan(height(T),1);

    for i = 1:height(T)
        if ~all(isfinite(dmat(i,:)))
            continue;
        end
        [dmin, ic] = min(dmat(i,:));
        dother = dmat(i, 3-ic);
        pred(i) = string(classes{ic});
        pred_margin(i) = dother - dmin;
        pred_conf(i) = 1 ./ (1 + exp(-(dother - dmin)));
    end

    T.pred_source_label = cellstr(pred);
    T.pred_margin = pred_margin;
    T.pred_conf = pred_conf;

    if srcOpt.do_leave_one_distance_out
        cvres = leave_one_distance_out_cv(T, featNames, srcOpt);
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
        if iscell(x)
            x = cellfun(@double, x);
        end
        X(:,j) = double(x);
    end
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

function cvres = leave_one_distance_out_cv(T, featNames, srcOpt)
    cvres = struct('overall_accuracy', NaN, 'n_eval', 0, 'by_distance', table());
    known = ismember(lower(string(T.source_label)), {'surface','edge'});
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
        if nnz(test) < 1 || nnz(train) < 2
            continue;
        end

        Xtr = table_to_numeric(T(train,:), featNames);
        Xte = table_to_numeric(T(test,:),  featNames);

        [Xtrz, scaler] = robust_standardize(Xtr);
        Xtez = apply_scaler(Xte, scaler);

        ytr = lower(string(T.source_label(train)));
        classes = {'surface','edge'};
        if ~all(ismember(classes, unique(cellstr(ytr))))
            continue;
        end

        dmat = nan(size(Xtez,1), 2);
        for c = 1:2
            idxc = strcmp(ytr, classes{c});
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

        yte = lower(string(T.source_label(test)));
        acc = mean(pred == yte);
        rows(end+1,:) = {d0, nnz(test), acc}; %#ok<AGROW>

        pred_all = [pred_all; pred(:)]; %#ok<AGROW>
        true_all = [true_all; yte(:)]; %#ok<AGROW>
        dist_log = [dist_log; repmat(d0, numel(yte), 1)]; %#ok<AGROW>
    end

    if ~isempty(rows)
        cvres.by_distance = cell2table(rows, 'VariableNames', {'distance_mm','n_test','accuracy'});
    end
    good = (pred_all == "surface" | pred_all == "edge") & (true_all == "surface" | true_all == "edge");
    if any(good)
        cvres.overall_accuracy = mean(pred_all(good) == true_all(good));
        cvres.n_eval = nnz(good);
        cvres.dist_log = dist_log(good);
        cvres.pred = cellstr(pred_all(good));
        cvres.truth = cellstr(true_all(good));
    end
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

function [rankTab, selectedTab] = rank_basic_features(T, basicOpt)
    featNames = basicOpt.feature_names;
    featNames = featNames(ismember(featNames, T.Properties.VariableNames));
    labels = lower(string(T.source_label));
    distv = double(T.distance_mm);

    rows = {};
    for j = 1:numel(featNames)
        x = double(T.(featNames{j}));
        good = isfinite(x) & isfinite(distv) & (labels=="surface" | labels=="edge");
        if nnz(good) < 4
            continue;
        end
        xg = x(good);
        yg = labels(good);
        dg = distv(good);

        xs = xg(yg=="surface");
        xe = xg(yg=="edge");

        if numel(xs) < 2 || numel(xe) < 2
            continue;
        end

        ms = mean(xs,'omitnan');
        me = mean(xe,'omitnan');
        sps = std(xs,0,'omitnan');
        spe = std(xe,0,'omitnan');
        sp = sqrt((sps.^2 + spe.^2)/2 + eps);
        source_effect_d = abs(me - ms) / sp;

        eta2_distance = oneway_eta2(xg, dg);
        drift_span = group_mean_span(xg, dg);
        simple_score = source_effect_d / (1 + 2*eta2_distance + 0.25*drift_span);

        rows(end+1,:) = {featNames{j}, source_effect_d, eta2_distance, drift_span, simple_score}; %#ok<AGROW>
    end

    if isempty(rows)
        rankTab = table(); selectedTab = table(); return;
    end

    rankTab = cell2table(rows, 'VariableNames', ...
        {'feature','source_effect_d','eta2_distance','drift_span','simple_score'});
    rankTab = sortrows(rankTab, {'simple_score','source_effect_d'}, {'descend','descend'});

    nSel = min(basicOpt.topN, height(rankTab));
    selectedTab = rankTab(1:nSel,:);
end

function eta2 = oneway_eta2(x, g)
    eta2 = NaN;
    x = x(:); g = g(:);
    good = isfinite(x) & isfinite(g);
    x = x(good); g = g(good);
    ug = unique(g);
    if numel(ug) < 2 || numel(x) < 4
        eta2 = 0; return;
    end
    mu = mean(x);
    ssb = 0; sst = sum((x - mu).^2);
    for k = 1:numel(ug)
        idx = (g == ug(k));
        nk = nnz(idx);
        if nk == 0, continue; end
        mk = mean(x(idx));
        ssb = ssb + nk * (mk - mu).^2;
    end
    eta2 = ssb / (sst + eps);
end

function sp = group_mean_span(x, g)
    x = x(:); g = g(:);
    good = isfinite(x) & isfinite(g);
    x = x(good); g = g(good);
    ug = unique(g);
    if isempty(ug), sp = NaN; return; end
    m = nan(numel(ug),1);
    for k = 1:numel(ug)
        m(k) = mean(x(g==ug(k)), 'omitnan');
    end
    den = std(x,0,'omitnan');
    if ~isfinite(den) || den < eps, den = 1; end
    sp = (max(m) - min(m)) / den;
end
