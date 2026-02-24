%% ================== Compare_EMD_EEMD_CEEMDAN_MAE.m ==================
% 目的：
%   单文件整合对比 EMD / EEMD / CEEMDAN：
%     1) 同一次选取 wave 数据（通道/事件 + 截取时间窗）
%     2) 三种方法分别分解 -> 自动 A0/S0 分组 -> 重构
%     3) 每种方法各自输出：
%           Figure: Decomposition (Original + IMFs, AGUlike)
%           Figure: Reconstruction (Original / S0 / A0, AGUlike)
%     4) 额外输出 1 张总对比图：
%           Figure: Comparison of Reconstructions (EMD vs EEMD vs CEEMDAN)
%
% 参考：
%   - Huang et al., 1998, Proc. R. Soc. A (EMD/HHT)
%   - Wu & Huang, 2009, Advances in Adaptive Data Analysis (EEMD)
%   - Torres et al., 2011, ICASSP (CEEMDAN)
%
% 依赖：
%   - waveReader（你已有）
%   - emd() + cwt()（Signal Processing Toolbox）
%
% 作者：为你当前三份脚本“合并重构”的单文件版本（便于对比与复现）
% =====================================================================

clear; clc; close all;

%% ================= 0) 用户可调参数（优先改这里） =================
% -------- (A) 时间窗：只分析前 focus_us 微秒 --------
timeOpt.focus_us = 500;       % 关注窗口长度 [us]
timeOpt.start_us = 0;         % 起始时刻 [us]

% -------- (B) 读 wave 的默认通道/事件（也会弹窗让你改）--------
dataOpt.sel_ch    = 1;
dataOpt.sel_event = 1;        % 若 wave 为 3D（N x ch x event）

% -------- (C) 色散曲线叠加（可选）--------
curveOpt.enable   = true;     % true: 叠加色散曲线；false: 不叠加
curveOpt.dist_mm  = 150;      % 距离 [mm]
curveOpt.shift_us = 44;       % 曲线整体平移 [us]（对齐触发/起点偏移）
curveOpt.lineWidth = 1.6;
curveOpt.A0_style  = '-';     % A0 实线
curveOpt.S0_style  = '--';    % S0 虚线
curveOpt.color     = [1 1 1]; % 白色

% -------- (D) CWT/显示参数（只影响显示清晰度）--------
dispOpt.fmax_kHz          = 400;
dispOpt.nFreqBins         = 400;
dispOpt.voicesPerOctave   = 32;
dispOpt.threshold         = 0.05;
dispOpt.gamma             = 1.3;
dispOpt.force_same_rawrow = true;   % Figure 里 Original 行归一化一致（跨图）
dispOpt.disable_local_norm = true; % false: 保留“弱行局部归一化”（更易看弱信号）
                                   % true : 强制所有行用全局归一化（更适合严格对比幅值）

% -------- (E) A0 / S0 自动分组 --------
% 说明：
%   - 优先使用“脊线匹配 + 理论色散曲线”进行 IMF 归类（推荐，需加载曲线文件）
%   - 若曲线不可用/匹配失败，则回退到 CWT 低频能量占比分组（旧策略）
sepOpt.method = 'ridge';       % 'ridge' / 'lowband'
sepOpt.split_freq_Hz  = 150e3; % 低频阈值（回退策略使用）
sepOpt.ratio_th       = 0.55;  % lowRatio >= ratio_th -> A0-like（回退策略）
sepOpt.sort_plot_by_centroid = true; % 仅影响 IMF 显示顺序（高频->低频），不改分解结果

% --- 脊线匹配参数（IMF dominant ridge vs A0/S0 curve）---
sepOpt.ridge.min_col_rel      = 0.08;   % 列峰值低于该比例则不参与脊线评分（相对该 IMF CWT 峰值）
sepOpt.ridge.smooth_cols      = 7;      % 脊线频率中值平滑窗口（列数，奇数更稳）
sepOpt.ridge.freq_tol_kHz     = 25;     % 频差容忍（kHz），用于 ridge match score
sepOpt.ridge.band_half_kHz    = 20;     % 曲线采样带宽（kHz），用于 curve-energy score
sepOpt.ridge.min_valid_cols   = 12;     % 至少有效列数，否则回退低频法
sepOpt.ridge.score_margin     = 0.03;   % A0/S0 分数过近时，回退低频法（避免硬判）
sepOpt.ridge.w_ridge          = 0.70;   % 最终分数 = w_ridge*ridgeMatch + (1-w_ridge)*curveEnergy
sepOpt.ridge.w_curve          = 0.30;   % 保留字段（便于结果记录）
% --- 分时段重构（解决“混合 IMF 被整条判给 S0/A0”的问题）---
sepOpt.ridge.segmented_recon   = true;   % true: 同一 IMF 允许分时段归入 S0/A0（推荐）
sepOpt.ridge.local_score_margin = 0.02;  % 局部列分数差阈值；过小更敏感，过大更保守
sepOpt.ridge.mask_smooth_cols   = 11;    % 时间掩膜平滑窗口（列数）
sepOpt.ridge.keep_unknown_to_global = true; % 局部不确定列按该 IMF 全局类别补齐
sepOpt.ridge.include_res_in_A0_diag = true; % 诊断用：是否把 residual 并到 A0（默认false）


% -------- (F) 三种方法的“输出 IMF 数”统一（便于对比）--------
Kmax = 8;

% -------- (G) 标准 EMD 参数（你原 EMD_MAE 里那套）--------
emdStdOpt.max_imf  = Kmax;
emdStdOpt.siftRelTol = 0.01;
emdStdOpt.interp  = 'pchip';     % 'spline'/'pchip'
emdStdOpt.maxEnergyRatio = 20;
emdStdOpt.siftMaxIter = 100;
emdStdOpt.maxNumExtrema = 1;
emdStdOpt.display = 0;
emdStdOpt.use_default_stop = false; % true: MATLAB 默认停止准则（更“标准”）
                                    % false: 使用你自定义的这些停止相关参数

% -------- (H) EEMD 参数（Wu & Huang 2009）--------
eemdOpt.max_imf         = Kmax;
eemdOpt.ensemble_size   = 100;     % 50~200 常用
eemdOpt.noise_std_ratio = 0.1;     % 0.1~0.3 常用
eemdOpt.rng_seed        = 1;       % 固定随机种子；设 [] 则每次不同
eemdOpt.paired_noise    = true;    % +noise/-noise 成对平均，残余噪声更小

% -------- (I) CEEMDAN 参数（Torres 2011）--------
ceemdanOpt.max_imf         = Kmax;
ceemdanOpt.ensemble_size   = 100;
ceemdanOpt.noise_std_ratio = 0.1;
ceemdanOpt.rng_seed        = 1;
ceemdanOpt.paired_noise    = true;
ceemdanOpt.early_stop      = false; % true: 残差极值太少提前停止
ceemdanOpt.min_extrema     = 2;

% -------- (J) EEMD/CEEMDAN 内部单次 EMD 参数 --------
% 重要：为公平比较，这里将“内部 EMD 参数”与标准 EMD 参数保持一致（同一套 sifting 参数）
%       实际生成 emdInnerOpt 的代码放在运行前（Section 4）自动构建。
emdInnerOpt = struct();  % placeholder（后面由 build_emd_inneropt_from_std() 填充）

% -------- (K) 输出控制 --------
plotOpt.per_method_figures = true;   % 每种方法各画 2 张图
plotOpt.comparison_figure  = true;   % 额外画 1 张“总对比重构图”
plotOpt.print_metrics      = true;   % 命令行输出对比指标
saveOpt.save_mat           = false;   % 保存结果 .mat（同目录）
saveOpt.prefix             = 'MAE_compare';

%% ================= 1) 选择 wave 文件 =================
fprintf('请选择数据文件（*.wave / *.tradb）...\n');
[fname, pname] = uigetfile({'*.wave;*.tradb;*.*'}, '选择 wave 文件');
if isequal(fname, 0), error('未选择 wave 文件'); end
wavefile = fullfile(pname, fname);

%% ================= 2) 读 wave -> raw_sig（并可交互选 ch/event） =================
[t0, wave, ~, ~, ~, ~, sf, ~, ~, ~, ~, ~, ~] = waveReader(wavefile);
fs = double(sf);

% 交互输入 ch / event（避免你每次手动改参数）
try
    if ndims(wave) == 3
        prompt = {'Channel index (1..ch):', 'Event index (1..event):'};
        defans = {num2str(dataOpt.sel_ch), num2str(dataOpt.sel_event)};
        answ = inputdlg(prompt, 'Select channel/event', 1, defans);
        if ~isempty(answ)
            dataOpt.sel_ch    = max(1, round(str2double(answ{1})));
            dataOpt.sel_event = max(1, round(str2double(answ{2})));
        end
        raw_sig = double(wave(:, dataOpt.sel_ch, dataOpt.sel_event));
    else
        prompt = {'Channel index (1..ch):'};
        defans = {num2str(dataOpt.sel_ch)};
        answ = inputdlg(prompt, 'Select channel', 1, defans);
        if ~isempty(answ)
            dataOpt.sel_ch = max(1, round(str2double(answ{1})));
        end
        raw_sig = double(wave(:, dataOpt.sel_ch));
    end
catch
    raw_sig = double(wave(:,1));
end

raw_sig(~isfinite(raw_sig)) = 0;
raw_sig = raw_sig(:) - mean(raw_sig);

% 截取时间窗
start_idx = max(1, round(timeOpt.start_us*1e-6*fs) + 1);
end_idx   = min(length(raw_sig), start_idx + round(timeOpt.focus_us*1e-6*fs) - 1);
raw_sig   = raw_sig(start_idx:end_idx);
N         = numel(raw_sig);
time_us   = (0:N-1)/fs * 1e6;

fprintf('Loaded: N=%d, fs=%.0f Hz, window=[%.1f, %.1f] us\n', ...
    N, fs, time_us(1), time_us(end));

%% ================= 3) 选择色散曲线文件（可选） =================
curve_data = struct();
if curveOpt.enable
    fprintf('请选择色散曲线文件（Excel/CSV，可取消）...\n');
    [cname, cpname] = uigetfile({'*.xlsx;*.xls;*.csv;*.txt;*.*'}, '选择色散曲线文件（可取消）');
    if isequal(cname, 0)
        warning('未选择曲线文件：将不叠加色散曲线。');
        curveOpt.enable = false;
    else
        curvefile = fullfile(cpname, cname);
        curve_data = load_dispersion_curve(curvefile, curveOpt.dist_mm, curveOpt.shift_us);
        if ~curve_data.ok
            warning('曲线文件读取失败：%s -> 不叠加曲线。', curve_data.msg);
            curveOpt.enable = false;
        end
    end
end

%% ================= 4) 三种方法：分解 + 分组 + 重构 =================
results = struct();

% 公平比较：EEMD/CEEMDAN 内部 EMD 与标准 EMD 使用同一组参数
emdInnerOpt = build_emd_inneropt_from_std(emdStdOpt, Kmax);

% ---- 4.1 标准 EMD ----
fprintf('\n=== Running Standard EMD ===\n');
[imfs_emd, res_emd] = emd_standard(raw_sig, emdStdOpt);
[resStruct_emd] = group_and_recon(raw_sig, imfs_emd, res_emd, fs, sepOpt, dispOpt, time_us, curve_data, curveOpt);
results.EMD = resStruct_emd;

% ---- 4.2 EEMD ----
fprintf('\n=== Running EEMD  ===\n');
[imfs_eemd, res_eemd, info_eemd] = eemd_decompose(raw_sig, eemdOpt, emdInnerOpt);
[resStruct_eemd] = group_and_recon(raw_sig, imfs_eemd, res_eemd, fs, sepOpt, dispOpt, time_us, curve_data, curveOpt);
resStruct_eemd.info = info_eemd;
results.EEMD = resStruct_eemd;

% ---- 4.3 CEEMDAN ----
fprintf('\n=== Running CEEMDAN ===\n');
[imfs_ce, res_ce, info_ce] = ceemdan_decompose(raw_sig, ceemdanOpt, emdInnerOpt);
[resStruct_ce] = group_and_recon(raw_sig, imfs_ce, res_ce, fs, sepOpt, dispOpt, time_us, curve_data, curveOpt);
resStruct_ce.info = info_ce;
results.CEEMDAN = resStruct_ce;

%% ================= 5) 输出对比指标（可选，但强烈建议） =================
if plotOpt.print_metrics
    fprintf('\n================= QUICK METRICS (for comparison) =================\n');
    print_method_metrics('EMD',     results.EMD,     raw_sig, fs);
    print_method_metrics('EEMD',    results.EEMD,    raw_sig, fs);
    print_method_metrics('CEEMDAN', results.CEEMDAN, raw_sig, fs);
    fprintf('==================================================================\n');
end

%% ================= 6) 画图：每种方法 2 张 + 总对比 1 张 =================
% 为了跨图保持 Original 行完全一致：共享 ref_max_raw / raw_ylim（只影响显示）
[ref_max_raw, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt);

if plotOpt.per_method_figures
    plot_method('EMD',     results.EMD,     time_us, fs, dispOpt, ref_max_raw, raw_ylim, curve_data, curveOpt, sepOpt);
    plot_method('EEMD',    results.EEMD,    time_us, fs, dispOpt, ref_max_raw, raw_ylim, curve_data, curveOpt, sepOpt);
    plot_method('CEEMDAN', results.CEEMDAN, time_us, fs, dispOpt, ref_max_raw, raw_ylim, curve_data, curveOpt, sepOpt);
end

if plotOpt.comparison_figure
    % 对比图建议禁用“弱行局部归一化”，让三方法在同一幅值基准下更可比
    dispCmp = dispOpt;
    dispCmp.disable_local_norm = true;

    list_cmp = cell(1, 9);
    % EMD
    list_cmp{1}.sig  = raw_sig;               list_cmp{1}.name  = 'EMD | Original (mix)';
    list_cmp{2}.sig  = results.EMD.S0_recon;  list_cmp{2}.name  = sprintf('EMD | S0 recon (IMFs %s)', int2str(results.EMD.idx_S0));
    list_cmp{3}.sig  = results.EMD.A0_recon;  list_cmp{3}.name  = sprintf('EMD | A0 recon (IMFs %s)', int2str(results.EMD.idx_A0));
    % EEMD
    list_cmp{4}.sig  = raw_sig;                list_cmp{4}.name = 'EEMD | Original (mix)';
    list_cmp{5}.sig  = results.EEMD.S0_recon;  list_cmp{5}.name = sprintf('EEMD | S0 recon (IMFs %s)', int2str(results.EEMD.idx_S0));
    list_cmp{6}.sig  = results.EEMD.A0_recon;  list_cmp{6}.name = sprintf('EEMD | A0 recon (IMFs %s)', int2str(results.EEMD.idx_A0));
    % CEEMDAN
    list_cmp{7}.sig  = raw_sig;                   list_cmp{7}.name = 'CEEMDAN | Original (mix)';
    list_cmp{8}.sig  = results.CEEMDAN.S0_recon;   list_cmp{8}.name = sprintf('CEEMDAN | S0 recon (IMFs %s)', int2str(results.CEEMDAN.idx_S0));
    list_cmp{9}.sig  = results.CEEMDAN.A0_recon;   list_cmp{9}.name = sprintf('CEEMDAN | A0 recon (IMFs %s)', int2str(results.CEEMDAN.idx_A0));

    plot_dual_view_AGU_curve(list_cmp, time_us, fs, dispCmp, ...
        'Comparison: Reconstructions (EMD vs EEMD vs CEEMDAN)', ref_max_raw, raw_ylim, curve_data, curveOpt);
end

%% ================= 7) 保存结果（可选） =================
if saveOpt.save_mat
    [~, base, ~] = fileparts(wavefile);
    outmat = fullfile(pname, sprintf('%s_%s_%s.mat', saveOpt.prefix, base, datestr(now,'yyyymmdd_HHMMSS')));
    save(outmat, 'results', 'fs', 'time_us', 'wavefile', 'dataOpt', 'timeOpt', 'curveOpt', 'sepOpt', 'emdStdOpt', 'eemdOpt', 'ceemdanOpt');
    fprintf('\nSaved results to:\n  %s\n', outmat);
end

fprintf('\nDone.\n');

%% =================================================================================
%% ================================== Functions ===================================
%% =================================================================================

function plot_method(tag, R, time_us, fs, dispOpt, ref_max_raw, raw_ylim, curve_data, curveOpt, sepOpt)
    % --- Decomposition figure ---
    imfs = R.imfs;
    K = size(imfs,2);
    [pkHz, cHz] = estimate_mode_freqs(imfs, fs);

    ord_plot = 1:K;
    if sepOpt.sort_plot_by_centroid
        [~, ord_plot] = sort(cHz, 'descend');
    end
    imfs_plot = imfs(:, ord_plot);
    c_plot    = cHz(ord_plot);
    pk_plot   = pkHz(ord_plot);

    list_decomp = cell(1, 1+K);
    list_decomp{1}.sig  = R.raw;
    list_decomp{1}.name = sprintf('%s | Original (mix)', tag);
    for k = 1:K
        orig_id = ord_plot(k);
        list_decomp{k+1}.sig = imfs_plot(:,k);
        list_decomp{k+1}.name = sprintf('%s | IMF %d (centroid %.1f kHz, peak %.1f kHz)', ...
            tag, orig_id, c_plot(k)/1000, pk_plot(k)/1000);
    end

    plot_dual_view_AGU_curve(list_decomp, time_us, fs, dispOpt, ...
        sprintf('%s: Decomposition (IMFs)', tag), ref_max_raw, raw_ylim, curve_data, curveOpt);

    % --- Reconstruction figure ---
    list_recon = cell(1,3);
    list_recon{1}.sig  = R.raw;
    list_recon{1}.name = sprintf('%s | Original (mix)', tag);
    list_recon{2}.sig  = R.S0_recon;
    list_recon{2}.name = sprintf('%s | S0 recon (IMFs %s)', tag, int2str(R.idx_S0));
    list_recon{3}.sig  = R.A0_recon;
    list_recon{3}.name = sprintf('%s | A0 recon (IMFs %s)', tag, int2str(R.idx_A0));

    plot_dual_view_AGU_curve(list_recon, time_us, fs, dispOpt, ...
        sprintf('%s: A0/S0 Reconstruction', tag), ref_max_raw, raw_ylim, curve_data, curveOpt);
end

function print_method_metrics(tag, R, raw_sig, fs)
    % 完备性误差：raw ?= sum(IMFs)+res
    xhat = sum(R.imfs,2) + R.res;
    e_rec = norm(raw_sig - xhat) / (norm(raw_sig) + eps);

    % A0/S0 的低频比（用同一个 lowband ratio 指标做“整体趋势”对比）
    rA0 = lowband_energy_ratio_cwt(R.A0_recon, fs, 150e3, 32);
    rS0 = lowband_energy_ratio_cwt(R.S0_recon, fs, 150e3, 32);

    % A0 与 S0 的谱重叠 SOI（越小越“分离”）
    soi = spectral_overlap_index(R.A0_recon, R.S0_recon, fs);

    if isfield(R,'scoreA0') && ~isempty(R.scoreA0)
        mA = mean(R.scoreA0(R.idx_A0), 'omitnan');
        mS = mean(R.scoreS0(R.idx_S0), 'omitnan');
        fprintf('%-8s | K=%d | idx_A0=%s | idx_S0=%s | e_rec=%.3e | rLow(A0)=%.3f | rLow(S0)=%.3f | SOI=%.3f | meanScore(A0)=%.3f | meanScore(S0)=%.3f\n', ...
            tag, size(R.imfs,2), mat2str(R.idx_A0), mat2str(R.idx_S0), e_rec, rA0, rS0, soi, mA, mS);
    else
        fprintf('%-8s | K=%d | idx_A0=%s | idx_S0=%s | e_rec=%.3e | rLow(A0)=%.3f | rLow(S0)=%.3f | SOI=%.3f\n', ...
            tag, size(R.imfs,2), mat2str(R.idx_A0), mat2str(R.idx_S0), e_rec, rA0, rS0, soi);
    end
end

function soi = spectral_overlap_index(xA, xS, fs)
    % SOI = ∫min(PA,PS)/∫mean(PA,PS)
    xA = xA(:); xS = xS(:);
    nfft = 2^nextpow2(min(4096, numel(xA)));
    [PA, f] = pwelch(xA, hamming(round(nfft/2)), [], nfft, fs);
    [PS, ~] = pwelch(xS, hamming(round(nfft/2)), [], nfft, fs);
    num = trapz(f, min(PA,PS));
    den = trapz(f, 0.5*(PA+PS)) + eps;
    soi = num/den;
end

function R = group_and_recon(raw_sig, imfs, res, fs, sepOpt, dispOpt, time_us, curve_data, curveOpt)
    % 统一输出结构：imfs/res + A0/S0 分组 + 重构
    % 支持两种重构方式：
    %   (1) 传统：整条 IMF -> A0 或 S0
    %   (2) 分时段：同一 IMF 按局部脊线匹配分数拆分到 A0/S0（推荐）
    if isempty(imfs)
        imfs = zeros(numel(raw_sig), 0);
    end
    if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
    K = size(imfs,2);

    lowRatio = zeros(1,K);

    % ridge 评分结构（全局）
    scoreA0 = nan(1,K); scoreS0 = nan(1,K);
    ridgeA0 = nan(1,K); ridgeS0 = nan(1,K);
    curveA0 = nan(1,K); curveS0 = nan(1,K);
    covA0   = nan(1,K); covS0   = nan(1,K);
    nValidCols = zeros(1,K);
    class_by = cell(1,K);

    % 分时段掩膜（每列对应一个 time sample；本脚本 CWT 列数与 time_us 长度一致）
    useSegmented = isfield(sepOpt,'ridge') && isfield(sepOpt.ridge,'segmented_recon') && sepOpt.ridge.segmented_recon;
    A0_masks = zeros(numel(raw_sig), K);
    S0_masks = zeros(numel(raw_sig), K);
    seg_share_A0 = nan(1,K);

    canUseRidge = strcmpi(sepOpt.method, 'ridge') && ...
        isfield(curve_data,'ok') && curve_data.ok && ...
        isfield(curve_data,'A0') && isfield(curve_data,'S0');

    idx_A0 = [];
    idx_S0 = [];

    for k = 1:K
        xk = imfs(:,k);
        lowRatio(k) = lowband_energy_ratio_cwt(xk, fs, sepOpt.split_freq_Hz, dispOpt.voicesPerOctave);

        usedRidge = false;
        globalClass = ''; % 'A0'/'S0'
        M = [];

        if canUseRidge
            try
                M = ridge_match_scores_imf(xk, time_us, fs, dispOpt, curve_data, sepOpt.ridge);

                scoreA0(k) = M.scoreA0;
                scoreS0(k) = M.scoreS0;
                ridgeA0(k) = M.ridgeMatchA0;
                ridgeS0(k) = M.ridgeMatchS0;
                curveA0(k) = M.curveEnergyA0;
                curveS0(k) = M.curveEnergyS0;
                covA0(k)   = M.coverageA0;
                covS0(k)   = M.coverageS0;
                nValidCols(k) = M.nValidCols;

                sA = M.scoreA0; sS = M.scoreS0;
                if M.nValidCols >= sepOpt.ridge.min_valid_cols && isfinite(sA) && isfinite(sS) ...
                        && abs(sA - sS) >= sepOpt.ridge.score_margin
                    if sA > sS
                        idx_A0(end+1) = k; %#ok<AGROW>
                        globalClass = 'A0';
                    else
                        idx_S0(end+1) = k; %#ok<AGROW>
                        globalClass = 'S0';
                    end
                    class_by{k} = 'ridge';
                    usedRidge = true;
                end
            catch ME
                warning('IMF %d ridge matching failed -> fallback lowband. (%s)', k, ME.message);
            end
        end

        if ~usedRidge
            if lowRatio(k) >= sepOpt.ratio_th
                idx_A0(end+1) = k; %#ok<AGROW>
                globalClass = 'A0';
            else
                idx_S0(end+1) = k; %#ok<AGROW>
                globalClass = 'S0';
            end
            class_by{k} = 'lowband';
        end

        % ---- 为该 IMF 构造 A0/S0 时间掩膜（用于分时段重构）----
        if useSegmented && canUseRidge && ~isempty(M) && isfield(M,'localScoreA0') && isfield(M,'localScoreS0')
            [mA, mS] = build_local_time_masks_from_scores(M, sepOpt.ridge, globalClass);
        else
            if strcmp(globalClass,'A0')
                mA = ones(size(xk)); mS = zeros(size(xk));
            else
                mA = zeros(size(xk)); mS = ones(size(xk));
            end
        end

        % 保底长度对齐
        mA = mA(:); mS = mS(:);
        if numel(mA) ~= numel(xk)
            nn = numel(xk);
            mA = interp1(linspace(0,1,numel(mA)), mA, linspace(0,1,nn), 'linear', 'extrap').';
            mA = max(0,min(1,mA));
            mS = 1 - mA;
        end

        A0_masks(:,k) = mA;
        S0_masks(:,k) = mS;
        seg_share_A0(k) = sum(abs(xk).*mA) / (sum(abs(xk))+eps);
    end

    % 保护：避免空组（基于全局分类索引）
    if isempty(idx_A0) && K>0
        if any(isfinite(scoreA0 - scoreS0))
            [~, imax] = max(scoreA0 - scoreS0);
            idx_A0 = imax;
        else
            [~, imax] = max(lowRatio);
            idx_A0 = imax;
        end
    end
    idx_S0 = setdiff(1:K, idx_A0);

    if isempty(idx_S0) && K>0
        if any(isfinite(scoreS0 - scoreA0))
            [~, imaxS] = max(scoreS0 - scoreA0);
            idx_S0 = imaxS;
        else
            [~, imin] = min(lowRatio);
            idx_S0 = imin;
        end
        idx_A0 = setdiff(1:K, idx_S0);
    end

    % ---- 重构：优先使用分时段掩膜；否则回退整条 IMF ----
    A0_recon = zeros(size(raw_sig));
    S0_recon = zeros(size(raw_sig));
    if useSegmented
        for k = 1:K
            xk = imfs(:,k);
            A0_recon = A0_recon + xk .* A0_masks(:,k);
            S0_recon = S0_recon + xk .* S0_masks(:,k);
        end
    else
        if ~isempty(idx_A0), A0_recon = sum(imfs(:,idx_A0),2); end
        if ~isempty(idx_S0), S0_recon = sum(imfs(:,idx_S0),2); end
    end

    if isfield(sepOpt,'ridge') && isfield(sepOpt.ridge,'include_res_in_A0_diag') && sepOpt.ridge.include_res_in_A0_diag
        A0_recon = A0_recon + res(:); % 仅诊断用途
    end

    R = struct();
    R.raw = raw_sig(:);
    R.imfs = imfs;
    R.res  = res(:);

    % 旧指标（保留）
    R.lowRatio = lowRatio;

    % 新：基于脊线/色散曲线的评分信息（便于对比）
    R.class_method = sepOpt.method;
    R.class_by     = {class_by{:}};
    R.scoreA0      = scoreA0;
    R.scoreS0      = scoreS0;
    R.ridgeMatchA0 = ridgeA0;
    R.ridgeMatchS0 = ridgeS0;
    R.curveEnergyA0 = curveA0;
    R.curveEnergyS0 = curveS0;
    R.coverageA0   = covA0;
    R.coverageS0   = covS0;
    R.nValidCols   = nValidCols;

    % 分时段信息
    R.segmented_recon_used = logical(useSegmented);
    R.A0_masks = A0_masks;
    R.S0_masks = S0_masks;
    R.seg_share_A0 = seg_share_A0;

    % 用于标题显示的“主导类别 IMF 列表”（仍保留）
    if useSegmented
        idx_A0_disp = find(seg_share_A0 >= 0.5);
        idx_S0_disp = setdiff(1:K, idx_A0_disp);
        if isempty(idx_A0_disp) && K>0, [~,ix]=max(seg_share_A0); idx_A0_disp=ix; idx_S0_disp=setdiff(1:K,ix); end
        if isempty(idx_S0_disp) && K>0, [~,ix]=min(seg_share_A0); idx_S0_disp=ix; idx_A0_disp=setdiff(1:K,ix); end
        R.idx_A0 = idx_A0_disp;
        R.idx_S0 = idx_S0_disp;
    else
        R.idx_A0 = idx_A0;
        R.idx_S0 = idx_S0;
    end

    R.idx_A0_global = idx_A0;
    R.idx_S0_global = idx_S0;

    R.A0_recon = A0_recon;
    R.S0_recon = S0_recon;
end

function [mA, mS] = build_local_time_masks_from_scores(M, ridgeOpt, globalClass)
    % 将 IMF 的局部 A0/S0 分数转成时间掩膜（0~1）
    sA = M.localScoreA0(:);
    sS = M.localScoreS0(:);

    d = sA - sS; % >0 倾向A0，<0 倾向S0
    valid = isfinite(d);

    mA = nan(size(d));
    if any(valid)
        % 不确定列：分差太小
        local_margin = ridgeOpt.local_score_margin;
        strongA = valid & (d >=  local_margin);
        strongS = valid & (d <= -local_margin);
        unsure  = valid & ~(strongA | strongS);

        % 先给确定列
        mA(strongA) = 1;
        mA(strongS) = 0;

        % 不确定列处理：按全局类补齐 or 最近邻
        if isfield(ridgeOpt,'keep_unknown_to_global') && ridgeOpt.keep_unknown_to_global
            if strcmpi(globalClass,'A0')
                mA(unsure) = 1;
            else
                mA(unsure) = 0;
            end
        end

        % 其余空值（含 valid但未赋值 + invalid）做插值/外推
        idxv = find(isfinite(mA));
        if ~isempty(idxv)
            mA = fillmissing(mA, 'nearest');
            mA = fillmissing(mA, 'linear', 'EndValues','nearest');
        else
            if strcmpi(globalClass,'A0'), mA(:)=1; else, mA(:)=0; end
        end
    else
        % 没有有效局部分数 -> 完全回退全局类别
        if strcmpi(globalClass,'A0'), mA = ones(size(d)); else, mA = zeros(size(d)); end
    end

    % 时间平滑（避免硬切造成振铃）
    win = max(1, round(ridgeOpt.mask_smooth_cols));
    if mod(win,2)==0, win = win + 1; end
    if win > 1 && numel(mA) >= win
        mA = movmean(mA, win, 'omitnan');
    end
    mA = max(0, min(1, mA));

    % 如果局部评分覆盖很少，则向全局类别更靠拢（防误分）
    if isfield(M,'nValidCols') && isfield(M,'nCols') && M.nCols > 0
        cov = M.nValidCols / M.nCols;
        if cov < 0.05
            if strcmpi(globalClass,'A0')
                mA = 0.8*mA + 0.2;
            else
                mA = 0.8*mA;
            end
            mA = max(0,min(1,mA));
        end
    end

    mS = 1 - mA;
end

%% ===================== (1) Standard EMD =====================
%% ===================== (1) Standard EMD =====================
function [imfs, res] = emd_standard(x, emdOpt)
    x = x(:);
    if exist('emd','file') ~= 2
        error('找不到 emd()。请确认安装 Signal Processing Toolbox。');
    end

    if emdOpt.use_default_stop
        try
            [imfs, res] = emd(x, 'MaxNumIMF', emdOpt.max_imf);
        catch
            imfs = emd(x);
            if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
            if size(imfs,2) > emdOpt.max_imf
                imfs = imfs(:,1:emdOpt.max_imf);
            end
            res = x - sum(imfs,2);
        end
    else
        try
            [imfs, res] = emd(x, ...
                'SiftRelativeTolerance', emdOpt.siftRelTol, ...
                'Interpolation',         emdOpt.interp, ...
                'MaxNumIMF',             emdOpt.max_imf, ...
                'MaxEnergyRatio',        emdOpt.maxEnergyRatio, ...
                'SiftMaxIterations',     emdOpt.siftMaxIter, ...
                'MaxNumExtrema',         emdOpt.maxNumExtrema, ...
                'Display',               emdOpt.display);
        catch
            % 兼容旧版本：SiftMaxIteration
            [imfs, res] = emd(x, ...
                'SiftRelativeTolerance', emdOpt.siftRelTol, ...
                'Interpolation',         emdOpt.interp, ...
                'MaxNumIMF',             emdOpt.max_imf, ...
                'MaxEnergyRatio',        emdOpt.maxEnergyRatio, ...
                'SiftMaxIteration',      emdOpt.siftMaxIter, ...
                'MaxNumExtrema',         emdOpt.maxNumExtrema, ...
                'Display',               emdOpt.display);
        end
    end

    if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
    if size(imfs,2) > emdOpt.max_imf
        imfs = imfs(:,1:emdOpt.max_imf);
    end
    if isempty(res)
        res = x - sum(imfs,2);
    end
end

%% ===================== (2) EEMD =====================
function [imfs_avg, res_avg, info] = eemd_decompose(x, eemdOpt, emdOpt)
    x = x(:);
    N = numel(x);
    K = eemdOpt.max_imf;

    if ~isempty(eemdOpt.rng_seed)
        rng(eemdOpt.rng_seed);
    end

    sigma = std(x); if sigma==0, sigma=1; end
    noise_std = eemdOpt.noise_std_ratio * sigma;

    imf_sum   = zeros(N, K);
    imf_count = zeros(1, K);
    res_sum   = zeros(N, 1);
    res_count = 0;

    runner = make_emd_runner(emdOpt, K);
    nBase = eemdOpt.ensemble_size;

    if eemdOpt.paired_noise
        for n = 1:nBase
            noise = noise_std * randn(N,1);
            [imfs_p, res_p] = runner(x + noise);
            [imfs_m, res_m] = runner(x - noise);

            [imfs_p, res_p] = pad_to_K(imfs_p, res_p, N, K);
            [imfs_m, res_m] = pad_to_K(imfs_m, res_m, N, K);

            imfs_nm = 0.5*(imfs_p + imfs_m);
            res_nm  = 0.5*(res_p  + res_m);

            for k = 1:K
                if any(imfs_nm(:,k) ~= 0)
                    imf_sum(:,k) = imf_sum(:,k) + imfs_nm(:,k);
                    imf_count(k) = imf_count(k) + 1;
                end
            end
            res_sum   = res_sum + res_nm;
            res_count = res_count + 1;
        end
    else
        for n = 1:nBase
            xn = x + noise_std * randn(N,1);
            [imfs_n, res_n] = runner(xn);
            [imfs_n, res_n] = pad_to_K(imfs_n, res_n, N, K);

            for k = 1:K
                if any(imfs_n(:,k) ~= 0)
                    imf_sum(:,k) = imf_sum(:,k) + imfs_n(:,k);
                    imf_count(k) = imf_count(k) + 1;
                end
            end
            res_sum   = res_sum + res_n;
            res_count = res_count + 1;
        end
    end

    imfs_avg = zeros(N, K);
    for k = 1:K
        if imf_count(k) > 0
            imfs_avg(:,k) = imf_sum(:,k) ./ imf_count(k);
        end
    end

    if res_count > 0
        res_avg = res_sum ./ res_count;
    else
        res_avg = zeros(N,1);
    end

    info = struct();
    info.imf_count = imf_count;
    info.res_count = res_count;
    info.noise_std = noise_std;
end

%% ===================== (3) CEEMDAN =====================
function [imfs_avg, res_avg, info] = ceemdan_decompose(x, ceemdanOpt, emdOpt)
    x = x(:);
    N = numel(x);
    K = ceemdanOpt.max_imf;
    Ne = ceemdanOpt.ensemble_size;

    if ~isfield(ceemdanOpt,'paired_noise'), ceemdanOpt.paired_noise = true; end
    if ~isfield(ceemdanOpt,'early_stop'),  ceemdanOpt.early_stop  = false; end
    if ~isfield(ceemdanOpt,'min_extrema'), ceemdanOpt.min_extrema = 2; end

    if ~isempty(ceemdanOpt.rng_seed)
        rng(ceemdanOpt.rng_seed);
    end

    runnerK = make_emd_runner(emdOpt, K); % 噪声用，取各阶 IMF
    runner1 = make_emd_runner(emdOpt, 1); % 只取 IMF1

    % 预生成噪声并计算其各阶 IMF（unit-std）
    noise_raw = randn(N, Ne);
    noise_raw = noise_raw ./ (std(noise_raw, 0, 1) + eps);

    noise_modes = zeros(N, K, Ne);
    for i = 1:Ne
        [w_imfs, ~] = runnerK(noise_raw(:,i));
        [w_imfs, ~] = pad_to_K(w_imfs, [], N, K);

        for k = 1:K
            wk = w_imfs(:,k);
            s  = std(wk);
            if s > 0
                noise_modes(:,k,i) = wk ./ s;
            else
                noise_modes(:,k,i) = 0;
            end
        end
    end

    imfs_avg = zeros(N, K);
    res      = x;

    noise_amp_k = zeros(1, K);
    imf_count   = zeros(1, K);
    stop_k      = K;

    for k = 1:K
        sig_r = std(res); if sig_r==0, sig_r=1; end
        noise_amp = ceemdanOpt.noise_std_ratio * sig_r;
        noise_amp_k(k) = noise_amp;

        imf_sum = zeros(N,1);

        for i = 1:Ne
            if k == 1
                nvec = noise_amp * noise_raw(:,i);
            else
                nvec = noise_amp * noise_modes(:,k,i);
            end

            if ceemdanOpt.paired_noise
                imf_p = first_imf_only(res + nvec, runner1, N);
                imf_m = first_imf_only(res - nvec, runner1, N);
                imf1  = 0.5*(imf_p + imf_m);
            else
                imf1  = first_imf_only(res + nvec, runner1, N);
            end

            imf_sum = imf_sum + imf1;
        end

        imfs_avg(:,k) = imf_sum ./ Ne;
        imf_count(k)  = Ne;

        res = res - imfs_avg(:,k);

        if ceemdanOpt.early_stop
            n_ext = count_extrema(res);
            if n_ext < ceemdanOpt.min_extrema
                stop_k = k;
                break;
            end
        end
    end

    if stop_k < K
        imfs_avg(:, stop_k+1:K) = 0;
    end

    res_avg = res;

    info = struct();
    info.ensemble_size   = Ne;
    info.max_imf         = K;
    info.noise_std_ratio = ceemdanOpt.noise_std_ratio;
    info.noise_amp_k     = noise_amp_k;
    info.paired_noise    = ceemdanOpt.paired_noise;
    info.imf_count       = imf_count;
    info.stop_k          = stop_k;
end

function imf1 = first_imf_only(x, runner1, N)
    x = x(:);
    try
        [imfs, ~] = runner1(x);
    catch
        imfs = [];
    end
    if isempty(imfs)
        imf1 = zeros(N,1); return;
    end
    if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
    imf1 = imfs(:,1);
    imf1 = imf1(:);
    if numel(imf1) > N
        imf1 = imf1(1:N);
    elseif numel(imf1) < N
        imf1 = [imf1; zeros(N-numel(imf1),1)];
    end
end

function n_ext = count_extrema(x)
    x = x(:);
    if numel(x) < 3, n_ext = 0; return; end
    dx = diff(x);
    s  = sign(dx);
    s(s==0) = 1;
    ds = diff(s);
    n_ext = sum(ds ~= 0);
end

%% ====== Shared helpers: EMD runner (version compatible) ======
%% ====== Build shared inner-EMD options from standard EMD (fair comparison) ======
function emdInnerOpt = build_emd_inneropt_from_std(emdStdOpt, K)
    % 将标准 EMD 的关键 sifting 参数复制给 EEMD/CEEMDAN 的内部单次 EMD
    % 说明：
    %   - make_emd_runner() 内部会再确保 MaxNumIMF=K
    %   - 同时兼容 MATLAB 某些版本使用 SiftMaxIteration/SiftMaxIterations 的差异
    emdInnerOpt = struct();
    emdInnerOpt.custom_args = { ...
        'SiftRelativeTolerance', emdStdOpt.siftRelTol, ...
        'Interpolation',         emdStdOpt.interp, ...
        'MaxEnergyRatio',        emdStdOpt.maxEnergyRatio, ...
        'SiftMaxIteration',      emdStdOpt.siftMaxIter, ...
        'MaxNumExtrema',         emdStdOpt.maxNumExtrema, ...
        'Display',               emdStdOpt.display, ...
        'MaxNumIMF',             K ...
        };
end

function runner = make_emd_runner(emdOpt, K)
    if exist('emd','file') ~= 2
        error('找不到 emd()。请确认安装 Signal Processing Toolbox。');
    end
    args = emdOpt.custom_args;
    hasMaxNum = any(strcmpi(args(1:2:end), 'MaxNumIMF'));
    if ~hasMaxNum
        args = [args, {'MaxNumIMF', K}];
    end

    xt = zeros(256,1);

    try
        [~,~] = emd(xt, args{:});   % 接收输出 -> 不会弹窗;
        runner = @(x) emd_call(x, args, K);
        return;
    catch
    end

    try
       [~,~] = emd(xt, args{:});   % 接收输出 -> 不会弹窗;
        runner = @(x) emd_call(x, {'MaxNumIMF', K}, K);
        warning('emd(): custom_args 部分不兼容，已降级为仅使用 MaxNumIMF=%d。', K);
        return;
    catch
    end

    runner = @(x) emd_call(x, {}, K);
    warning('emd(): name-value 不兼容，已退化为 emd(x) 并手动截断到 K=%d。', K);
end

function [imfs, res] = emd_call(x, args, K)
    x = x(:);
    try
        [imfs, res] = emd(x, args{:});
    catch
        imfs = emd(x, args{:});
        res  = x - sum(imfs,2);
    end
    if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
    if size(imfs,2) > K
        imfs = imfs(:,1:K);
    end
end

function [imfs, res] = pad_to_K(imfs, res, N, K)
    if isempty(imfs)
        imfs = zeros(N,K);
    else
        if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
        if size(imfs,1) ~= N
            imfs = imfs(1:min(end,N),:);
            if size(imfs,1) < N
                imfs = [imfs; zeros(N-size(imfs,1), size(imfs,2))];
            end
        end
        if size(imfs,2) < K
            imfs = [imfs, zeros(N, K-size(imfs,2))];
        elseif size(imfs,2) > K
            imfs = imfs(:,1:K);
        end
    end

    if isempty(res)
        res = zeros(N,1);
    else
        res = res(:);
        if numel(res) ~= N
            res = res(1:min(end,N));
            if numel(res) < N
                res = [res; zeros(N-numel(res),1)];
            end
        end
    end
end

%% ====== Ridge matching classification (IMF vs theoretical A0/S0 ridges) ======
function M = ridge_match_scores_imf(x, t_us, fs, dispOpt, curve_data, ridgeOpt)
    % IMF 与理论 A0/S0 曲线的匹配评分（全局 + 局部）
    % 输出：
    %   全局分数：scoreA0/scoreS0（用于“整条 IMF”主导类判断）
    %   局部分数：localScoreA0/localScoreS0（用于分时段重构掩膜）
    x = x(:);
    t_us = t_us(:).';

    [img, f_kHz] = cwt_mag_linear(x, fs, dispOpt);
    if isempty(img) || all(img(:)==0)
        M = empty_ridge_score_struct(numel(t_us)); return;
    end

    nCols = size(img,2);

    % dominant ridge：每列取最大幅值行
    [colMax, idxMax] = max(img, [], 1);
    gmax = max(colMax);
    if gmax <= 0
        M = empty_ridge_score_struct(nCols); return;
    end

    ridge_f = f_kHz(idxMax);
    ridge_w = colMax ./ (gmax + eps);       % 0~1 权重
    validCol = ridge_w >= ridgeOpt.min_col_rel;

    % 脊线平滑（仅在有效列数量足够时）
    win = max(1, round(ridgeOpt.smooth_cols));
    if mod(win,2)==0, win = win + 1; end
    if win > 1 && sum(validCol) >= 3
        tmp = ridge_f;
        tmp(~validCol) = nan;
        tmp = fillmissing(tmp, 'linear', 'EndValues', 'nearest');
        tmp = movmedian(tmp, win, 'omitnan');
        ridge_f = tmp;
    end

    [fA, okA] = interp_curve_freq_at_t(curve_data.A0.t_us, curve_data.A0.f_kHz, t_us);
    [fS, okS] = interp_curve_freq_at_t(curve_data.S0.t_us, curve_data.S0.f_kHz, t_us);

    vA = validCol & okA;
    vS = validCol & okS;
    nValid = sum(validCol);

    [ridgeMatchA0, covA0, localRidgeA0] = ridge_distance_score_with_profile(ridge_f, ridge_w, fA, vA, ridgeOpt.freq_tol_kHz);
    [ridgeMatchS0, covS0, localRidgeS0] = ridge_distance_score_with_profile(ridge_f, ridge_w, fS, vS, ridgeOpt.freq_tol_kHz);

    [curveEnergyA0, localCurveA0] = curve_band_energy_score_with_profile(img, f_kHz, t_us, curve_data.A0.t_us, curve_data.A0.f_kHz, ridgeOpt.band_half_kHz);
    [curveEnergyS0, localCurveS0] = curve_band_energy_score_with_profile(img, f_kHz, t_us, curve_data.S0.t_us, curve_data.S0.f_kHz, ridgeOpt.band_half_kHz);

    wr = ridgeOpt.w_ridge;
    if isfield(ridgeOpt,'w_curve'), wc = ridgeOpt.w_curve; else, wc = 1-wr; end
    ssum = wr + wc;
    if ssum <= 0, wr = 0.7; wc = 0.3; ssum = 1; end
    wr = wr/ssum; wc = wc/ssum;

    % 局部分数（用于分时段掩膜）
    localScoreA0 = wr*localRidgeA0 + wc*localCurveA0;
    localScoreS0 = wr*localRidgeS0 + wc*localCurveS0;

    % 全局分数（忽略 NaN）
    scoreA0 = local_profile_robust_mean(localScoreA0, ridge_w);
    scoreS0 = local_profile_robust_mean(localScoreS0, ridge_w);

    % 若局部汇总失败，回退旧式全局分数组合
    if ~isfinite(scoreA0), scoreA0 = wr*ridgeMatchA0 + wc*curveEnergyA0; end
    if ~isfinite(scoreS0), scoreS0 = wr*ridgeMatchS0 + wc*curveEnergyS0; end

    M = struct();
    M.scoreA0 = scoreA0;
    M.scoreS0 = scoreS0;
    M.ridgeMatchA0 = ridgeMatchA0;
    M.ridgeMatchS0 = ridgeMatchS0;
    M.curveEnergyA0 = curveEnergyA0;
    M.curveEnergyS0 = curveEnergyS0;
    M.coverageA0 = covA0;
    M.coverageS0 = covS0;
    M.nValidCols = nValid;
    M.nCols = nCols;

    % 局部信息（分时段重构用）
    M.localScoreA0 = localScoreA0(:);
    M.localScoreS0 = localScoreS0(:);
    M.localRidgeA0 = localRidgeA0(:);
    M.localRidgeS0 = localRidgeS0(:);
    M.localCurveA0 = localCurveA0(:);
    M.localCurveS0 = localCurveS0(:);
    M.localValidCol = validCol(:);
    M.ridgeWeight = ridge_w(:);
end

function M = empty_ridge_score_struct(nCols)
    if nargin < 1, nCols = 0; end
    M = struct('scoreA0',nan,'scoreS0',nan, ...
        'ridgeMatchA0',nan,'ridgeMatchS0',nan, ...
        'curveEnergyA0',nan,'curveEnergyS0',nan, ...
        'coverageA0',nan,'coverageS0',nan,'nValidCols',0,'nCols',nCols, ...
        'localScoreA0',nan(nCols,1),'localScoreS0',nan(nCols,1), ...
        'localRidgeA0',nan(nCols,1),'localRidgeS0',nan(nCols,1), ...
        'localCurveA0',nan(nCols,1),'localCurveS0',nan(nCols,1), ...
        'localValidCol',false(nCols,1),'ridgeWeight',zeros(nCols,1));
end

function [img, f_kHz] = cwt_mag_linear(x, fs, dispOpt)
    [cfs, f] = cwt(x, fs, 'amor', 'VoicesPerOctave', dispOpt.voicesPerOctave);
    img = abs(cfs);
    if isempty(img)
        f_kHz = []; return;
    end
    if f(1) > f(end), f = flipud(f); img = flipud(img); end
    f_lin = linspace(0, dispOpt.fmax_kHz*1000, dispOpt.nFreqBins);
    img = interp1(f, img, f_lin, 'linear', 0);
    f_kHz = f_lin(:)'/1000;
    mx = max(img(:));
    if mx > 0, img = img ./ mx; end
end

function [f_curve, ok] = interp_curve_freq_at_t(t_curve, f_curve_kHz, t_query_us)
    tq = t_query_us(:);
    tc = t_curve(:);
    fc = f_curve_kHz(:);

    idx = isfinite(tc) & isfinite(fc);
    tc = tc(idx); fc = fc(idx);
    if numel(tc) < 2
        f_curve = nan(size(tq)); ok = false(size(tq)); return;
    end

    [tc, ia] = unique(tc, 'stable');
    fc = fc(ia);

    [tc, is] = sort(tc);
    fc = fc(is);

    f_curve = interp1(tc, fc, tq, 'linear', nan);
    ok = isfinite(f_curve);
    f_curve = f_curve(:).';
    ok = ok(:).';
end

function [score, coverage, localScore] = ridge_distance_score_with_profile(ridge_f_kHz, ridge_w, f_curve_kHz, validMask, tol_kHz)
    validMask = validMask(:).';
    ridge_f_kHz = ridge_f_kHz(:).';
    ridge_w = ridge_w(:).';
    f_curve_kHz = f_curve_kHz(:).';

    localScore = nan(size(ridge_f_kHz));
    idx = validMask & isfinite(ridge_f_kHz) & isfinite(f_curve_kHz) & isfinite(ridge_w);

    coverage = mean(validMask);
    if ~any(idx)
        score = nan; return;
    end

    err = abs(ridge_f_kHz(idx) - f_curve_kHz(idx));
    s = exp(-(err./max(tol_kHz, eps)).^2);   % 0~1
    localScore(idx) = s;

    w = max(ridge_w(idx), 0);
    if all(w==0), w = ones(size(w)); end
    score = sum(w .* s) / (sum(w) + eps);
end

function [score, localVals] = curve_band_energy_score_with_profile(img, f_kHz_axis, t_us_axis, t_curve_us, f_curve_kHz, bandHalf_kHz)
    [f_curve_q, ok] = interp_curve_freq_at_t(t_curve_us, f_curve_kHz, t_us_axis);
    nT = numel(t_us_axis);
    localVals = nan(1, nT);

    if ~any(ok)
        score = nan; return;
    end

    for j = 1:nT
        if ~ok(j), continue; end
        fj = f_curve_q(j);
        idxf = abs(f_kHz_axis - fj) <= bandHalf_kHz;
        if ~any(idxf), continue; end
        col = img(idxf, j);
        if isempty(col), continue; end
        localVals(j) = max(col);
    end

    if all(~isfinite(localVals))
        score = nan;
    else
        score = median(localVals(isfinite(localVals))); % 抗孤立峰
    end
end

function score = local_profile_robust_mean(v, w)
    v = v(:); 
    if nargin < 2 || isempty(w), w = ones(size(v)); else, w = w(:); end
    idx = isfinite(v) & isfinite(w);
    if ~any(idx), score = nan; return; end
    vv = v(idx); ww = max(w(idx),0);
    if all(ww==0), ww = ones(size(ww)); end
    % 裁剪极端值后加权均值（更稳）
    lo = quantile(vv, 0.05);
    hi = quantile(vv, 0.95);
    vv = min(max(vv, lo), hi);
    score = sum(vv.*ww)/(sum(ww)+eps);
end

%% ====== Frequency estimation & lowband energy ratio ======

function [pkHz, cHz] = estimate_mode_freqs(modes, fs)
    [N, K] = size(modes);
    L = N;
    f = fs*(0:(L/2))/L;
    pkHz = zeros(1,K);
    cHz  = zeros(1,K);

    for k = 1:K
        x = modes(:,k);
        Y = fft(x);
        P2 = abs(Y/L).^2;
        P1 = P2(1:floor(L/2)+1);

        [~, idx] = max(P1);
        pkHz(k) = f(idx);

        num = sum(f(:).*P1(:));
        den = sum(P1(:)) + eps;
        cHz(k) = num/den;
    end
end

function rLow = lowband_energy_ratio_cwt(x, fs, split_freq_Hz, voicesPerOctave)
    x = x(:);
    [cfs, f] = cwt(x, fs, 'amor', 'VoicesPerOctave', voicesPerOctave);
    mag2 = abs(cfs).^2;
    if f(1) > f(end), f = flipud(f); mag2 = flipud(mag2); end
    idxLow = (f <= split_freq_Hz);
    Elow = sum(mag2(idxLow,:), 'all');
    Etot = sum(mag2, 'all') + eps;
    rLow = Elow / Etot;
end

function [ref_max, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt)
    f_lin = linspace(0, dispOpt.fmax_kHz*1000, dispOpt.nFreqBins);
    [cfs, f] = cwt(raw_sig, fs, 'amor', 'VoicesPerOctave', dispOpt.voicesPerOctave);
    mag = abs(cfs);
    if f(1) > f(end), f = flipud(f); mag = flipud(mag); end
    mag_lin = interp1(f, mag, f_lin, 'linear', 0);

    ref_max = max(mag_lin(:));
    if ref_max == 0, ref_max = 1; end

    a = max(abs(raw_sig));
    if a == 0, a = 1; end
    raw_ylim = [-1.05*a, 1.05*a];
end

%% ====== AGUlike plot (time + CWT) with optional dispersion curves ======
function plot_dual_view_AGU_curve(plotList, t_us, fs, dispOpt, figTitle, ref_max_raw, raw_ylim, curve_data, curveOpt)
    nRows = numel(plotList);
    f_lin = linspace(0, dispOpt.fmax_kHz*1000, dispOpt.nFreqBins);

    stored_cwt = cell(1,nRows);
    all_max = 0;
    for i = 1:nRows
        [cfs, f] = cwt(plotList{i}.sig, fs, 'amor', 'VoicesPerOctave', dispOpt.voicesPerOctave);
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
        % ---- left: time ----
        nexttile;
        plot(t_us, plotList{i}.sig, 'b', 'LineWidth', 1);
        grid on; xlim([0, t_us(end)]);
        ylabel('Amp [mV]');
        title(plotList{i}.name, 'FontWeight','bold', 'FontSize', 10);
        if i == 1 && dispOpt.force_same_rawrow && ~isempty(raw_ylim)
            ylim(raw_ylim);
        end
        if i < nRows, xticklabels([]); else, xlabel('Time [\mus]'); end

        % ---- right: CWT ----
        nexttile;
        img = stored_cwt{i};
        local_max = max(img(:));

        is_S0_like = contains(plotList{i}.name, 'S0');
        is_weak    = (local_max < 0.2*all_max);

        % 归一化基准
        if dispOpt.disable_local_norm
            norm_base = all_max; % 强制全局
        else
            if (is_S0_like || is_weak) && local_max > 0 && i > 1
                norm_base = local_max; % 弱行局部归一化（看得更清）
            else
                norm_base = all_max;
            end
        end

        if i == 1 && dispOpt.force_same_rawrow
            norm_base = ref_max_raw; % 保证跨图 Original 行一致
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
            plot(curve_data.A0.t_us, curve_data.A0.f_kHz, ...
                'LineStyle', curveOpt.A0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
            plot(curve_data.S0.t_us, curve_data.S0.f_kHz, ...
                'LineStyle', curveOpt.S0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
            if i == 1
                add_curve_legend(gca, curveOpt);
            end
            hold off;
        end
    end

    cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'Norm |CWT|';
    title(tl, figTitle, 'FontSize', 12, 'FontWeight','bold');
end

function add_curve_legend(ax, curveOpt)
    if isempty(ax) || ~isvalid(ax), return; end
    axes(ax); %#ok<LAXES>
    hold(ax, 'on');

    xl = xlim(ax); yl = ylim(ax);
    dx = xl(2) - xl(1); dy = yl(2) - yl(1);
    if dx <= 0 || dy <= 0, return; end

    x1 = xl(1) + 0.78*dx;
    x2 = xl(1) + 0.92*dx;
    yA = yl(1) + 0.93*dy;
    yS = yl(1) + 0.86*dy;

    plot(ax, [x1 x2], [yA yA], '-',  'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth, 'Clipping', 'on');
    text(ax, x2 + 0.01*dx, yA, 'A0', 'Color','w', 'FontWeight','bold', 'VerticalAlignment','middle');

    plot(ax, [x1 x2], [yS yS], '--', 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth, 'Clipping', 'on');
    text(ax, x2 + 0.01*dx, yS, 'S0', 'Color','w', 'FontWeight','bold', 'VerticalAlignment','middle');
end

%% ====== Dispersion curve loader (Vallen-like + generic table) ======
function curve = load_dispersion_curve(curvefile, dist_mm, shift_us)
    curve = struct('ok',false,'msg','', 'A0',struct(),'S0',struct());
    if ~isfile(curvefile)
        curve.msg = '文件不存在'; return;
    end

    try
        raw = readcell(curvefile);
        [ok, A0, S0, msg] = parse_vallen_like_cell(raw, dist_mm, shift_us);
        if ok
            curve.ok = true; curve.A0 = A0; curve.S0 = S0; curve.msg = msg; return;
        end
    catch
    end

    try
        T = readtable(curvefile, 'VariableNamingRule','preserve');

        freq_col = find_col(T, {'freq','frequency','f'}, {});
        if isempty(freq_col), curve.msg='未找到频率列（需含 freq/frequency/f）'; return; end
        f_raw = T.(freq_col{1}); f_raw = f_raw(:);
        f_mean = mean(f_raw(~isnan(f_raw)));
        if f_mean < 10,      f_kHz = f_raw*1000;     % MHz->kHz
        elseif f_mean < 1e4, f_kHz = f_raw;          % kHz
        else,               f_kHz = f_raw/1000;      % Hz->kHz
        end

        tA_col = find_col(T, {'time','t'}, {'A0'});
        tS_col = find_col(T, {'time','t'}, {'S0'});
        if ~isempty(tA_col) && ~isempty(tS_col)
            curve.A0.f_kHz = f_kHz; curve.A0.t_us = T.(tA_col{1})(:) + shift_us;
            curve.S0.f_kHz = f_kHz; curve.S0.t_us = T.(tS_col{1})(:) + shift_us;
            curve.ok = true; curve.msg = 'table time-columns parsed'; return;
        end

        vA_col = find_col(T, {'v','vel','velocity'}, {'A0'});
        vS_col = find_col(T, {'v','vel','velocity'}, {'S0'});
        if isempty(vA_col) || isempty(vS_col)
            curve.msg='未找到 A0/S0 速度列（列名需含 A0/S0 且含 v/vel/velocity），或未找到 time 列。';
            return;
        end

        vA = T.(vA_col{1}); vS = T.(vS_col{1});
        vA = vA(:); vS = vS(:);

        vA_mean = mean(vA(~isnan(vA)));
        if vA_mean > 50, vA_mm_per_us = vA/1000; else, vA_mm_per_us = vA; end
        vS_mean = mean(vS(~isnan(vS)));
        if vS_mean > 50, vS_mm_per_us = vS/1000; else, vS_mm_per_us = vS; end

        curve.A0.f_kHz = f_kHz;
        curve.A0.t_us  = dist_mm ./ (vA_mm_per_us + eps) + shift_us;

        curve.S0.f_kHz = f_kHz;
        curve.S0.t_us  = dist_mm ./ (vS_mm_per_us + eps) + shift_us;

        curve.ok = true; curve.msg = 'table velocity-columns parsed';
    catch ME
        curve.msg = sprintf('readtable 解析失败：%s', ME.message);
    end
end

function [ok, A0, S0, msg] = parse_vallen_like_cell(raw, dist_mm, shift_us)
    ok = false; msg = '';
    A0 = struct(); S0 = struct();

    [rS0, cS0] = find_cell(raw, 'Group S0');
    [rA0, cA0] = find_cell(raw, 'Group A0');
    if isempty(rS0) || isempty(rA0)
        msg = '未找到 Group S0 / Group A0'; return;
    end

    r_data = rS0(1) + 2;

    fS = cell2num(raw(r_data:end, cS0(1)));
    vS = cell2num(raw(r_data:end, cS0(1)+1));
    fA = cell2num(raw(r_data:end, cA0(1)));
    vA = cell2num(raw(r_data:end, cA0(1)+1));

    [fS, vS] = trim_nan_pairs(fS, vS);
    [fA, vA] = trim_nan_pairs(fA, vA);

    if numel(fS) < 5 || numel(fA) < 5
        msg = 'S0/A0 数据点太少（表结构可能不标准）'; return;
    end

    S0.f_kHz = fS * 1000;
    A0.f_kHz = fA * 1000;

    % v[m/ms] == mm/us
    S0.t_us  = dist_mm ./ (vS + eps) + shift_us;
    A0.t_us  = dist_mm ./ (vA + eps) + shift_us;

    ok = true;
    msg = 'parsed as Vallen-like dispersion export';
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
    x = nan(size(col,1),1);
    for i = 1:size(col,1)
        v = col{i,1};
        if isnumeric(v) && isfinite(v)
            x(i) = double(v);
        elseif ischar(v) || isstring(v)
            vv = str2double(v);
            if isfinite(vv), x(i) = vv; end
        end
    end
end

function [a,b] = trim_nan_pairs(a,b)
    n = min(numel(a), numel(b));
    a = a(1:n); b = b(1:n);
    idx = isfinite(a) & isfinite(b);
    if ~any(idx), return; end
    last = find(idx, 1, 'last');
    a = a(1:last); b = b(1:last);
end

function colname = find_col(T, mustContainAny, mustContainAll)
    vars = T.Properties.VariableNames;
    vars_low = lower(string(vars));

    any_ok = false(size(vars_low));
    for k = 1:numel(mustContainAny)
        any_ok = any_ok | contains(vars_low, lower(string(mustContainAny{k})));
    end

    all_ok = true(size(vars_low));
    for k = 1:numel(mustContainAll)
        all_ok = all_ok & contains(vars_low, lower(string(mustContainAll{k})));
    end

    idx = find(any_ok & all_ok, 1, 'first');
    if isempty(idx), colname = {}; else, colname = {vars{idx}}; end
end