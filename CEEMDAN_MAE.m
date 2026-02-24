%% ================== CEEMDAN_curve.m ==================
% 目的：
%   (1) 对单个 AE burst 信号进行 **CEEMDAN** 分解（Torres et al., 2011：Complete Ensemble EMD with Adaptive Noise）；
%   (2) 将每个 CEEMDAN “IMF” 以 AGUlike 风格显示：左侧时域，右侧 CWT 幅值图；
%   (3) 基于每个 IMF 的 CWT 低频能量占比，将 IMF 分为 A0-like / S0-like 并重构；
%   (4) 在所有 CWT 图上叠加 A0/S0 色散到达时间曲线（Vallen Excel/CSV 或自定义表格）。
%
% 输出：
%   Figure 1：CEEMDAN 分解（Original + 每个 IMF 的时域 + CWT，叠加 A0/S0 曲线）
%   Figure 2：A0/S0 重构（Original / S0_recon / A0_recon 的时域 + CWT，叠加 A0/S0 曲线）
%
% 参考（建议在论文/摘要中引用）：
%   - M. E. Torres, M. A. Colominas, G. Schlotthauer, P. Flandrin, “A Complete Ensemble Empirical Mode Decomposition with Adaptive Noise,” ICASSP 2011.
%   - Z. Wu, N. E. Huang, “Ensemble Empirical Mode Decomposition: A Noise-Assisted Data Analysis Method,” Advances in Adaptive Data Analysis, 2009.
%
% 说明：
%   - CEEMDAN 的关键差异在于：每阶 IMF 的提取都在“当前残差 + 与该尺度对应的噪声 IMF”上进行，
%     并通过集合平均得到唯一 IMF，从而恢复 EMD 的“完备性（completeness）”并减轻 mode mixing。

clear; clc; close all;

%% ================= 0) 用户可调参数（优先改这里） =================
% -------- (A) 只分析前 focus_us 微秒（最小改动窗口）--------
timeOpt.focus_us = 500;     % 关注窗口长度 [us]
timeOpt.start_us = 0;       % 从信号起点开始 [us]

% -------- (B) 色散曲线叠加：距离与时移 --------
curveOpt.enable   = true;   % true: 叠加色散曲线；false: 不叠加
curveOpt.dist_mm  = 200;    % 距离 [mm]
curveOpt.shift_us = 18;     % 曲线整体平移 [us]（对齐触发/起点偏移）

% -------- (C) 读 wave 的通道/事件 --------
dataOpt.sel_ch    = 1;
dataOpt.sel_event = 4;

% -------- (D) CEEMDAN 参数（核心：集合数 & 噪声幅值）--------
ceemdanOpt.max_imf         = 5;      % 输出 IMF 数（越大越“细”，也更容易把低频趋势拆出来）
ceemdanOpt.ensemble_size   = 100;     % 集合次数 N（越大越稳健，但越慢；常用 50~200）
ceemdanOpt.noise_std_ratio = 0.2;   % 每阶噪声幅值 = ratio * std(当前残差)（常用 0.1~0.3）
ceemdanOpt.rng_seed        = 1;      % 固定随机种子，保证可复现（想“每次不同”就设 []）
ceemdanOpt.paired_noise    = true;   % true：每次用 +noise 和 -noise 成对平均（更稳、更少残余噪声）
                                 % 说明：这是常见的“互补噪声（paired noise）”做法，能明显降低集合平均后的残余噪声偏置

% -------- (E) 单次 EMD 的参数（尽量按你以前那套写法尝试；版本不支持会自动降级）--------
emdOpt.custom_args = { ...
    'SiftRelativeTolerance', 0.1, ...   % 你之前用的
    'Interpolation', 'pchip', ...       % 你之前用的
    'MaxEnergyRatio', 20, ...           % 你之前用的
    'SiftMaxIteration', 100, ...        % 你之前用的（有的版本叫 MaxNumSift）
    'MaxNumExtrema', 1, ...             % 你之前用的
    'Display', 0 ...                    % 关显示（避免刷屏）
    };
% 注意：不同 MATLAB 版本 emd() 支持的 name-value 不同：
%   - 若某些参数不被支持，本代码会“自动去掉不支持的项”并降级运行
%   - 你不会再因为版本差异而直接报错停住

% -------- (F) CWT/显示参数：影响“图的清晰度与对比度”（只影响显示）--------
dispOpt.fmax_kHz          = 400;
dispOpt.nFreqBins         = 400;
dispOpt.voicesPerOctave   = 32;
dispOpt.threshold         = 0.05;
dispOpt.gamma             = 1.5;
dispOpt.force_same_rawrow = true;   % Figure1/2 的 Original 行共享同一显示归一化

% -------- (G) A0 / S0 自动分组参数（基于 CWT 低频能量占比）--------
sepOpt.split_freq_Hz  = 150e3; % 低/高频分界（<= 认为低频）
sepOpt.ratio_th       = 0.55;  % lowRatio >= ratio_th -> A0-like，否则 S0-like
sepOpt.sort_plot_by_centroid = true; % 仅影响 Figure1 的 IMF 显示顺序（高频→低频）

% -------- (H) 曲线绘图样式 --------
curveOpt.lineWidth = 1.6;
curveOpt.A0_style  = '-';      % A0 实线
curveOpt.S0_style  = '--';     % S0 虚线
curveOpt.color     = [1 1 1];  % 白色

%% ================= 1) 选择 wave 文件（信号） =================
fprintf('请选择数据文件（*.wave / *.tradb）...\n');
[fname, pname] = uigetfile({'*.wave;*.tradb;*.*'}, '选择 wave 文件');
if isequal(fname, 0), error('未选择 wave 文件'); end
wavefile = fullfile(pname, fname);

%% ================= 2) 读 wave -> 原始信号 raw_sig =================
[t0, wave, ~, ~, ~, ~, sf, ~, ~, ~, ~, ~, ~] = waveReader(wavefile);
fs = double(sf);

try
    if ndims(wave) == 3
        raw_sig = double(wave(:, dataOpt.sel_ch, dataOpt.sel_event));
    else
        raw_sig = double(wave(:, dataOpt.sel_ch));
    end
catch
    raw_sig = double(wave(:, 1));
end

raw_sig(~isfinite(raw_sig)) = 0;
raw_sig = raw_sig(:) - mean(raw_sig);

% ===== 只关注前 focus_us 微秒（你要的“开头最小修改”本质就是这段）=====
start_idx = max(1, round(timeOpt.start_us*1e-6*fs) + 1);
end_idx   = min(length(raw_sig), start_idx + round(timeOpt.focus_us*1e-6*fs) - 1);
raw_sig   = raw_sig(start_idx:end_idx);
N         = numel(raw_sig);
time_us   = (0:N-1)/fs * 1e6;
% =====================================================================

%% ================= 3) 选择并读取色散曲线文件（可选） =================
curve_data = struct();
if curveOpt.enable
    fprintf('请选择色散曲线文件（Vallen 导出 Excel/CSV，或自定义表格）...\n');
    [cname, cpname] = uigetfile({'*.xlsx;*.xls;*.csv;*.txt;*.*'}, '选择色散曲线文件（可取消）');
    if isequal(cname, 0)
        warning('未选择曲线文件：将不叠加色散曲线。');
        curveOpt.enable = false;
    else
        curvefile = fullfile(cpname, cname);
        curve_data = load_dispersion_curve(curvefile, curveOpt.dist_mm, curveOpt.shift_us);
        if ~curve_data.ok
            warning('曲线文件读取失败：%s\n将不叠加曲线。', curve_data.msg);
            curveOpt.enable = false;
        end
    end
end

%% ================= 4) CEEMDAN 分解 =================
fprintf('Running CEEMDAN: max_imf=%d, N=%d, noise_ratio=%.3f, paired=%d ...\n', ...
    ceemdanOpt.max_imf, ceemdanOpt.ensemble_size, ceemdanOpt.noise_std_ratio, ceemdanOpt.paired_noise);

[imfs, res, ceemdanInfo] = ceemdan_decompose(raw_sig, fs, ceemdanOpt, emdOpt);

% imfs: N x K
K = size(imfs,2);

% 估计每个 IMF 的频率特征（用于标题）
[pkHz, cHz] = estimate_mode_freqs(imfs, fs);
fprintf('CEEMDAN IMF 质心频率 (kHz)：%s\n', mat2str(round(cHz/1000,1)));

%% ================= 5) 基于 CWT 低频能量占比，自动分组并重构 A0/S0 =================
lowRatio = zeros(1,K);
for k = 1:K
    lowRatio(k) = lowband_energy_ratio_cwt(imfs(:,k), fs, sepOpt.split_freq_Hz, dispOpt.voicesPerOctave);
end

idx_A0 = find(lowRatio >= sepOpt.ratio_th);
idx_S0 = setdiff(1:K, idx_A0);

% 简化保护：避免空组（不做“按频率一分为二”的复杂 fallback）
if isempty(idx_A0)
    [~, imax] = max(lowRatio);
    idx_A0 = imax;
end
idx_S0 = setdiff(1:K, idx_A0);
if isempty(idx_S0)
    [~, imin] = min(lowRatio);
    idx_S0 = imin;
    idx_A0 = setdiff(1:K, idx_S0);
end

fprintf('CEEMDAN 自动分组：split=%.0f kHz, ratio_th=%.2f\n', sepOpt.split_freq_Hz/1000, sepOpt.ratio_th);
fprintf('  A0-like IMFs = %s\n', mat2str(idx_A0));
fprintf('  S0-like IMFs = %s\n', mat2str(idx_S0));

A0_recon = sum(imfs(:, idx_A0), 2);
S0_recon = sum(imfs(:, idx_S0), 2);

%% ================= 6) 构造绘图列表 =================
% Figure 1：Original + IMFs（可选：按质心频率排序显示，但不改变 imfs 本身）
ord_plot = 1:K;
if sepOpt.sort_plot_by_centroid
    [~, ord_plot] = sort(cHz, 'descend'); % 高频 -> 低频
end
imfs_plot = imfs(:, ord_plot);
c_plot    = cHz(ord_plot);
pk_plot   = pkHz(ord_plot);

list_decomp = cell(1, 1+K);
list_decomp{1}.sig  = raw_sig;
list_decomp{1}.name = 'Original signal (mix)';
for k = 1:K
    orig_id = ord_plot(k);
    list_decomp{k+1}.sig  = imfs_plot(:,k);
    list_decomp{k+1}.name = sprintf('IMF %d (centroid %.1f kHz, peak %.1f kHz)', ...
        orig_id, c_plot(k)/1000, pk_plot(k)/1000);
end

% Figure 2：Original + S0 + A0（保持你现在图里的顺序：Original / S0 / A0）
list_recon = cell(1,3);
list_recon{1}.sig  = raw_sig;
list_recon{1}.name = 'Original signal (mix)';
list_recon{2}.sig  = S0_recon;
list_recon{2}.name = sprintf('S0 reconstruction (IMFs %s)', int2str(idx_S0));
list_recon{3}.sig  = A0_recon;
list_recon{3}.name = sprintf('A0 reconstruction (IMFs %s)', int2str(idx_A0));

%% ================= 7) 绘图（AGUlike 风格 + 曲线叠加） =================
[ref_max_raw, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt);

plot_dual_view_AGU_curve(list_decomp, time_us, fs, dispOpt, ...
    'Figure 1: CEEMDAN Decomposition (IMFs)', ref_max_raw, raw_ylim, curve_data, curveOpt);

plot_dual_view_AGU_curve(list_recon, time_us, fs, dispOpt, ...
    'Figure 2: A0/S0 Reconstruction (CEEMDAN)', ref_max_raw, raw_ylim, curve_data, curveOpt);

fprintf('Done.\n');

%% =================================================================================
%% ================================== 函数区 ======================================
%% =================================================================================


function [imfs_avg, res_avg, info] = ceemdan_decompose(x, fs, ceemdanOpt, emdOpt)
% ceemdan_decompose：CEEMDAN 核心（Complete Ensemble EMD with Adaptive Noise）
%
% 输入：
%   x            : 原始信号（列向量）
%   fs           : 采样率（此处不参与分解，但保留以与其它方法脚本的输入保持一致）
%   ceemdanOpt   : 参数结构体
%       .max_imf         : 输出 IMF 数 K
%       .ensemble_size   : 集合次数 Ne
%       .noise_std_ratio : 噪声幅值比例（每阶噪声幅值 = ratio * std(当前残差)）
%       .rng_seed        : 随机种子（[] 表示不固定）
%       .paired_noise    : 是否用互补噪声（+noise 与 -noise）做平均（可明显减小残余噪声偏置）
%       .early_stop      : true: 当残差极值数很少时提前停止（剩余 IMF 置零）
%       .min_extrema     : early_stop 的极值阈值（默认 2）
%   emdOpt       : 单次 EMD 的参数（custom_args 等，复用你其它脚本里的版本兼容逻辑）
%
% 输出：
%   imfs_avg : N x K 的 CEEMDAN IMF（每阶 IMF 唯一、集合平均）
%   res_avg  : N x 1 残差（x - sum(IMF,2)，保证“完备性”数值上成立）
%   info     : 记录噪声幅值、集合次数、停止信息等（便于你写摘要/补充材料）
%
% 算法要点（Torres et al., 2011）：
%   - 第 1 阶：对 x + 噪声 做 EMD，只取 IMF1；集合平均得到 IMF1。
%   - 第 k 阶(k>=2)：对残差 r_{k-1} + 噪声的第 k 阶 IMF(white noise) 做 EMD，只取 IMF1；
%                   集合平均得到 IMF_k；更新残差 r_k = r_{k-1} - IMF_k。
%
% 备注：
%   - 这里采用“每阶噪声幅值自适应于当前残差标准差”的常见实现，
%     在 AE burst（幅值快速衰减）情况下更稳健。
%   - 若你要严格复现某个开源实现的噪声标定方式，可在 noise_amp_k 的计算处替换。

    %#ok<NASGU>
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

    % --- 构造 EMD 调用器 ---
    runnerK = make_emd_runner(emdOpt, K);  % 给噪声做 EMD（取各阶 IMF）
    runner1 = make_emd_runner(emdOpt, 1);  % 只取 IMF1（提速）

    % --- 预生成 Ne 组白噪声，并预计算其各阶 IMF（用于自适应噪声注入） ---
    noise_raw = randn(N, Ne);                 % N x Ne
    % 列归一化到 unit-std（保证噪声幅值标定可控）
    noise_raw = noise_raw ./ (std(noise_raw, 0, 1) + eps);

    noise_modes = zeros(N, K, Ne);            % 噪声的各阶 IMF（每阶 unit-std）
    for i = 1:Ne
        [w_imfs, ~] = runnerK(noise_raw(:,i));
        [w_imfs, ~] = pad_to_K(w_imfs, [], N, K);

        for k = 1:K
            wk = w_imfs(:,k);
            s  = std(wk);
            if s > 0
                noise_modes(:,k,i) = wk ./ s; % unit-std
            else
                noise_modes(:,k,i) = 0;
            end
        end
    end

    % --- CEEMDAN 主循环 ---
    imfs_avg = zeros(N, K);
    res      = x;

    noise_amp_k = zeros(1, K);
    imf_count   = zeros(1, K);
    stop_k      = K;  % 实际停止阶

    for k = 1:K
        % 每阶噪声幅值：与当前残差 std 成比例（自适应）
        sig_r = std(res);
        if sig_r == 0, sig_r = 1; end
        noise_amp = ceemdanOpt.noise_std_ratio * sig_r;
        noise_amp_k(k) = noise_amp;

        imf_sum = zeros(N,1);

        for i = 1:Ne
            if k == 1
                nvec = noise_amp * noise_raw(:,i);           % 第1阶：直接加白噪声
            else
                nvec = noise_amp * noise_modes(:,k,i);       % 第k阶：加噪声的第k阶 IMF
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

        % 更新残差（保证数值上完备：x = sum(IMF) + res）
        res = res - imfs_avg(:,k);

        % 可选提前停止：残差极值太少 -> 近似趋势项
        if ceemdanOpt.early_stop
            n_ext = count_extrema(res);
            if n_ext < ceemdanOpt.min_extrema
                stop_k = k;
                break;
            end
        end
    end

    % 若提前停止，后续 IMF 保持为 0（与其它脚本“输出固定 K 阶”的习惯一致）
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
% first_imf_only：对输入 x 做一次 EMD，仅取 IMF1，保证输出为 N x 1
    x = x(:);
    try
        [imfs, ~] = runner1(x);
    catch
        imfs = [];
    end

    if isempty(imfs)
        imf1 = zeros(N,1);
        return;
    end
    if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
    imf1 = imfs(:,1);
    imf1 = imf1(:);

    % 长度保护（极少数实现可能发生）
    if numel(imf1) > N
        imf1 = imf1(1:N);
    elseif numel(imf1) < N
        imf1 = [imf1; zeros(N-numel(imf1),1)];
    end
end

function n_ext = count_extrema(x)
% count_extrema：估计极值数（通过一阶差分符号变化计数）
    x = x(:);
    if numel(x) < 3
        n_ext = 0; return;
    end
    dx = diff(x);
    s  = sign(dx);
    s(s==0) = 1;          % 将 0 当作 +1，避免出现大量 0 导致误判
    ds = diff(s);
    n_ext = sum(ds ~= 0); % 符号变化次数 ~ 极值数量
end

function runner = make_emd_runner(emdOpt, K)
% make_emd_runner：构造一个“尽量兼容版本差异”的 emd 调用器
% 策略：
%   1) 优先尝试 emd(x, custom_args{:}, 'MaxNumIMF',K)（若你 custom_args 里没写 MaxNumIMF）
%   2) 若失败，逐级降级：只保留 MaxNumIMF
%   3) 再失败：直接 emd(x) 然后截断到 K
    if exist('emd','file') ~= 2
        error('当前 MATLAB 环境找不到 emd()。请确认安装 Signal Processing Toolbox 或将 EMD 实现加入路径。');
    end

    % 预处理 custom_args：如果你没显式给 MaxNumIMF，就补上
    args = emdOpt.custom_args;
    hasMaxNum = any(strcmpi(args(1:2:end), 'MaxNumIMF'));
    if ~hasMaxNum
        args = [args, {'MaxNumIMF', K}];
    end

    % 试一次探测（用短向量即可）
    xt = zeros(256,1);

    % 1) 全量 custom args
    try
        emd(xt, args{:});
        runner = @(x) emd_call(x, args, K);
        return;
    catch
    end

    % 2) 只用 MaxNumIMF
    try
        emd(xt, 'MaxNumIMF', K);
        runner = @(x) emd_call(x, {'MaxNumIMF', K}, K);
        warning('emd(): 你的版本不支持 custom_args 中的部分参数，已自动降级为仅使用 MaxNumIMF=%d。', K);
        return;
    catch
    end

    % 3) 最后退化：emd(x) 然后手动截断
    runner = @(x) emd_call(x, {}, K);
    warning('emd(): 当前版本不支持常用 name-value，已退化为 emd(x) 并手动截断到 K=%d。', K);
end

function [imfs, res] = emd_call(x, args, K)
% emd_call：真正执行一次 EMD，并尽量兼容输出个数差异
    x = x(:);
    try
        % 新版本可能输出 [imf,residual,info]
        [imfs, res] = emd(x, args{:});
    catch
        % 老版本可能只输出 imf
        imfs = emd(x, args{:});
        res  = x - sum(imfs,2);
    end

    % 手动截断到 K（避免输出过多）
    if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
    if size(imfs,2) > K
        imfs = imfs(:,1:K);
    end
end

function [imfs, res] = pad_to_K(imfs, res, N, K)
% pad_to_K：将一次分解的 imfs/res 对齐到 N x K，缺的补零
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

function [ref_max, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt)
    f_lin = linspace(0, dispOpt.fmax_kHz*1000, dispOpt.nFreqBins);
    [cfs, f] = cwt(raw_sig, fs, 'amor', 'VoicesPerOctave', dispOpt.voicesPerOctave);
    mag = abs(cfs);
    if f(1) > f(end), f = flipud(f); mag = flipud(mag); end
    mag_lin = interp1(f, mag, f_lin, 'linear', 0);
    ref_max = max(mag_lin(:));
    if ref_max == 0, ref_max = 1; end
    a = max(abs(raw_sig)); if a == 0, a = 1; end
    raw_ylim = [-1.05*a, 1.05*a];
end

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

%% ================= 曲线读取与绘图：与你现有脚本一致（复制版，单文件可运行） =================

function curve = load_dispersion_curve(curvefile, dist_mm, shift_us)
    curve = struct('ok',false,'msg','', 'A0',struct(),'S0',struct());
    if ~isfile(curvefile)
        curve.msg = '文件不存在'; return;
    end

    [~,~,ext] = fileparts(curvefile);
    ext = lower(ext);

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
        if f_mean < 10,      f_kHz = f_raw*1000;
        elseif f_mean < 1e4, f_kHz = f_raw;
        else,               f_kHz = f_raw/1000;
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
    ok = false; msg=''; A0=struct(); S0=struct();
    [rS0, cS0] = find_cell(raw, 'Group S0');
    [rA0, cA0] = find_cell(raw, 'Group A0');
    if isempty(rS0) || isempty(rA0)
        msg='未找到 Group S0 / Group A0'; return;
    end
    r_data = rS0(1) + 2;
    fS = cell2num(raw(r_data:end, cS0(1)));
    vS = cell2num(raw(r_data:end, cS0(1)+1));
    fA = cell2num(raw(r_data:end, cA0(1)));
    vA = cell2num(raw(r_data:end, cA0(1)+1));
    [fS,vS]=trim_nan_pairs(fS,vS);
    [fA,vA]=trim_nan_pairs(fA,vA);
    if numel(fS)<5 || numel(fA)<5
        msg='S0/A0 数据点太少（表结构可能不标准）'; return;
    end
    S0.f_kHz = fS*1000; A0.f_kHz = fA*1000;
    S0.t_us  = dist_mm ./ (vS + eps) + shift_us; % v[m/ms]=mm/us
    A0.t_us  = dist_mm ./ (vA + eps) + shift_us;
    ok=true; msg='parsed as Vallen-like dispersion export';
end

function [r, c] = find_cell(raw, key)
    r=[]; c=[]; key=lower(string(key));
    for i=1:size(raw,1)
        for j=1:size(raw,2)
            v=raw{i,j};
            if ischar(v) || isstring(v)
                s=lower(string(v));
                if contains(s,key)
                    r(end+1)=i; c(end+1)=j; %#ok<AGROW>
                end
            end
        end
    end
end

function x = cell2num(col)
    x = nan(size(col,1),1);
    for i=1:size(col,1)
        v=col{i,1};
        if isnumeric(v) && isfinite(v)
            x(i)=double(v);
        elseif ischar(v) || isstring(v)
            vv=str2double(v);
            if isfinite(vv), x(i)=vv; end
        end
    end
end

function [a,b] = trim_nan_pairs(a,b)
    n=min(numel(a),numel(b));
    a=a(1:n); b=b(1:n);
    idx=isfinite(a)&isfinite(b);
    if ~any(idx), return; end
    last=find(idx,1,'last');
    a=a(1:last); b=b(1:last);
end

function colname = find_col(T, mustContainAny, mustContainAll)
    vars = T.Properties.VariableNames;
    vars_low = lower(string(vars));
    any_ok = false(size(vars_low));
    for k=1:numel(mustContainAny)
        any_ok = any_ok | contains(vars_low, lower(string(mustContainAny{k})));
    end
    all_ok = true(size(vars_low));
    for k=1:numel(mustContainAll)
        all_ok = all_ok & contains(vars_low, lower(string(mustContainAll{k})));
    end
    idx=find(any_ok & all_ok, 1, 'first');
    if isempty(idx), colname={}; else, colname={vars{idx}}; end
end

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
        % --- 左：时域 ---
        nexttile;
        plot(t_us, plotList{i}.sig, 'b', 'LineWidth', 1);
        grid on; xlim([0, t_us(end)]);
        ylabel('Amp [mV]');
        title(plotList{i}.name, 'FontWeight','bold', 'FontSize', 10);
        if i == 1 && dispOpt.force_same_rawrow && ~isempty(raw_ylim)
            ylim(raw_ylim);
        end
        if i < nRows, xticklabels([]); else, xlabel('Time [\mus]'); end

        % --- 右：CWT ---
        nexttile;
        img = stored_cwt{i};
        local_max = max(img(:));

        is_S0_like = contains(plotList{i}.name, 'S0');
        is_weak = (local_max < 0.2*all_max);

        if (is_S0_like || is_weak) && local_max > 0 && i > 1
            norm_base = local_max;  % 弱行局部归一化
        else
            norm_base = all_max;    % 全局归一化
        end
        if i == 1 && dispOpt.force_same_rawrow
            norm_base = ref_max_raw; % 原始行跨图一致
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

        % --- 叠加色散曲线 ---
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
    dx = xl(2)-xl(1); dy = yl(2)-yl(1);
    x1 = xl(1)+0.78*dx; x2 = xl(1)+0.92*dx;
    yA = yl(1)+0.93*dy; yS = yl(1)+0.86*dy;
    plot(ax, [x1 x2],[yA yA], '-',  'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth, 'Clipping','on');
    text(ax, x2+0.01*dx, yA, 'A0', 'Color','w', 'FontWeight','bold', 'VerticalAlignment','middle');
    plot(ax, [x1 x2],[yS yS], '--', 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth, 'Clipping','on');
    text(ax, x2+0.01*dx, yS, 'S0', 'Color','w', 'FontWeight','bold', 'VerticalAlignment','middle');
end