%% ================== AGUlike_curve.m ==================
% 目的：
%   (1) 对单个 AE burst 信号进行 **标准 EMD** 分解（MATLAB emd）；
%   (2) 将每个 IMF 以“AGUlike”风格显示：左侧时域，右侧 CWT 幅值图；
%   (3) 基于每个 IMF 的 CWT 低频能量占比，自动将 IMF 分为 A0-like / S0-like 两组并重构；
%   (4) 在所有 CWT 图上叠加 **A0 / S0 色散到达时间曲线**（来自 Vallen 导出的 Excel/CSV 或自定义表格）；
%
% 输出：
%   Figure 1：EMD 分解（Original + 每个 IMF 的时域 + CWT，叠加 A0/S0 曲线）
%   Figure 2：A0/S0 重构（Original / A0_recon / S0_recon 的时域 + CWT，叠加 A0/S0 曲线）
%
% 你可以在“0) 用户可调参数”里修改距离 dist_mm 和曲线时移 shift_us，
% 以及 EMD/显示/CWT/自动分组阈值等参数（均有中文注释说明影响）。

clear; clc; close all;

%% ================= 0) 用户可调参数（优先改这里） =================
% -------- (A) 色散曲线叠加：距离与时移 --------
curveOpt.enable   = true;    % true: 叠加色散曲线；false: 不叠加
curveOpt.dist_mm  = 100;     % 传播距离 [mm] —— 直接影响曲线的到达时间 t_us = dist_mm / v(m/ms)
curveOpt.shift_us = 66;      % 曲线整体水平平移 [us] —— 用于对齐你当前信号的起始/触发偏移

% -------- (B) 读 wave 的通道/事件 --------
dataOpt.sel_ch    = 1;       % 选择通道索引
dataOpt.sel_event = 1;       % 若 wave 为 3D（N x ch x event），选择事件索引

% -------- (C) EMD 分解参数（标准 EMD 的停止/输出控制）--------
emdOpt.max_imf  = 5;          % 对应 MaxNumIMF
emdOpt.siftRelTol = 0.1;      % 对应 SiftRelativeTolerance
emdOpt.interp  = 'PCHIP';     % 对应 Interpolation: 'spline'/'pchip'
emdOpt.maxEnergyRatio = 20;   % 对应 MaxEnergyRatio
emdOpt.siftMaxIter = 100;     % 对应 SiftMaxIterations (或 SiftMaxIteration)
emdOpt.maxNumExtrema = 1;     % 对应 MaxNumExtrema
emdOpt.display = 1;           % 建议 0；你要 emd 自带显示就设 1
emdOpt.use_default_stop = false; % true: 使用 MATLAB emd 默认停止准则（推荐，等价“标准 EMD”）
                                % false: 你可以在下面 emdOpt.stop_* 中指定停止条件（更可控，但不是“默认标准”）

% -------- (D) CWT/显示参数：影响“图的清晰度与对比度”--------
dispOpt.fmax_kHz          = 400;  % CWT 显示最高频率 [kHz]（越大 -> 频率范围更广，但低频细节相对变小）
dispOpt.nFreqBins         = 400;  % 频率插值点数（越大 -> 频率轴更平滑，但计算更慢）
dispOpt.voicesPerOctave   = 32;   % 小波每倍频程声音数（越大 -> 频率分辨率更高，但时域分辨率更差 & 更慢）
dispOpt.threshold         = 0.01; % 归一化后幅值阈值（越大 -> 背景更干净，但弱成分更容易被抹掉）
dispOpt.gamma             = 1;  % γ增强（>1 增强强信号对比；越大 -> 强分量更突出，弱分量更暗）
dispOpt.force_same_rawrow = false; % true: 强制 Figure1/2 的 Original 行（时域+TF）完全一致（便于论文对比）

% -------- (E) A0 / S0 自动分组参数（基于 CWT 低频能量占比）--------
sepOpt.split_freq_Hz  = 150e3; % 低/高频分界 [Hz]：<= 该频率认为是“低频”
                               %   增大 -> 更多能量被算作低频 -> 更容易判为 A0-like
sepOpt.ratio_th       = 0.50;  % A0 判别阈值：lowRatio >= ratio_th -> A0-like，否则 S0-like
                               %   增大 -> 更“严格”判 A0（A0 组变少，S0 组变多）
sepOpt.use_sort_by_centroid = true; % true: 按质心频率从高到低排序 IMF（仅影响显示顺序，不改分解结果）

% -------- (F) 曲线绘图样式 --------
curveOpt.lineWidth = 1.6;     % 曲线线宽
curveOpt.A0_style  = '-';     % A0 曲线线型
curveOpt.S0_style  = '--';    % S0 曲线线型
curveOpt.color     = [1 1 1]; % 曲线颜色（白色，便于叠加在 CWT 上）

%% ================= 1) 选择 wave 文件（信号） =================
fprintf('请选择数据文件（*.wave / *.tradb）...\n');
[fname, pname] = uigetfile({'*.wave;*.tradb;*.*'}, '选择 wave 文件');
if isequal(fname, 0), error('未选择 wave 文件'); end
wavefile = fullfile(pname, fname);

%% ================= 2) 读 wave -> 原始信号 raw_sig =================
% waveReader 为你已有的函数（MiSTRAS/Hexagon 导出 wave）
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
% 
% N = numel(raw_sig);
% time_us = (0:N-1)/fs * 1e6;
% ================== 只关注前 500 us（最小改动） ==================
focus_us  = 500;          % 关注窗口长度（us）
start_us  = 0;            % 从信号起点开始（us）。如果想从某个时刻开始，改这里

start_idx = max(1, round(start_us*1e-6*fs) + 1);
end_idx   = min(length(raw_sig), start_idx + round(focus_us*1e-6*fs) - 1);

raw_sig   = raw_sig(start_idx:end_idx);   % 直接截断信号（后续分解/CWT/重构都只算这段）
N         = numel(raw_sig);
time_us   = (0:N-1)/fs * 1e6;             % 时间轴重置为 0~focus_us
% ================================================================

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

%% ================= 4) EMD 分解（标准 EMD） =================
% 说明：这里使用 MATLAB 的 emd()，其基本思想与 Huang et al. (1998) 提出的 EMD 一致。
% 若 use_default_stop=true，则采用 MATLAB 默认停止准则（最接近“标准 EMD 的常用实现”）。
fprintf('Running EMD (max_imf=%d)...\n', emdOpt.max_imf);

[imfs, res] = emd_compat(raw_sig, emdOpt);

% 保证 imfs 维度为 N x K
if size(imfs,1) < size(imfs,2), imfs = imfs.'; end
K = size(imfs,2);

% 估计每个 IMF 的频率特征（用于标题）
[pkHz, cHz] = estimate_mode_freqs(imfs, fs);

% 可选：按质心频率排序【仅用于绘图顺序】
% 说明：
%   - EMD 输出的 IMF 编号有其原始意义（例如 IMF1/IMF2…）
%   - 这里如果直接重排 imfs，会导致“IMF 编号”在图里被重定义，后续重构/讨论容易混淆
%   - 因此我们只生成用于绘图的 imfs_plot / c_plot / pk_plot，不修改 imfs 本身
ord_plot = 1:K;
% if sepOpt.use_sort_by_centroid
%     [~, ord_plot] = sort(cHz, 'descend'); % 高频 -> 低频
% end
imfs_plot = imfs(:, ord_plot);
c_plot    = cHz(ord_plot);
pk_plot   = pkHz(ord_plot);

fprintf('IMF 质心频率 (kHz)：%s\n', mat2str(round(cHz/1000,1)));

%% ================= 5) 基于 CWT 低频能量占比，自动分组并重构 A0/S0 =================
lowRatio = zeros(1,K);
for k = 1:K
    lowRatio(k) = lowband_energy_ratio_cwt(imfs(:,k), fs, sepOpt.split_freq_Hz, dispOpt.voicesPerOctave);
end

idx_A0 = find(lowRatio >= sepOpt.ratio_th);
idx_S0 = setdiff(1:K, idx_A0);

% 【简化保护】若阈值导致某一组为空：
%   - A0 组空：取 lowRatio 最大的 IMF 作为 A0
%   - S0 组空：取 lowRatio 最小的 IMF 作为 S0
% 这样既避免报错/空重构，又不引入额外“按频率一分为二”的复杂逻辑。
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

fprintf('EMD 自动分组：split=%.0f kHz, ratio_th=%.2f\n', sepOpt.split_freq_Hz/1000, sepOpt.ratio_th);
fprintf('  A0-like IMFs = %s\n', mat2str(idx_A0));
fprintf('  S0-like IMFs = %s\n', mat2str(idx_S0));

A0_recon = sum(imfs(:, idx_A0), 2);
S0_recon = sum(imfs(:, idx_S0), 2);

%% ================= 6) 构造绘图列表 =================
% Figure 1：Original + IMFs
list_decomp = cell(1, 1+K);
list_decomp{1}.sig  = raw_sig;
list_decomp{1}.name = 'Original signal (mix)';
for k = 1:K
    list_decomp{k+1}.sig  = imfs_plot(:,k);
    orig_id = ord_plot(k);
    list_decomp{k+1}.name = sprintf('IMF %d (centroid %.1f kHz, peak %.1f kHz)', ...
        orig_id, c_plot(k)/1000, pk_plot(k)/1000);
end

% Figure 2：Original + A0 + S0
list_recon = cell(1,3);
list_recon{1}.sig  = raw_sig;
list_recon{1}.name = 'Original signal (mix)';
list_recon{2}.sig  = S0_recon;
list_recon{2}.name = sprintf('S0 reconstruction (IMFs %s)', int2str(idx_S0));
list_recon{3}.sig  = A0_recon;
list_recon{3}.name = sprintf('A0 reconstruction (IMFs %s)', int2str(idx_A0));

%% ================= 7) 绘图（AGUlike 风格 + 曲线叠加） =================
% 为了让 Figure1/2 第一行（Original）完全一致：共享 CWT 最大值和时域 y 轴范围（仅显示层面，不改计算）
[ref_max_raw, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt);

plot_dual_view_AGU_curve(list_decomp, time_us, fs, dispOpt, ...
    'Figure 1: EMD Decomposition (IMFs)', ref_max_raw, raw_ylim, curve_data, curveOpt);

plot_dual_view_AGU_curve(list_recon, time_us, fs, dispOpt, ...
    'Figure 2: A0/S0 Reconstruction', ref_max_raw, raw_ylim, curve_data, curveOpt);

fprintf('Done.\n');

%% =================================================================================
%% ================================== 函数区 ======================================
%% =================================================================================

function [imfs, res] = emd_compat(x, emdOpt)
% emd_compat：兼容不同 MATLAB 版本 emd() 的参数写法
%   x      : 输入信号（列向量）
%   emdOpt : 结构体，见主脚本 0) 参数区
% 输出：
%   imfs : IMF 矩阵（N x K）
%   res  : 残差项（N x 1）
    x = x(:);
    if exist('emd','file') ~= 2
        error('当前 MATLAB 环境找不到 emd()。请确认安装 Signal Processing Toolbox，或将 EMD 实现加入路径。');
    end

    if emdOpt.use_default_stop
        % “标准/默认” EMD：只控制输出数量（最常用、最接近默认实现）
        try
            [imfs, res] = emd(x, 'MaxNumIMF', emdOpt.max_imf);
        catch
            % 某些版本参数名可能不同：尝试不带 name-value
            [imfs, res] = emd(x);
            if size(imfs,2) > emdOpt.max_imf
                imfs = imfs(:,1:emdOpt.max_imf);
            end
        end
    else
        % 自定义停止准则：更可控，但会偏离“默认停止”
        % --- 自定义参数调用（优先用官方名字 SiftMaxIterations）---
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
    % --- 兼容少数版本/第三方实现用的别名：SiftMaxIteration（少一个 s）---
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
end

function [ref_max, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt)
% compute_raw_display_refs：
%   计算“Original 行”统一显示所需的参考量：
%   (1) ref_max：raw_sig 的 CWT 经过插值后的全局最大值（用于统一 TF 归一化）
%   (2) raw_ylim：raw_sig 时域的统一 y 轴范围（用于统一时域幅值显示）
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

function [pkHz, cHz] = estimate_mode_freqs(modes, fs)
% estimate_mode_freqs：
%   用 FFT 能量谱估计每个 mode/IMF 的：
%   - peak frequency（峰值频率）
%   - centroid frequency（频谱质心）
% 返回单位：Hz
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
% lowband_energy_ratio_cwt：
%   在 CWT 时频域计算低频能量占比：
%       rLow = Elow / Etot
%   - split_freq_Hz：低频/高频分界（Hz）
%   - rLow 越大：越偏低频 -> 更可能 A0-like
    x = x(:);
    [cfs, f] = cwt(x, fs, 'amor', 'VoicesPerOctave', voicesPerOctave);
    mag2 = abs(cfs).^2;
    if f(1) > f(end), f = flipud(f); mag2 = flipud(mag2); end
    idxLow = (f <= split_freq_Hz);
    Elow = sum(mag2(idxLow,:), 'all');
    Etot = sum(mag2, 'all') + eps;
    rLow = Elow / Etot;
end

function curve = load_dispersion_curve(curvefile, dist_mm, shift_us)
% load_dispersion_curve：
%   读取色散曲线文件，并转换为“到达时间 t_us vs 频率 f_kHz”的两条曲线（A0/S0）。
%
% 支持两类输入：
%   (A) Vallen Dispersion 导出表（Excel/CSV）：
%       表头包含 "Group S0" / "Group A0"，其下两列为 F[MHz] 与 v[m/ms]
%   (B) 自定义表格（readtable 可读）：
%       - 需要包含频率列（freq/frequency/f 等）
%       - 以及 A0/S0 的速度列（列名含 A0 或 S0，且含 v/vel/velocity）
%       或者直接包含 A0/S0 的到达时间列（列名含 time/t 且含 A0/S0）
%
% 物理换算：
%   Vallen 的 v 单位通常为 m/ms，而 1 m/ms = 1 mm/us，
%   因此到达时间：t_us = dist_mm / v(m/ms) + shift_us
%
% 输出：
%   curve.ok = true/false
%   curve.A0.t_us, curve.A0.f_kHz
%   curve.S0.t_us, curve.S0.f_kHz

    curve = struct('ok',false,'msg','', 'A0',struct(),'S0',struct());
    if ~isfile(curvefile)
        curve.msg = '文件不存在';
        return;
    end

    [~,~,ext] = fileparts(curvefile);
    ext = lower(ext);

    try
        % 优先尝试 Vallen 导出格式（按单元格扫描 "Group S0/A0"）
        if any(strcmp(ext, {'.xlsx','.xls'}))
            raw = readcell(curvefile);
        else
            % csv/txt 也用 readcell（MATLAB R2020+）
            raw = readcell(curvefile);
        end

        [ok, A0, S0, msg] = parse_vallen_like_cell(raw, dist_mm, shift_us);
        if ok
            curve.ok = true; curve.A0 = A0; curve.S0 = S0; curve.msg = 'Vallen-like parsed';
            return;
        end
    catch
        % 忽略，进入 readtable 流程
    end

    % 退化：readtable + 列名匹配
    try
        T = readtable(curvefile, 'VariableNamingRule','preserve');

        % --- 1) 找频率列 ---
        freq_col = find_col(T, {'freq','frequency','f'}, {});
        if isempty(freq_col)
            curve.msg = '未找到频率列（列名需含 freq/frequency/f）';
            return;
        end
        f_raw = T.(freq_col{1});
        f_raw = f_raw(:);

        % 单位猜测：如果均值 < 10，可能是 MHz；如果均值 < 1e4 可能是 kHz；否则可能是 Hz
        f_mean = mean(f_raw(~isnan(f_raw)));
        if f_mean < 10
            f_kHz = f_raw * 1000;      % MHz -> kHz
        elseif f_mean < 1e4
            f_kHz = f_raw;             % kHz
        else
            f_kHz = f_raw / 1000;      % Hz -> kHz
        end

        % --- 2) A0/S0：优先找 time 列，否则找 velocity 列 ---
        tA_col = find_col(T, {'time','t'}, {'A0'});
        tS_col = find_col(T, {'time','t'}, {'S0'});

        if ~isempty(tA_col) && ~isempty(tS_col)
            tA = T.(tA_col{1}); tS = T.(tS_col{1});
            curve.A0.f_kHz = f_kHz; curve.A0.t_us = tA(:) + shift_us;
            curve.S0.f_kHz = f_kHz; curve.S0.t_us = tS(:) + shift_us;
            curve.ok = true; curve.msg = 'table time-columns parsed';
            return;
        end

        vA_col = find_col(T, {'v','vel','velocity'}, {'A0'});
        vS_col = find_col(T, {'v','vel','velocity'}, {'S0'});
        if isempty(vA_col) || isempty(vS_col)
            curve.msg = '未找到 A0/S0 速度列（列名需含 A0/S0 且含 v/vel/velocity），或未找到 A0/S0 time 列。';
            return;
        end

        vA = T.(vA_col{1}); vS = T.(vS_col{1});
        vA = vA(:); vS = vS(:);

        % v 单位猜测：Vallen 常为 m/ms（=mm/us），若值大约 0.1~10 合理；
        % 若值在 100~5000 可能是 m/s，需要换算 mm/us = v_mps / 1e3
        vA_mean = mean(vA(~isnan(vA)));
        if vA_mean > 50
            vA_mm_per_us = vA / 1000;  % m/s -> mm/us
        else
            vA_mm_per_us = vA;         % m/ms -> mm/us
        end

        vS_mean = mean(vS(~isnan(vS)));
        if vS_mean > 50
            vS_mm_per_us = vS / 1000;
        else
            vS_mm_per_us = vS;
        end

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
% parse_vallen_like_cell：
%   从 readcell 得到的 raw 单元格中，定位 "Group S0"/"Group A0" 区块并读取两列数据：
%     F [MHz] | v [m/ms]
%   输出：
%     A0.f_kHz, A0.t_us
%     S0.f_kHz, S0.t_us
    ok = false; msg = '';
    A0 = struct(); S0 = struct();

    if isempty(raw) || size(raw,1) < 10
        msg = 'raw 太小';
        return;
    end

    % 搜索包含 "Group S0" / "Group A0" 的单元格位置
    [rS0, cS0] = find_cell(raw, 'Group S0');
    [rA0, cA0] = find_cell(raw, 'Group A0');

    if isempty(rS0) || isempty(rA0)
        msg = '未在单元格中找到 Group S0 / Group A0';
        return;
    end

    % 频率行通常在 label 下一行（例如 row+1），数据从再下一行开始
    r_head = rS0(1) + 1;
    r_data = rS0(1) + 2;

    % S0 两列：freq = col, vel = col+1
    fS = cell2num(raw(r_data:end, cS0(1)));
    vS = cell2num(raw(r_data:end, cS0(1)+1));

    % A0 两列
    fA = cell2num(raw(r_data:end, cA0(1)));
    vA = cell2num(raw(r_data:end, cA0(1)+1));

    % 清理 NaN 尾部
    [fS, vS] = trim_nan_pairs(fS, vS);
    [fA, vA] = trim_nan_pairs(fA, vA);

    if numel(fS) < 5 || numel(fA) < 5
        msg = 'S0/A0 数据点太少（可能表结构不标准）';
        return;
    end

    % 单位：F 通常为 MHz；v 通常为 m/ms (=mm/us)
    S0.f_kHz = fS * 1000;
    A0.f_kHz = fA * 1000;

    vS_mm_per_us = vS;  % m/ms == mm/us
    vA_mm_per_us = vA;

    S0.t_us = dist_mm ./ (vS_mm_per_us + eps) + shift_us;
    A0.t_us = dist_mm ./ (vA_mm_per_us + eps) + shift_us;

    ok = true;
    msg = 'parsed as Vallen-like dispersion export';
end

function [r, c] = find_cell(raw, key)
% 在单元格 raw 中查找包含 key 的位置（不区分大小写）
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
% 将单元格列转换为数值列（无法转换的置 NaN）
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
% 删除末尾连续 NaN（或非有限）对
    n = min(numel(a), numel(b));
    a = a(1:n); b = b(1:n);
    idx = isfinite(a) & isfinite(b);
    if ~any(idx), return; end
    last = find(idx, 1, 'last');
    a = a(1:last); b = b(1:last);
end

function colname = find_col(T, mustContainAny, mustContainAll)
% 在表 T 的列名中查找：
%   - mustContainAny：列名需包含其中任意一个关键词（不区分大小写）
%   - mustContainAll：列名还需同时包含这些关键词（不区分大小写），例如 {'A0'}
% 返回：第一个匹配到的列名（cell 形式），若无则 {}
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

function plot_dual_view_AGU_curve(plotList, t_us, fs, dispOpt, figTitle, ref_max_raw, raw_ylim, curve_data, curveOpt)
% plot_dual_view_AGU_curve：
%   “AGUlike”显示：每一行一个信号（左：时域；右：CWT），并可叠加 A0/S0 曲线。
%
% 关键显示策略（与你之前 AGUlike 逻辑一致）：
%   1) 对每张图，先计算所有行 CWT 的全局最大值 all_max；
%   2) 对“弱信号行”或名称含 S0 的行，用 local_max 做局部归一化（更容易看清弱成分）；
%   3) 对第 1 行（Original），如果 force_same_rawrow=true，则强制使用 ref_max_raw（保证两张图完全一致）；
%   4) threshold + gamma 用于压背景、增强对比（只影响显示，不影响计算）。
    nRows = numel(plotList);
    f_lin = linspace(0, dispOpt.fmax_kHz*1000, dispOpt.nFreqBins);

    % 1) 预计算每行 CWT（插值到同一频率轴），并统计全局最大值
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

    % 2) 建图
    fig_h = min(1100, 170*nRows);
    figure('Name', figTitle, 'Color','w', 'Position', [60, 60, 1250, fig_h]);
    tl = tiledlayout(nRows, 2, 'TileSpacing','compact', 'Padding','compact');

    % 曲线频率（用于与 f_lin 对应的 kHz 轴）
    f_kHz_axis = f_lin/1000;

    for i = 1:nRows
        % ========== 左：时域 ==========
        nexttile;
        plot(t_us, plotList{i}.sig, 'b', 'LineWidth', 1);
        grid on; xlim([0, t_us(end)]);
        ylabel('Amp [mV]');
        title(plotList{i}.name, 'FontWeight','bold', 'FontSize', 10);
        if i == 1 && dispOpt.force_same_rawrow && ~isempty(raw_ylim)
            ylim(raw_ylim);
        end
        if i < nRows, xticklabels([]); else, xlabel('Time [\mus]'); end

        % ========== 右：CWT ==========
        nexttile;
        img = stored_cwt{i};
        local_max = max(img(:));

        % --- 归一化基准选择（AGUlike 逻辑） ---
        is_S0_like = contains(plotList{i}.name, 'S0');
        is_weak = (local_max < 0.2*all_max);

        if (is_S0_like || is_weak) && local_max > 0 && i > 1
            norm_base = local_max;   % 局部归一化：让弱行更清晰
        else
            norm_base = all_max;     % 全局归一化：保持同图可比性
        end

        % --- Original 行强制一致（跨 Figure1/2） ---
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

        % --- 叠加色散曲线（A0/S0） ---
        if curveOpt.enable && isfield(curve_data,'ok') && curve_data.ok
            hold on;
            hA = plot(curve_data.A0.t_us, curve_data.A0.f_kHz, ...
                'LineStyle', curveOpt.A0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);
            hS = plot(curve_data.S0.t_us, curve_data.S0.f_kHz, ...
                'LineStyle', curveOpt.S0_style, 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth);

            % 线型图例：A0 实线 / S0 虚线（只在第一行显示一次，避免重复/拥挤）
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
% 在当前 TF 轴右上角画“线型图例”，避免依赖曲线末端位置（末端可能超出 0~fmax 导致标签不显示）
%
% 需求：用户希望明确区分 A0 / S0 两条曲线的线型（实线/虚线），且固定在右上角。
% 实现：在坐标轴范围内画两段短线，并在其右侧写上 'A0'/'S0'。
%
% 输入:
%   ax       : 目标坐标轴句柄
%   curveOpt : 曲线绘图选项（包含 lineWidth / color 等）
%
% 注意：这里不使用 legend()，避免与 tiledlayout/colormap 交互导致位置不稳定。

    if isempty(ax) || ~isvalid(ax), return; end
    axes(ax); %#ok<LAXES>
    hold(ax, 'on');

    xl = xlim(ax); yl = ylim(ax);
    dx = xl(2) - xl(1); dy = yl(2) - yl(1);
    if dx <= 0 || dy <= 0, return; end

    % 右上角的“安全内边距”
    x1 = xl(1) + 0.78*dx;
    x2 = xl(1) + 0.92*dx;
    yA = yl(1) + 0.93*dy;
    yS = yl(1) + 0.86*dy;

    % A0: 实线
    plot(ax, [x1 x2], [yA yA], '-',  'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth, 'Clipping', 'on');
    text(ax, x2 + 0.01*dx, yA, 'A0', 'Color','w', 'FontWeight','bold', 'VerticalAlignment','middle');

    % S0: 虚线
    plot(ax, [x1 x2], [yS yS], '--', 'Color', curveOpt.color, 'LineWidth', curveOpt.lineWidth, 'Clipping', 'on');
    text(ax, x2 + 0.01*dx, yS, 'S0', 'Color','w', 'FontWeight','bold', 'VerticalAlignment','middle');

end

function place_label(t, f, txt, curveOpt)
% 在曲线的某个“合适点”放置文本，带黑色背景增强可读性
    t = t(:); f = f(:);
    idx = find(isfinite(t) & isfinite(f), 1, 'last');
    if isempty(idx), return; end
    % 取 80% 位置点（如果太少就取最后一个）
    idx2 = max(1, round(0.8*idx));
    x = t(idx2); y = f(idx2);
    text(x, y, txt, 'Color', curveOpt.color, 'FontWeight','bold', ...
        'BackgroundColor', [0 0 0], 'Margin', 2, 'Clipping','on');
end
