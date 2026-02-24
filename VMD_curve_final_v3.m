%% ================== VMD_curve.m ==================
% 目的：
%   (1) 对单个 AE burst 信号进行 VMD 分解（优先用 MATLAB vmd；若不可用则尝试开源 VMD.m）；
%   (2) 每个 VMD mode 以“AGUlike”风格显示：左侧时域，右侧 CWT 幅值图；
%   (3) 基于每个 mode 的 CWT 低频能量占比，自动将 mode 分为 A0-like / S0-like 两组并重构；
%   (4) 在所有 CWT 图上叠加 **A0 / S0 色散到达时间曲线**（来自 Vallen 导出的 Excel/CSV 或自定义表格）；
%
% 输出：
%   Figure 1：VMD 分解（Original + 每个 Mode 的时域 + CWT）
%   Figure 2：A0/S0 重构（Original / A0_recon / S0_recon 的时域 + CWT）
%
% 说明：
%   - 所有 “threshold/gamma/归一化” 仅影响显示效果，不改变分解或重构结果。
%   - 曲线叠加支持先选 wave，再选曲线文件；dist_mm 与 shift_us 可在开头快速调参。

clear; clc; close all;

%% ================= 0) 用户可调参数（优先改这里） =================
% -------- (A) 色散曲线叠加：距离与时移 --------
curveOpt.enable   = true;
curveOpt.dist_mm  = 100;      % [mm]
curveOpt.shift_us = 66;       % [us]

% -------- (B) 读 wave 的通道/事件 --------
dataOpt.sel_ch    = 1;
dataOpt.sel_event = 7;

% -------- (C) VMD 参数（影响分解结果的“频带数量/带宽惩罚”）--------
vmdOpt.K     = 3;             % 模态数 K：越大 -> 频带划分更细；过大可能把同一物理模态拆碎
vmdOpt.alpha = 5000;         % 带宽惩罚：越大 -> 每个模态更“窄带”；过大可能造成欠拟合/能量分散
vmdOpt.tol   = 1e-5;          % 收敛阈值：越小 -> 收敛更严格（更慢但更稳定）

% -------- (D) CWT/显示参数：影响“图的清晰度与对比度”--------
dispOpt.fmax_kHz          = 400;
dispOpt.nFreqBins         = 400;
dispOpt.voicesPerOctave   = 32;
dispOpt.threshold         = 0.05;
dispOpt.gamma             = 1.3;
dispOpt.force_same_rawrow = true;

% -------- (E) A0 / S0 自动分组参数（基于 CWT 低频能量占比）--------
sepOpt.split_freq_Hz  = 150e3;
sepOpt.ratio_th       = 0.60;
sepOpt.use_sort_by_centroid = true;  % true: 显示顺序按质心频率排序（不改分解）

% -------- (F) 曲线绘图样式 --------
curveOpt.lineWidth = 1.6;
curveOpt.A0_style  = '-';
curveOpt.S0_style  = '--';
curveOpt.color     = [1 1 1];

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

%% ================= 4) VMD 分解 =================
fprintf('Running VMD (K=%d, alpha=%g)...\n', vmdOpt.K, vmdOpt.alpha);
modes = vmd_compat(raw_sig, vmdOpt.K, vmdOpt.alpha, vmdOpt.tol);

% 确保 modes 为 N x K
if size(modes,1) < size(modes,2), modes = modes.'; end
K = size(modes,2);

% 估计每个 mode 的频率特征（用于标题）
[pkHz, cHz] = estimate_mode_freqs(modes, fs);

% 可选：按质心频率排序【仅用于绘图顺序】
% 说明：
%   - 只调整“绘图显示”的 mode 顺序（高频 -> 低频）
%   - 不直接改写 modes / cHz / pkHz，避免影响后续自动分组、重构和打印的索引一致性
ord_plot = 1:K;
if sepOpt.use_sort_by_centroid
    [~, ord_plot] = sort(cHz, 'descend');
end
modes_plot = modes(:, ord_plot);
c_plot     = cHz(ord_plot);
pk_plot    = pkHz(ord_plot);

fprintf('Mode 质心频率 (kHz)：%s\n', mat2str(round(cHz/1000,1)));

%% ================= 5) 基于 CWT 低频能量占比，自动分组并重构 A0/S0 =================
lowRatio = zeros(1,K);
for k = 1:K
    lowRatio(k) = lowband_energy_ratio_cwt(modes(:,k), fs, sepOpt.split_freq_Hz, dispOpt.voicesPerOctave);
end

idx_A0 = find(lowRatio >= sepOpt.ratio_th);
idx_S0 = setdiff(1:K, idx_A0);

% 【简化保护】阈值导致某组为空时，取 lowRatio 的极值来保证两组非空（不再做“按频率一分为二”的 fallback）。
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

fprintf('VMD 自动分组：split=%.0f kHz, ratio_th=%.2f\n', sepOpt.split_freq_Hz/1000, sepOpt.ratio_th);
fprintf('  A0-like modes = %s\n', mat2str(idx_A0));
fprintf('  S0-like modes = %s\n', mat2str(idx_S0));

A0_recon = sum(modes(:, idx_A0), 2);
S0_recon = sum(modes(:, idx_S0), 2);

%% ================= 6) 构造绘图列表 =================
% Figure 1：Original + modes
list_decomp = cell(1, 1+K);
list_decomp{1}.sig  = raw_sig;
list_decomp{1}.name = 'Original signal (mix)';
for k = 1:K
    list_decomp{k+1}.sig  = modes_plot(:,k);
    orig_id = ord_plot(k);
    list_decomp{k+1}.name = sprintf('Mode %d (centroid %.1f kHz, peak %.1f kHz)', ...
        orig_id, c_plot(k)/1000, pk_plot(k)/1000);
end

% Figure 2：Original + A0 + S0
list_recon = cell(1,3);
list_recon{1}.sig  = raw_sig;
list_recon{1}.name = 'Original signal (mix)';
list_recon{2}.sig  = S0_recon;
list_recon{2}.name = sprintf('S0 reconstruction (modes %s)', int2str(idx_S0));
list_recon{3}.sig  = A0_recon;
list_recon{3}.name = sprintf('A0 reconstruction (modes %s)', int2str(idx_A0));


%% ================= 7) 绘图（AGUlike 风格 + 曲线叠加） =================
[ref_max_raw, raw_ylim] = compute_raw_display_refs(raw_sig, fs, dispOpt);

plot_dual_view_AGU_curve(list_decomp, time_us, fs, dispOpt, ...
    'Figure 1: VMD Decomposition (Modes)', ref_max_raw, raw_ylim, curve_data, curveOpt);

plot_dual_view_AGU_curve(list_recon, time_us, fs, dispOpt, ...
    'Figure 2: A0/S0 Reconstruction (VMD-only)', ref_max_raw, raw_ylim, curve_data, curveOpt);

fprintf('Done.\n');

%% =================================================================================
%% ================================== 函数区 ======================================
%% =================================================================================

function modes = vmd_compat(x, K, alpha, tol)
% vmd_compat：尽量兼容两种 VMD 实现
%   1) MATLAB Signal Processing Toolbox 的 vmd()
%   2) 开源 VMD.m（Dragomiretskiy & Zosso）函数名通常为 VMD()
%
% 输出 modes 统一为 N x K
    x = x(:);

    % --- 1) 优先使用 MATLAB 内置 vmd ---
    try
        [modes, ~, ~] = vmd(x, 'NumIMFs', K, 'PenaltyFactor', alpha, 'AbsoluteTolerance', tol);
        return;
    catch
        % 继续尝试开源 VMD()
    end

    % --- 2) 尝试开源 VMD.m（函数名 VMD）---
    if exist('VMD','file') == 2
        try
            tau = 0; DC = 0; init = 1;
            [u, ~, ~] = VMD(x, alpha, tau, K, DC, init, tol);
            % u 常见输出为 K x N
            if size(u,1) == K && size(u,2) == numel(x)
                modes = u.'; % N x K
            else
                modes = u;
                if size(modes,1) < size(modes,2), modes = modes.'; end
            end
            return;
        catch ME
            error('尝试调用开源 VMD() 失败：%s', ME.message);
        end
    end

    error('无法调用 vmd() 或 VMD()。请确认：1) 安装 Signal Processing Toolbox；或 2) 路径中有开源 VMD.m。');
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

% --------- 曲线读取与绘图函数：与 AGUlike_curve.m 完全一致（复制一份，便于单文件运行） ---------

function curve = load_dispersion_curve(curvefile, dist_mm, shift_us)
    curve = struct('ok',false,'msg','', 'A0',struct(),'S0',struct());
    if ~isfile(curvefile)
        curve.msg = '文件不存在';
        return;
    end

    [~,~,ext] = fileparts(curvefile);
    ext = lower(ext);

    try
        if any(strcmp(ext, {'.xlsx','.xls'}))
            raw = readcell(curvefile);
        else
            raw = readcell(curvefile);
        end

        [ok, A0, S0, msg] = parse_vallen_like_cell(raw, dist_mm, shift_us);
        if ok
            curve.ok = true; curve.A0 = A0; curve.S0 = S0; curve.msg = 'Vallen-like parsed';
            return;
        end
    catch
    end

    try
        T = readtable(curvefile, 'VariableNamingRule','preserve');

        freq_col = find_col(T, {'freq','frequency','f'}, {});
        if isempty(freq_col)
            curve.msg = '未找到频率列（列名需含 freq/frequency/f）';
            return;
        end
        f_raw = T.(freq_col{1}); f_raw = f_raw(:);

        f_mean = mean(f_raw(~isnan(f_raw)));
        if f_mean < 10
            f_kHz = f_raw * 1000;
        elseif f_mean < 1e4
            f_kHz = f_raw;
        else
            f_kHz = f_raw / 1000;
        end

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

        vA_mean = mean(vA(~isnan(vA)));
        if vA_mean > 50
            vA_mm_per_us = vA / 1000;
        else
            vA_mm_per_us = vA;
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
    ok = false; msg = '';
    A0 = struct(); S0 = struct();

    if isempty(raw) || size(raw,1) < 10
        msg = 'raw 太小';
        return;
    end

    [rS0, cS0] = find_cell(raw, 'Group S0');
    [rA0, cA0] = find_cell(raw, 'Group A0');

    if isempty(rS0) || isempty(rA0)
        msg = '未在单元格中找到 Group S0 / Group A0';
        return;
    end

    r_data = rS0(1) + 2;

    fS = cell2num(raw(r_data:end, cS0(1)));
    vS = cell2num(raw(r_data:end, cS0(1)+1));
    fA = cell2num(raw(r_data:end, cA0(1)));
    vA = cell2num(raw(r_data:end, cA0(1)+1));

    [fS, vS] = trim_nan_pairs(fS, vS);
    [fA, vA] = trim_nan_pairs(fA, vA);

    if numel(fS) < 5 || numel(fA) < 5
        msg = 'S0/A0 数据点太少（可能表结构不标准）';
        return;
    end

    S0.f_kHz = fS * 1000;
    A0.f_kHz = fA * 1000;

    vS_mm_per_us = vS;
    vA_mm_per_us = vA;

    S0.t_us = dist_mm ./ (vS_mm_per_us + eps) + shift_us;
    A0.t_us = dist_mm ./ (vA_mm_per_us + eps) + shift_us;

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
        nexttile;
        plot(t_us, plotList{i}.sig, 'b', 'LineWidth', 1);
        grid on; xlim([0, t_us(end)]);
        ylabel('Amp [mV]');
        title(plotList{i}.name, 'FontWeight','bold', 'FontSize', 10);
        if i == 1 && dispOpt.force_same_rawrow && ~isempty(raw_ylim)
            ylim(raw_ylim);
        end
        if i < nRows, xticklabels([]); else, xlabel('Time [\mus]'); end

        nexttile;
        img = stored_cwt{i};
        local_max = max(img(:));

        is_S0_like = contains(plotList{i}.name, 'S0');
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
    t = t(:); f = f(:);
    idx = find(isfinite(t) & isfinite(f), 1, 'last');
    if isempty(idx), return; end
    idx2 = max(1, round(0.8*idx));
    x = t(idx2); y = f(idx2);
    text(x, y, txt, 'Color', curveOpt.color, 'FontWeight','bold', ...
        'BackgroundColor', [0 0 0], 'Margin', 2, 'Clipping','on');
end
