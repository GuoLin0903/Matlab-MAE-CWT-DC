%% ===================== AE_AIC_CWT_basic_auto_multiEvent.m =====================
% 功能概述（精简版）：
%  1) 读取 .wave / .tradb（兼容不同 waveReader 输出）
%  2) 支持 1D / 2D / 3D 波形；3D时可选通道、可多选事件
%  3) 每个事件绘制 Raw + CWT（EMD_MAE 风格）
%  4) 自动找两个波包，计算 Δt（AIC可开关；关闭时用包络峰时刻）
%  5) CWT 可叠加频散曲线（Vallen 导出的 Excel）
%  6) 频散曲线样式：A0白色实线、S0白色虚线；图例白字、透明背景、无边框、右上角
%
% 依赖：
%   - waveReader
%   - Signal Processing Toolbox (butter / filtfilt / hilbert)
%   - Wavelet Toolbox (cwt)

clear; clc; close all;

%% =============== 0) 参数区（主要改这里） ===============

% --- 文件 ---
wavefile = "";                 % 留空=弹窗选择
% wavefile = "D:\Lin\xxx.wave";

% --- 截取范围（可调）---
start_us = 0;                  % 截取起点（us）
sigLen_us = 500;               % 截取长度（us）；设为 inf 表示到信号末尾

% --- 包络/AIC 的带通 ---
fband = [80e3 250e3];          % Hz

% --- AIC 搜索窗（围绕峰值）---
pre_us  = 40;
post_us = 15;

% --- 自动找两波包（稳健参数）---
det.noise_us      = 70;                 % 前段噪声估计长度（us）
det.thrK_list     = [8 6 5 4 3 2.5];    % 阈值逐步放宽
det.mergeGap_us   = 8;                  % 分段合并 gap（us）
det.minDur_us     = 10;                 % 最短段长（us）
det.minSep_us     = 60;                 % 两包最小间隔（us）
det.pkt2_minDelay_us  = 80;             % 第二包最早相对延迟（us）
det.pkt2_maxDelay_us  = 260;            % 第二包最晚相对延迟（us）
det.pkt2_minPeakRatio = 0.18;           % 第二包峰值至少为第一包峰值的比例
det.verbose = false;                    % true = 输出更多调试信息

% --- CWT显示 ---
dispOpt.fmax_kHz        = 500;
dispOpt.nFreqBins       = 700;
dispOpt.voicesPerOctave = 48;
dispOpt.threshold       = 0.05;
dispOpt.gamma           = 1.3;

% --- 绘图控制 ---
plotOpt.doPlot      = true;
plotOpt.maxShow     = inf;              % 只画前N个事件；inf=全画
plotOpt.showMarkers = true;             % AIC开时是否显示拾取虚线

% --- AIC 开关 ---
pickOpt.useAIC = false;                 % true=AIC onset；false=包络峰时刻
pickOpt.showMethodInTitle = false;      % CWT标题显示 [AIC]/[Peak]

% --- 频散曲线叠加（CWT）---
curveOpt.enable       = true;
curveOpt.curveFile    = "";             % 留空=弹窗选择 Excel
% curveOpt.curveFile = "D:\Lin\3mmPMMA_Dispersion Curve.xlsx";

curveOpt.type         = 'group';        % 'group' 或 'phase'
curveOpt.modes        = {'A0','S0'};    % 叠加模态
curveOpt.distance_mm  = 100;            % 源-传感器距离（mm）
curveOpt.timeShift_us = 40;             % 水平平移补偿（us）
curveOpt.lineWidth    = 2.0;            % 频散曲线线宽（建议>=2）
curveOpt.showLegend   = true;           % 白字图例（透明、无框、右上角）
curveOpt.verbose      = false;

%% =============== 1) 选文件 ===============
if strlength(wavefile) == 0
    [fname, pname] = uigetfile({'*.wave;*.tradb;*.*'}, '选择 wave 文件');
    if isequal(fname,0), error('未选择文件'); end
    wavefile = fullfile(pname, fname);
end

%% =============== 2) 读 wave ===============
[wave, fs] = readWaveRobust(wavefile);
wave = double(wave);
fs   = double(fs);

%% =============== 3) 选择通道 + 多事件 ===============
[chSel, evList, sigList] = pickMultiEvents_from_wave(wave);
nEv = numel(evList);

%% =============== 3.1) 读取频散曲线（可选） ===============
curveDB = [];
if curveOpt.enable
    if strlength(curveOpt.curveFile) == 0
        [cfname, cpname] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx,*.xls)'}, ...
                                     '选择频散曲线Excel（Vallen导出）');
        if isequal(cfname,0)
            warning('未选择频散曲线文件，已关闭频散曲线叠加。');
            curveOpt.enable = false;
        else
            curveOpt.curveFile = fullfile(cpname, cfname);
        end
    end

    if curveOpt.enable
        curveDB = readDispersionExcelVallen(curveOpt.curveFile, curveOpt.verbose);
        if isempty(curveDB)
            warning('频散曲线文件读取失败或未识别到有效曲线，已关闭频散曲线叠加。');
            curveOpt.enable = false;
        end
    end
end

%% =============== 4) 结果容器 ===============
% 主输出（最终采用的时刻：AIC或Peak）
t1_us_all = nan(nEv,1);
t2_us_all = nan(nEv,1);
dt_us_all = nan(nEv,1);

% 诊断/辅助（始终保存峰值时刻）
t1_peak_us_all = nan(nEv,1);
t2_peak_us_all = nan(nEv,1);
dt_peak_us_all = nan(nEv,1);

% 分段质量指标
pkt2_segmented = false(nEv,1);
peak2_ratio    = nan(nEv,1);

% AIC关闭时不显示拾取虚线
showPickLines = plotOpt.showMarkers && pickOpt.useAIC;

%% =============== 5) 主循环：逐事件处理 ===============
for i = 1:nEv
    evSel = evList(i);
    sig   = sigList{i};
    sig(~isfinite(sig)) = 0;
    sig = sig(:) - mean(sig(:));

    % ---------- 截断（长度可调） ----------
    start_idx = max(1, round(start_us*1e-6*fs) + 1);
    if isinf(sigLen_us)
        end_idx = numel(sig);
    else
        end_idx = min(numel(sig), start_idx + round(sigLen_us*1e-6*fs) - 1);
    end

    if end_idx <= start_idx
        warning('event=%d: 截取范围无效（start_us/sigLen_us导致长度<=0），跳过。', evSel);
        continue;
    end

    x = sig(start_idx:end_idx);
    N = numel(x);
    t_us = (0:N-1)/fs * 1e6;

    % ---------- 带通 + 包络 ----------
    Wn = fband/(fs/2);
    Wn(1) = max(Wn(1), 1e-6);
    Wn(2) = min(Wn(2), 0.999999);

    [b,a] = butter(4, Wn, 'bandpass');
    x_bp  = filtfilt(b,a,x);

    env   = abs(hilbert(x_bp));
    env_s = movmean(env, max(3, round(2e-6*fs)));

    % ---------- 自动找两个波包 ----------
    [seg1, k1, seg2, k2, pkt2Ok, p2ratio] = autoFindTwoPackets_anchored(env_s, fs, det);
    pkt2_segmented(i,1) = pkt2Ok;
    peak2_ratio(i,1)    = p2ratio;

    % ---------- 峰值时刻（始终计算） ----------
    t1_peak_us = (k1 - 1)/fs * 1e6;
    t2_peak_us = (k2 - 1)/fs * 1e6;

    t1_peak_us_all(i,1) = t1_peak_us;
    t2_peak_us_all(i,1) = t2_peak_us;
    dt_peak_us_all(i,1) = t2_peak_us - t1_peak_us;

    % ---------- AIC 或 Peak（开关） ----------
    if pickOpt.useAIC
        preN  = round(pre_us  *1e-6*fs);
        postN = round(post_us *1e-6*fs);

        t1_on_s = aicOnsetTime_np_limitedRise(x_bp, k1, preN, postN, fs, seg1(1), seg1(2));
        t2_on_s = aicOnsetTime_np_limitedRise(x_bp, k2, preN, postN, fs, seg2(1), seg2(2));

        t1_on_us = t1_on_s * 1e6;
        t2_on_us = t2_on_s * 1e6;
    else
        % AIC关闭：直接用包络峰值时刻
        t1_on_us = t1_peak_us;
        t2_on_us = t2_peak_us;
    end

    t1_us_all(i,1) = t1_on_us;
    t2_us_all(i,1) = t2_on_us;
    dt_us_all(i,1) = t2_on_us - t1_on_us;

    % ---------- 绘图 ----------
    if plotOpt.doPlot && i <= plotOpt.maxShow
        [~, f_kHz_axis, img_show] = cwtImage_EMDstyle(x, fs, dispOpt);

        figure('Color','w','Position',[60 60 1600 720]);
        set(gcf,'Renderer','opengl');  % 复制到PPT更稳（imagesc+line+legend组合）

        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        % ===== Raw =====
        nexttile;
        plot(t_us, x, 'b-', 'LineWidth', 1); grid on;
        xlim([0 t_us(end)]);
        ylabel('Amplitude');
        title(sprintf('Raw (ch=%d, event=%d)  truncated %.0f–%.0f \\mus', ...
            chSel, evSel, start_us, start_us + t_us(end)), 'FontWeight','bold');

        if showPickLines
            xline(t1_on_us,'r--','LineWidth',1);
            xline(t2_on_us,'r--','LineWidth',1);
        end

        % ===== CWT =====
        nexttile;
        imagesc(t_us, f_kHz_axis, img_show);
        axis xy;
        colormap(jet(256)); caxis([0 1]);
        xlim([0 t_us(end)]); ylim([0 dispOpt.fmax_kHz]);
        xlabel('Time [\mus]');
        ylabel('Freq [kHz]');

        if pickOpt.showMethodInTitle
            methodStr = iff(pickOpt.useAIC, 'AIC', 'Peak');
            title(sprintf('CWT   \\Delta t = %.3f \\mus   [%s]', dt_us_all(i), methodStr), ...
                'FontWeight','bold');
        else
            title(sprintf('CWT   \\Delta t = %.3f \\mus', dt_us_all(i)), 'FontWeight','bold');
        end

        hold on;

        % 叠加频散曲线（A0白实线 / S0白虚线；白字图例透明无框）
        if curveOpt.enable && ~isempty(curveDB)
            overlayDispersionOnCWT(gca, curveDB, curveOpt, start_us, t_us(end), dispOpt.fmax_kHz);
        end

        % AIC拾取线（仅AIC开启时显示；用白线与CWT风格一致）
        if showPickLines
            xline(t1_on_us,'w--','LineWidth',1);
            xline(t2_on_us,'w--','LineWidth',1);
        end

        hold off;
    end
end

%% =============== 6) 输出结果表 ===============
pickMethod = repmat({iff(pickOpt.useAIC,'AIC','Peak')}, nEv, 1);

T = table(evList(:), repmat(chSel,nEv,1), ...
    t1_us_all, t2_us_all, dt_us_all, ...
    t1_peak_us_all, t2_peak_us_all, dt_peak_us_all, ...
    pkt2_segmented, peak2_ratio, pickMethod, ...
    'VariableNames', {'event','ch','t1_us','t2_us','dt_us', ...
                      't1_peak_us','t2_peak_us','dt_peak_us', ...
                      'pkt2_segmented','peak2_ratio','pick_method'});

fprintf('\n========== AIC / Peak Δt results ==========\n');
fprintf('Timing method used for t1/t2/dt: %s\n', iff(pickOpt.useAIC,'AIC onset','Envelope peak'));
disp(T);

%% ===================== local functions =====================

function [wave, fs] = readWaveRobust(wavefile)
% 功能：兼容 waveReader 不同输出个数（7或13个）
    try
        out = cell(1,13);
        [out{:}] = waveReader(wavefile);
        wave = out{2};
        fs   = out{7};
        if isempty(fs), error('fs empty'); end
    catch
        out = cell(1,7);
        [out{:}] = waveReader(wavefile);
        wave = out{2};
        fs   = out{7};
    end
end

function [chSel, evList, sigList] = pickMultiEvents_from_wave(wave)
% 功能：根据 wave 维度选择通道/事件，输出每个事件的单列波形
% 支持：
%   1D: N
%   2D: N x ch
%   3D: N x ch x event
    chSel = 1;
    sigList = {};

    if isvector(wave)
        evList = 1;
        sigList{1} = wave(:);
        return;
    end

    nd = ndims(wave);

    if nd == 2
        [~,C] = size(wave);
        if C > 1
            listC = compose('ch %d', 1:C);
            [chSel, ok] = listdlg('ListString',listC,'SelectionMode','single', ...
                'PromptString','选择通道 channel（单选）', 'InitialValue',1,'ListSize',[220 260]);
            if isempty(ok) || ok==0, chSel = 1; end
        end
        evList = 1;
        sigList{1} = wave(:,chSel);
        return;
    end

    % nd >= 3
    C  = size(wave,2);
    Ev = size(wave,3);

    if C > 1
        listC = compose('ch %d', 1:C);
        [chSel, ok] = listdlg('ListString',listC,'SelectionMode','single', ...
            'PromptString','选择通道 channel（单选）', 'InitialValue',1,'ListSize',[220 260]);
        if isempty(ok) || ok==0, chSel = 1; end
    end

    if Ev > 1
        listE = compose('event %d', 1:Ev);
        [evList, ok] = listdlg('ListString',listE,'SelectionMode','multiple', ...
            'PromptString','选择事件 event（可多选：Ctrl/Shift）', ...
            'InitialValue',1,'ListSize',[240 320]);
        if isempty(ok) || ok==0
            evList = 1;
        end
    else
        evList = 1;
    end

    sigList = cell(numel(evList),1);
    for i = 1:numel(evList)
        ev = evList(i);
        sig = squeeze(wave(:, chSel, ev));
        sigList{i} = sig(:);
    end
end

function [seg1, k1, seg2, k2, pkt2_segmented, peak2_ratio] = autoFindTwoPackets_anchored(env_s, fs, det)
% 功能：自动找两个波包（第二包限定在第一包后的时间窗内）
% 输出：
%   seg1/seg2 = [start end] 段
%   k1/k2     = 段内峰值点索引
    env_s = env_s(:);
    N = numel(env_s);

    noiseN = max(5, round(det.noise_us*1e-6*fs));
    noiseN = min(noiseN, N);

    med0 = median(env_s(1:noiseN));
    mad0 = median(abs(env_s(1:noiseN) - med0)) + eps;

    % ---- 找第一包（最早显著段）----
    idxSearch1 = (noiseN+1):N;
    seg1 = []; k1 = [];

    for K = det.thrK_list
        thr = med0 + K*mad0;
        segs = segmentAboveThr(env_s, idxSearch1, thr, fs, det.mergeGap_us, det.minDur_us);
        if ~isempty(segs)
            seg1 = segs(1,:);
            [~,r] = max(env_s(seg1(1):seg1(2)));
            k1 = seg1(1) + r - 1;
            break;
        end
    end

    if isempty(seg1)
        [~,k1] = max(env_s);  % 兜底：全局峰
        seg1 = expandLocalSegment(env_s, k1, fs, det.minDur_us);
    end

    % ---- 第二包搜索窗 ----
    minSepN   = round(det.minSep_us*1e-6*fs);
    minDelayN = round(det.pkt2_minDelay_us*1e-6*fs);
    maxDelayN = round(det.pkt2_maxDelay_us*1e-6*fs);

    s2 = min(N, k1 + max(minSepN, minDelayN));
    e2 = min(N, k1 + maxDelayN);
    if e2 <= s2
        s2 = min(N, k1 + minSepN);
        e2 = N;
    end
    idxSearch2 = s2:e2;

    p1 = max(env_s(seg1(1):seg1(2)));  % 第一包峰值（用于ratio）
    seg2 = []; k2 = [];
    pkt2_segmented = false;

    for K = det.thrK_list
        thr = med0 + K*mad0;
        segs2 = segmentAboveThr(env_s, idxSearch2, thr, fs, det.mergeGap_us, det.minDur_us);
        if isempty(segs2), continue; end

        nS = size(segs2,1);
        peaks = zeros(nS,1);
        kpk   = zeros(nS,1);

        for j = 1:nS
            sj = segs2(j,1); ej = segs2(j,2);
            [peaks(j), rr] = max(env_s(sj:ej));
            kpk(j) = sj + rr - 1;
        end

        % 优先“最早且足够强”的段
        [~,ord] = sort(kpk,'ascend');
        found = false;
        for jj = 1:numel(ord)
            j = ord(jj);
            if (kpk(j)-k1) < minSepN, continue; end
            if peaks(j) >= det.pkt2_minPeakRatio * p1
                seg2 = segs2(j,:);
                k2   = kpk(j);
                pkt2_segmented = true;
                found = true;
                break;
            end
        end
        if found, break; end

        % 否则按峰值×能量评分选最佳段
        E = zeros(nS,1);
        for j = 1:nS
            sj = segs2(j,1); ej = segs2(j,2);
            E(j) = sum(env_s(sj:ej).^2);
        end
        score = peaks .* sqrt(E + eps);
        [~,jBest] = max(score);
        seg2 = segs2(jBest,:);
        k2   = kpk(jBest);
        pkt2_segmented = true;
        break;
    end

    if isempty(seg2)
        % 兜底：搜索窗内最大峰
        [~,rr] = max(env_s(idxSearch2));
        k2 = idxSearch2(1) + rr - 1;
        seg2 = expandLocalSegment(env_s, k2, fs, det.minDur_us);
        pkt2_segmented = false;
    end

    p2 = max(env_s(seg2(1):seg2(2)));
    peak2_ratio = p2 / (p1 + eps);
end

function segs = segmentAboveThr(x, idxRange, thr, fs, mergeGap_us, minDur_us)
% 功能：在 idxRange 内对 x>thr 分段，并进行段合并/短段剔除
    N = numel(x);
    mask = false(N,1);
    mask(idxRange) = x(idxRange) > thr;

    d = diff([false; mask; false]);
    st = find(d==1);
    en = find(d==-1)-1;

    if isempty(st)
        segs = [];
        return;
    end

    % 合并近邻段
    mergeGapN = round(mergeGap_us*1e-6*fs);
    k = 1;
    while k < numel(st)
        gap = st(k+1) - en(k) - 1;
        if gap <= mergeGapN
            en(k) = en(k+1);
            st(k+1) = [];
            en(k+1) = [];
        else
            k = k + 1;
        end
    end

    % 删除过短段
    minDurN = round(minDur_us*1e-6*fs);
    dur = en - st + 1;
    keep = dur >= max(6, minDurN);

    st = st(keep);
    en = en(keep);

    segs = [st(:) en(:)];
end

function seg = expandLocalSegment(env_s, k0, fs, minDur_us)
% 功能：兜底时，以峰值为中心扩一个最小长度段
    N = numel(env_s);
    halfN = round(0.5*minDur_us*1e-6*fs);
    halfN = max(6, halfN);
    s = max(1, k0-halfN);
    e = min(N, k0+halfN);
    seg = [s e];
end

function t_on_s = aicOnsetTime_np_limitedRise(x, k_peak, preN, postN, fs, segStart, segEnd)
% 功能：AIC onset（非参数）
% 限制：只在“峰值之前的上升段”搜索，避免跑到波包之后
    x = x(:);
    N = numel(x);

    % 初始窗口（围绕峰值）
    i1 = max(1, k_peak - preN);
    i2 = min(N, k_peak + postN);

    % 限制在该段内
    i1 = max(i1, segStart);
    i2 = min(i2, segEnd);

    % 关键：只允许到峰值（rise only）
    i2 = min(i2, k_peak);

    if i2 - i1 + 1 < 12
        t_on_s = (max(1,i1)-1)/fs;
        return;
    end

    seg = x(i1:i2);
    M = numel(seg);

    AIC = inf(M,1);
    for m = 2:(M-2)
        v1 = var(seg(1:m),1);      v1 = max(v1, eps);
        v2 = var(seg(m+1:end),1);  v2 = max(v2, eps);
        AIC(m) = m*log(v1) + (M-m-1)*log(v2);
    end

    [~,m0] = min(AIC);
    k_on = i1 + m0 - 1;
    t_on_s = (k_on-1)/fs;
end

function [f_lin_Hz, f_kHz_axis, img_show] = cwtImage_EMDstyle(sig, fs, dispOpt)
% 功能：CWT(amor) + 线性频轴插值 + 归一化 + 阈值 + gamma
    sig = sig(:);

    f_lin_Hz   = linspace(0, dispOpt.fmax_kHz*1000, dispOpt.nFreqBins);
    f_kHz_axis = f_lin_Hz / 1000;

    [cfs, f] = cwt(sig, fs, 'amor', 'VoicesPerOctave', dispOpt.voicesPerOctave);
    mag = abs(cfs);

    % 频率升序
    if f(1) > f(end)
        f = flipud(f);
        mag = flipud(mag);
    end

    [f_u, iu] = unique(f, 'stable');
    mag_u = mag(iu,:);

    % 保证 0Hz 可插值
    if f_u(1) > 0
        f_u   = [0; f_u];
        mag_u = [zeros(1,size(mag_u,2)); mag_u];
    end

    mag_lin = interp1(f_u, mag_u, f_lin_Hz, 'linear', 0);

    mmax = max(mag_lin(:));
    if mmax <= 0, mmax = 1; end

    img_show = mag_lin ./ (mmax + eps);
    img_show(img_show < dispOpt.threshold) = 0;
    img_show(img_show > 1) = 1;
    img_show = img_show .^ dispOpt.gamma;
end

function out = iff(cond, a, b)
% 功能：简易三目运算
    if cond
        out = a;
    else
        out = b;
    end
end

function curveDB = readDispersionExcelVallen(xlsFile, verboseFlag)
% 功能：读取 Vallen 导出的频散曲线 Excel
% 常见格式（例如 3mmPMMA_Dispersion Curve.xlsx）：
%   - 某一行出现 "Group A0" / "Group S0" / "Phase A0" ...
%   - 后续行为数值列（频率+速度）
%
% 输出 struct 数组字段：
%   .kind  = 'group' / 'phase'
%   .mode  = 'A0','S0',...
%   .f_MHz
%   .v_m_per_ms   (数值上等于 mm/us)
    if nargin < 2, verboseFlag = false; end
    curveDB = struct('kind',{},'mode',{},'f_MHz',{},'v_m_per_ms',{});

    C = [];
    try
        C = readcell(xlsFile, 'Sheet', 1);
    catch
        try
            [~,~,raw] = xlsread(xlsFile, 1);
            C = raw;
        catch ME
            warning('读取 Excel 失败: %s', ME.message);
            return;
        end
    end

    if isempty(C) || size(C,1) < 5
        warning('Excel 内容为空或格式不符合预期。');
        return;
    end

    headerRowsTry = [8, 7, 6, 5, 4, 3, 2, 1];  % 常见Vallen导出标题行
    parsed = false;

    for hRow = headerRowsTry
        if hRow > size(C,1), continue; end
        tmp = parseCurveFromHeaderRow(C, hRow, verboseFlag);
        if ~isempty(tmp)
            curveDB = tmp;
            parsed = true;
            break;
        end
    end

    if ~parsed
        warning('未识别到 "Group/Phase + 模态名" 列头，请检查 Excel 格式。');
    end
end

function curveDB = parseCurveFromHeaderRow(C, headerRow, verboseFlag)
% 功能：从指定标题行扫描 "Group A0" / "Phase S0" 等列头
    curveDB = struct('kind',{},'mode',{},'f_MHz',{},'v_m_per_ms',{});
    nrow = size(C,1);
    ncol = size(C,2);

    dataRow = min(nrow, headerRow + 2);  % 通常下一行是单位，再下一行起为数值

    for c = 1:(ncol-1)
        h = C{headerRow, c};
        if ~(ischar(h) || isstring(h)), continue; end
        h = strtrim(char(h));

        tok = regexp(h, '^(Group|Phase)\s+([SA]\d+)$', 'tokens', 'once', 'ignorecase');
        if isempty(tok), continue; end

        kind = lower(tok{1});
        mode = upper(tok{2});

        % 默认假设成对列：(F, v)
        fcol = cell2numCol(C(dataRow:end, c));
        vcol = cell2numCol(C(dataRow:end, min(c+1,ncol)));

        % 如果点数太少，尝试从下一行再读一次
        if nnz(isfinite(fcol) & isfinite(vcol)) < 3 && (dataRow+1)<=nrow
            fcol = cell2numCol(C((dataRow+1):end, c));
            vcol = cell2numCol(C((dataRow+1):end, min(c+1,ncol)));
        end

        keep = isfinite(fcol) & isfinite(vcol) & (fcol > 0) & (vcol > 0);
        fcol = fcol(keep);
        vcol = vcol(keep);

        if isempty(fcol), continue; end

        [fcol, ord] = sort(fcol(:), 'ascend');
        vcol = vcol(ord);

        [fcol_u, iu] = unique(fcol, 'stable');
        vcol_u = vcol(iu);

        curveDB(end+1).kind = kind; %#ok<AGROW>
        curveDB(end).mode = mode;
        curveDB(end).f_MHz = fcol_u;
        curveDB(end).v_m_per_ms = vcol_u;

        if verboseFlag
            fprintf('[DispCurve] row=%d  %-5s %-3s : %d points\n', ...
                headerRow, upper(kind), mode, numel(fcol_u));
        end
    end
end

function x = cell2numCol(cc)
% 功能：将 readcell/xlsread 返回的 cell 列安全转为 double 列
    n = numel(cc);
    x = nan(n,1);

    for i = 1:n
        v = cc{i};
        if isnumeric(v) && isscalar(v)
            x(i) = double(v);
        elseif ischar(v) || isstring(v)
            t = str2double(strtrim(string(v)));
            if isfinite(t), x(i) = t; end
        end
    end
end

function overlayDispersionOnCWT(ax, curveDB, curveOpt, start_us, xEnd_us, fmax_kHz)
% 功能：在 CWT 图上叠加频散曲线（白线 + 白字图例）
% 时间换算（group/phase 均用）：
%   t_us = distance_mm / v(mm/us) + timeShift_us - start_us
%
% 图例样式（按你的要求）：
%   - 白字
%   - 透明背景（无黑底）
%   - 无边框（避免PPT出现白框）
%   - 固定右上角
    if isempty(curveDB), return; end

    modeList = normalizeModeList(curveOpt.modes);
    if isempty(modeList), return; end

    kindWant = lower(strtrim(curveOpt.type));
    if ~ismember(kindWant, {'group','phase'})
        warning('curveOpt.type 必须是 ''group'' 或 ''phase''。当前=%s', kindWant);
        return;
    end

    lgdHandles = gobjects(0);
    lgdNames   = {};

    for k = 1:numel(curveDB)
        if ~strcmpi(curveDB(k).kind, kindWant), continue; end
        if ~any(strcmpi(curveDB(k).mode, modeList)), continue; end

        f_kHz = curveDB(k).f_MHz(:) * 1000;      % MHz -> kHz
        v     = curveDB(k).v_m_per_ms(:);        % [m/ms] == [mm/us]（数值等价）

        t_theory_us = curveOpt.distance_mm ./ v; % us
        t_plot_us   = t_theory_us + curveOpt.timeShift_us - start_us;

        keep = isfinite(t_plot_us) & isfinite(f_kHz) & ...
               (t_plot_us >= 0) & (t_plot_us <= xEnd_us) & ...
               (f_kHz >= 0) & (f_kHz <= fmax_kHz);

        if nnz(keep) < 2, continue; end

        t_plot = t_plot_us(keep);
        f_plot = f_kHz(keep);

        % 按频率排序，视觉更稳定
        [f_plot, ord] = sort(f_plot, 'ascend');
        t_plot = t_plot(ord);

        [clr, ls] = modeLineStyle(curveDB(k).mode, kindWant);

        h = plot(ax, t_plot, f_plot, ...
                 'LineStyle', ls, ...
                 'Color', clr, ...
                 'LineWidth', max(curveOpt.lineWidth, 2.0), ...
                 'HandleVisibility', 'on');

        lgdHandles(end+1) = h; %#ok<AGROW>
        lgdNames{end+1}   = upper(curveDB(k).mode); %#ok<AGROW>  % 只显示 A0 / S0
    end

    if curveOpt.showLegend && ~isempty(lgdHandles)
        [lgdNamesU, ia] = unique(lgdNames, 'stable');

        % 固定图例顺序：A0（实线）在前，S0（虚线）在后
        prefOrder = {'A0','S0'};
        rank = inf(size(lgdNamesU));
        for ii = 1:numel(lgdNamesU)
            jj = find(strcmpi(lgdNamesU{ii}, prefOrder), 1);
            if ~isempty(jj), rank(ii) = jj; end
        end
        [~, so] = sort(rank, 'ascend');
        lgdNamesU = lgdNamesU(so);
        ia = ia(so);

        % 创建图例
        try
            lgd = legend(ax, lgdHandles(ia), lgdNamesU, ...
                'Location','northeast', ...
                'Interpreter','none', ...
                'AutoUpdate','off');
        catch
            lgd = legend(ax, lgdHandles(ia), lgdNamesU, ...
                'Location','northeast', ...
                'Interpreter','none');
        end

        % 白字 + 透明背景 + 无边框（避免PPT白框）
        try
            lgd.TextColor = 'w';
            lgd.Color     = 'none';
            lgd.EdgeColor = 'none';
            lgd.Box       = 'off';
            lgd.FontSize  = 10;
        catch
            % 老版本兼容
            try
                set(lgd, 'TextColor','w', 'Color','none', 'EdgeColor','none', 'Box','off');
            catch
                % 再不支持时至少保证白字 + 去框
                set(lgd, 'TextColor','w', 'Box','off');
            end
        end
    end
end

function modeList = normalizeModeList(modesIn)
% 功能：将 char/string/cellstr 统一转成 cellstr（大写）
    if ischar(modesIn)
        modeList = {upper(strtrim(modesIn))};
    elseif isstring(modesIn)
        modeList = cellstr(upper(strtrim(modesIn(:))));
    elseif iscell(modesIn)
        modeList = cell(size(modesIn));
        for i = 1:numel(modesIn)
            modeList{i} = upper(strtrim(char(modesIn{i})));
        end
    else
        modeList = {};
    end
end

function [clr, ls] = modeLineStyle(modeName, kindName)
% 功能：频散曲线样式（按你的要求）
%   A0 = 白色实线
%   S0 = 白色虚线
%   其他 = 白色点划线/点线（兜底）
    %#ok<INUSD>
    modeName = upper(strtrim(char(modeName)));
    clr = [1 1 1];   % 全白

    switch modeName
        case 'A0'
            ls = '-';
        case 'S0'
            ls = '--';
        case 'A1'
            ls = '-.';
        case 'S1'
            ls = ':';
        otherwise
            ls = '-';
    end
end