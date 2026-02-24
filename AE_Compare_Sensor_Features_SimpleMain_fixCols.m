%% AE_Compare_Sensor_FeatureSpace_Full_withBrowser.m
% -------------------------------------------------------------------------
% 功能（完整整合版）
%  1) 读取一个或多个 .wave 文件（不同距离）
%  2) 用你的工具链计算特征（AIC + 99%能量终点 + AE_Calc_Features）
%  3) 幅值异常剔除（简单 robust MAD）
%  4) 绘制二维特征空间图：
%       - 颜色 = 距离
%       - 标记形状 = 传感器（Ch1=B1025 -> ○, Ch2=B454 -> △）
%       - 横纵坐标 = 两个不同特征值
%  5) 绘制总览拼图（已兼容旧版 MATLAB，修复 legend(tlo,...) 报错）
%  6) 启动交互浏览器（自定义 X/Y、坐标范围、距离勾选、连线开关）
%
% 依赖（需自行 addpath）：
%   - waveReader
%   - AE_Tools_V3
%   - AE_Calc_Features
%   - ExtractFeati_OPT
%
% 注意：
%   - 按你的要求：不再验证 waveReader 路径
%   - AE_Calc_Features 输出列名若不同，脚本已做别名兼容
% -------------------------------------------------------------------------

clear; clc; close all;
rng(11);

%% ====================== 用户配置 ======================
cfg = struct();

% ---- 文件 ----
cfg.use_uigetfile = true;   % true=交互选择 wave 文件；false=使用 cfg.files
cfg.files = {
    % 'F:\Vallen\Distance effect\2026-2-12-2sensor\center_100mm-1.wave'
    % 'F:\Vallen\Distance effect\2026-2-12-2sensor\center_150mm-1.wave'
    % 'F:\Vallen\Distance effect\2026-2-12-2sensor\center_200mm-1.wave'
};

% 若留空则自动从文件名提取 xxmm
cfg.distance_mm = [];   % e.g. [100 150 200]

% ---- 传感器与通道 ----
cfg.ch1 = 1;              % B1025
cfg.ch2 = 2;              % B454
cfg.sensor1 = 'B1025';
cfg.sensor2 = 'B454';

% ---- AIC / 分段参数 ----
cfg.prepad_us = 20;
cfg.endEnergyFrac = 0.99;
cfg.AIC = struct('thr_main_factor',7,'dt_rel_margin',0.2,'deep_all',true,'RefineFcn',[]);
cfg.GEO = [];

% ---- 主要特征对（横轴,纵轴）----
cfg.plotFeaturePairs = { ...
    {'PF_kHz','FC_kHz'}; ...
    {'Amp_dB','PF_kHz'}; ...
    {'Amp_dB','FC_kHz'}; ...
    {'RiseTime_us','Amp_dB'}; ...
    % {'Dur_us','Energy_V2'};  % 可选（Energy 适合 log）
};

% ---- 幅值异常剔除（简单版）----
flt = struct();
flt.abs_amp_range_dB = [25, 95];
flt.amp_mad_k        = 4.0;
flt.pairdiff_mad_k   = 4.5;
flt.min_win_us       = 20;
cfg.filter = flt;

% ---- 绘图样式 ----
cfg.sensorMarkers     = {'o','^'};   % B1025=○, B454=△
cfg.sensorMarkerSize  = [34, 42];
cfg.sensorLineWidth   = 1.0;

cfg.showPairedLine    = false;       % 同事件双传感器点连线
cfg.pairLineWidth     = 0.5;
cfg.pairLineLightness = 0.35;

cfg.useLogX_forEnergy = true;
cfg.useLogY_forEnergy = true;

cfg.figPosSingle = [100 100 980 650];
cfg.figPosPanel  = [80 60 1450 900];
cfg.panelCols    = 2;

% ---- 输出 ----
cfg.outDir = fullfile(pwd, 'AE_sensor_feature_space_out');
cfg.saveFig = true;
cfg.figExts = {'.png','.fig'};  % 可加 '.pdf'
cfg.exportDPI = 220;

% ---- 运行结束后自动打开交互浏览器 ----
cfg.launchBrowser = true;

%% ====================== 依赖检查（不检查 waveReader） ======================
need = {'AE_Tools_V3','AE_Calc_Features','ExtractFeati_OPT'};
for i = 1:numel(need)
    if exist(need{i},'file') ~= 2 && exist(need{i},'class') ~= 8
        error('缺少依赖：%s（请确认已 addpath）', need{i});
    end
end
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

%% ====================== 选择文件 ======================
filePaths = cfg.files;
if cfg.use_uigetfile
    [f,p] = uigetfile({'*.wave','Hexagon wave (*.wave)'}, ...
        '选择一个或多个 .wave 文件（可多选不同距离）', 'MultiSelect','on');
    if isequal(f,0), error('未选择文件。'); end
    if ischar(f), f = {f}; end
    filePaths = cellfun(@(x) fullfile(p,x), f(:), 'UniformOutput', false);
end
if isempty(filePaths), error('没有输入文件。'); end
nFiles = numel(filePaths);

% 距离：优先从文件名提取 xxmm
if isempty(cfg.distance_mm)
    d = nan(nFiles,1);
    for i = 1:nFiles
        [~,nm,~] = fileparts(filePaths{i});
        tok = regexp(nm,'(?i)(\d+(?:\.\d+)?)\s*mm','tokens','once');
        if ~isempty(tok), d(i) = str2double(tok{1}); end
    end
    if any(~isfinite(d))
        prompt = arrayfun(@(k)sprintf('文件 #%d 的距离 (mm):',k), 1:nFiles, 'UniformOutput', false);
        answ = inputdlg(prompt, '输入距离', [1 40], repmat({''},1,nFiles));
        if isempty(answ), error('已取消距离输入。'); end
        d = cellfun(@str2double, answ(:));
    end
    cfg.distance_mm = d(:)';
end
if numel(cfg.distance_mm) ~= nFiles
    error('cfg.distance_mm 数量与文件数量不一致。');
end

fprintf('\n===== 文件列表 =====\n');
for i = 1:nFiles
    fprintf('[%d/%d] %s  (distance=%g mm)\n', i, nFiles, filePaths{i}, cfg.distance_mm(i));
end

%% ====================== 主循环：读取 + 分段 + 特征提取 ======================
Tall = table();
Tfile = table();

for iFile = 1:nFiles
    fpath = filePaths{iFile};
    distMM = cfg.distance_mm(iFile);
    [~,fname,~] = fileparts(fpath);

    fprintf('\n[%d/%d] Processing: %s\n', iFile, nFiles, fpath);

    % 读取 wave（按你的 waveReader 输出顺序）
    [~, wave, ~, ~, eventTime, ~, fs, ~, ~, ~, ~, ~, preTrigPoints] = waveReader(fpath);
    wave = double(wave);
    eventTime = eventTime(:);

    [Ns,Nch,Nev] = size(wave);
    if max([cfg.ch1 cfg.ch2]) > Nch
        error('文件 %s 只有 %d 个通道，但请求了 [%d %d]。', fname, Nch, cfg.ch1, cfg.ch2);
    end

    % AIC 起点（双通道）
    CH = struct('S1',cfg.ch1,'S2',cfg.ch2);
    if isempty(cfg.GEO), GEO = struct(); else, GEO = cfg.GEO; end
    A = AE_Tools_V3.pick_all_events_hybrid(wave, fs, CH, preTrigPoints, GEO, cfg.AIC);

    % 统一窗口：min(AIC)前补 + 99%能量终点
    prepad = round(cfg.prepad_us * 1e-6 * fs);

    k0_all = min(A.idx, [], 2);
    k0_all(~isfinite(k0_all)) = preTrigPoints + 1;
    k0_all = max(1, k0_all - prepad);

    xsum_all = squeeze(sum(wave(:, [cfg.ch1 cfg.ch2], :).^2, 2));  % [Ns x Nev]
    if isvector(xsum_all), xsum_all = xsum_all(:); end

    kE_all = nan(Nev,1);
    for ev = 1:Nev
        kE_all(ev) = AE_Tools_V3.find_end_by_energy(xsum_all(:,ev), fs, k0_all(ev), cfg.endEnergyFrac);
    end
    kE_all = min(max(kE_all, k0_all), Ns);

    % 构造基础表 E，并调用 AE_Calc_Features
    E = table((1:Nev)','VariableNames',{'eventID'});
    E.hittime  = eventTime;
    E.k0       = k0_all;
    E.kE       = kE_all;
    E.t0_S1_us = A.t0_S1_us(:);
    E.t0_S2_us = A.t0_S2_us(:);
    E.dt_us    = A.dt_us(:);
    if isfield(A,'valid')
        E.loc_valid = logical(A.valid(:));
    else
        E.loc_valid = false(Nev,1);
    end

    [~, Efeat] = AE_Calc_Features.run_and_merge( ...
        wave, fs, cfg.ch1, cfg.ch2, eventTime, k0_all, kE_all, E, true);

    if isempty(Efeat) || ~istable(Efeat)
        warning('AE_Calc_Features 未返回有效表，跳过文件：%s', fpath);
        continue;
    end

    % 元信息
    Efeat.distance_mm = repmat(distMM, height(Efeat), 1);
    Efeat.file_id     = repmat(iFile, height(Efeat), 1);
    Efeat.file_name   = repmat(string(fname), height(Efeat), 1);
    Efeat.fs_Hz       = repmat(fs, height(Efeat), 1);

    % 窗口长度过滤（可选）
    if ismember('win_us', Efeat.Properties.VariableNames) && ~isempty(cfg.filter.min_win_us)
        Efeat = Efeat(isfinite(Efeat.win_us) & Efeat.win_us >= cfg.filter.min_win_us, :);
    end

    % 幅值异常过滤（只做简单排除）
    featAliasMap = get_feature_alias_map();
    [amp1col, amp2col] = resolve_pair_cols(Efeat, 'Amp_dB', featAliasMap);
    if ~isempty(amp1col) && ~isempty(amp2col)
        Efeat = local_filter_by_amp(Efeat, cfg.filter, amp1col, amp2col);
    else
        warning('未找到幅值列（Amp_dB），本文件不做幅值异常过滤：%s', fname);
    end

    % 文件级统计
    Tfile = [Tfile; table(iFile, string(fname), distMM, fs, Nev, height(Efeat), ...
        'VariableNames',{'file_id','file_name','distance_mm','fs_Hz','nEvents_raw','nEvents_kept'})]; %#ok<AGROW>

    Tall = [Tall; Efeat]; %#ok<AGROW>
end

if isempty(Tall)
    error('过滤后没有可用事件。请放宽幅值过滤阈值或检查数据。');
end

writetable(Tfile, fullfile(cfg.outDir,'file_summary.csv'));

fprintf('\n合并后有效事件数：%d\n', height(Tall));
disp('部分列名预览：');
disp(Tall.Properties.VariableNames(1:min(30, numel(Tall.Properties.VariableNames)))');

%% ====================== 单图 & 总览拼图：特征空间 ======================
featAliasMap = get_feature_alias_map();

dListGlobal = unique(Tall.distance_mm(isfinite(Tall.distance_mm)))';
if isempty(dListGlobal), dListGlobal = []; end
cmap = lines(max(3, numel(dListGlobal)));

panelInfo = struct('featX',{},'featY',{},'T',{},'colNames',{});

for ip = 1:size(cfg.plotFeaturePairs,1)
    pairCell = cfg.plotFeaturePairs{ip};
    featX = pairCell{1};
    featY = pairCell{2};

    [c1x, c2x] = resolve_pair_cols(Tall, featX, featAliasMap);
    [c1y, c2y] = resolve_pair_cols(Tall, featY, featAliasMap);

    if isempty(c1x) || isempty(c2x) || isempty(c1y) || isempty(c2y)
        warning('跳过特征对 (%s, %s)：找不到列。', featX, featY);
        continue;
    end

    needCols = {'distance_mm','file_id','file_name', c1x,c1y,c2x,c2y};
    if ismember('eventID', Tall.Properties.VariableNames)
        needCols = [{'eventID'}, needCols];
    end
    T = Tall(:, needCols);

    % 四列都有效
    m = isfinite(T.(c1x)) & isfinite(T.(c1y)) & isfinite(T.(c2x)) & isfinite(T.(c2y));
    T = T(m,:);
    if height(T) < 3
        warning('跳过特征对 (%s, %s)：有效点太少。', featX, featY);
        continue;
    end

    % log轴前过滤能量非正值
    if contains(c1x, 'Energy', 'IgnoreCase', true) || contains(c2x, 'Energy', 'IgnoreCase', true)
        m2 = T.(c1x)>0 & T.(c2x)>0;
        T = T(m2,:);
    end
    if contains(c1y, 'Energy', 'IgnoreCase', true) || contains(c2y, 'Energy', 'IgnoreCase', true)
        m2 = T.(c1y)>0 & T.(c2y)>0;
        T = T(m2,:);
    end
    if height(T) < 3
        warning('跳过特征对 (%s, %s)：过滤后点太少。', featX, featY);
        continue;
    end

    % ---------- 单图 ----------
    f = figure('Color','w', ...
        'Name', sprintf('FeatureSpace_%s_vs_%s', featY, featX), ...
        'Position', cfg.figPosSingle);
    ax = axes(f); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    % 可选：同事件双传感器连线
    if cfg.showPairedLine
        for ii = 1:height(T)
            di = T.distance_mm(ii);
            ic = find(dListGlobal == di, 1, 'first');
            if isempty(ic), ic = 1; end
            col = cmap(ic,:);
            colLight = lighten_color(col, cfg.pairLineLightness);
            plot(ax, [T.(c1x)(ii), T.(c2x)(ii)], [T.(c1y)(ii), T.(c2y)(ii)], '-', ...
                'Color', colLight, 'LineWidth', cfg.pairLineWidth, 'HandleVisibility','off');
        end
    end

    dList = unique(T.distance_mm(isfinite(T.distance_mm)))';
    if isempty(dList), dList = NaN; end
    hDist = gobjects(numel(dList),1);

    for idd = 1:numel(dList)
        di = dList(idd);
        idx = (T.distance_mm == di);
        if ~any(idx), continue; end

        ic = find(dListGlobal == di, 1, 'first');
        if isempty(ic), ic = 1; end
        col = cmap(ic,:);

        % Ch1 / B1025 : ○
        scatter(ax, T.(c1x)(idx), T.(c1y)(idx), cfg.sensorMarkerSize(1), ...
            cfg.sensorMarkers{1}, ...
            'MarkerEdgeColor', col, 'MarkerFaceColor', 'none', ...
            'LineWidth', cfg.sensorLineWidth, 'HandleVisibility','off');

        % Ch2 / B454 : △
        scatter(ax, T.(c2x)(idx), T.(c2y)(idx), cfg.sensorMarkerSize(2), ...
            cfg.sensorMarkers{2}, ...
            'MarkerEdgeColor', col, 'MarkerFaceColor', 'none', ...
            'LineWidth', cfg.sensorLineWidth, 'HandleVisibility','off');

        % 距离图例（只显示颜色）
        hDist(idd) = plot(ax, nan, nan, '-', 'Color', col, 'LineWidth', 2.2, ...
            'DisplayName', sprintf('%g mm (n=%d)', di, sum(idx)));
    end

    % 传感器图例（只显示形状）
    hS1 = plot(ax, nan, nan, cfg.sensorMarkers{1}, ...
        'LineStyle','none', 'MarkerSize',7, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor','none', ...
        'DisplayName', sprintf('%s (Ch%d)', cfg.sensor1, cfg.ch1));
    hS2 = plot(ax, nan, nan, cfg.sensorMarkers{2}, ...
        'LineStyle','none', 'MarkerSize',8, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor','none', ...
        'DisplayName', sprintf('%s (Ch%d)', cfg.sensor2, cfg.ch2));

    xlabel(ax, featX, 'Interpreter','none');
    ylabel(ax, featY, 'Interpreter','none');

    if cfg.useLogX_forEnergy && contains(featX, 'Energy', 'IgnoreCase', true)
        set(ax,'XScale','log');
    end
    if cfg.useLogY_forEnergy && contains(featY, 'Energy', 'IgnoreCase', true)
        set(ax,'YScale','log');
    end

    title(ax, sprintf('%s vs %s in feature space (%s vs %s)\nColor = distance, Marker = sensor', ...
        cfg.sensor1, cfg.sensor2, featY, featX), ...
        'FontWeight','bold', 'Interpreter','none');

    lgHandles = [hS1, hS2, hDist(:).'];
    lgHandles = lgHandles(isgraphics(lgHandles));
    lgd = legend(ax, lgHandles, 'Location','best', 'FontSize',9, 'Interpreter','none');
    try
        lgd.Title.String = 'Marker / Distance';
        lgd.Title.FontWeight = 'normal';
    catch
    end

    set(ax, 'LineWidth', 1.0, 'FontSize', 10);

    if cfg.saveFig
        base = fullfile(cfg.outDir, sprintf('FeatureSpace_%s_vs_%s_%s_%s', ...
            sanitize_filename(featY), sanitize_filename(featX), ...
            cfg.sensor1, cfg.sensor2));
        save_figure_multi(f, base, cfg.figExts, cfg.exportDPI);
    end

    panelInfo(end+1).featX = featX; %#ok<SAGROW>
    panelInfo(end).featY   = featY;
    panelInfo(end).T       = T;
    panelInfo(end).colNames = struct('c1x',c1x,'c1y',c1y,'c2x',c2x,'c2y',c2y);
end

%% ====================== 总览拼图（修复 legend(tlo,...) 报错） ======================
if ~isempty(panelInfo)
    nP = numel(panelInfo);
    nCol = cfg.panelCols;
    nRow = ceil(nP / nCol);

    fpanel = figure('Color','w', 'Name','FeatureSpace_Panel', 'Position', cfg.figPosPanel);
    tlo = tiledlayout(nRow, nCol, 'Padding','compact', 'TileSpacing','compact');
    title(tlo, sprintf('Sensor comparison in feature space (%s vs %s)\nColor = distance, Marker = sensor', ...
        cfg.sensor1, cfg.sensor2), 'FontWeight','bold', 'Interpreter','none');

    hS1_global = gobjects(1);
    hS2_global = gobjects(1);
    hDist_global = gobjects(max(1, numel(dListGlobal)),1);

    for ip = 1:nP
        ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

        T = panelInfo(ip).T;
        featX = panelInfo(ip).featX;
        featY = panelInfo(ip).featY;
        cN = panelInfo(ip).colNames;

        if cfg.showPairedLine
            for ii = 1:height(T)
                di = T.distance_mm(ii);
                ic = find(dListGlobal == di, 1, 'first');
                if isempty(ic), ic = 1; end
                col = cmap(ic,:);
                colLight = lighten_color(col, cfg.pairLineLightness);
                plot(ax, [T.(cN.c1x)(ii), T.(cN.c2x)(ii)], [T.(cN.c1y)(ii), T.(cN.c2y)(ii)], '-', ...
                    'Color', colLight, 'LineWidth', cfg.pairLineWidth, 'HandleVisibility','off');
            end
        end

        dList = unique(T.distance_mm(isfinite(T.distance_mm)))';
        for idd = 1:numel(dList)
            di = dList(idd);
            idx = (T.distance_mm == di);
            if ~any(idx), continue; end
            ic = find(dListGlobal == di, 1, 'first');
            if isempty(ic), ic = 1; end
            col = cmap(ic,:);

            scatter(ax, T.(cN.c1x)(idx), T.(cN.c1y)(idx), cfg.sensorMarkerSize(1), ...
                cfg.sensorMarkers{1}, 'MarkerEdgeColor', col, 'MarkerFaceColor','none', ...
                'LineWidth', cfg.sensorLineWidth, 'HandleVisibility','off');

            scatter(ax, T.(cN.c2x)(idx), T.(cN.c2y)(idx), cfg.sensorMarkerSize(2), ...
                cfg.sensorMarkers{2}, 'MarkerEdgeColor', col, 'MarkerFaceColor','none', ...
                'LineWidth', cfg.sensorLineWidth, 'HandleVisibility','off');

            if ip == 1
                if ~isgraphics(hDist_global(idd))
                    hDist_global(idd) = plot(ax, nan, nan, '-', 'Color', col, 'LineWidth', 2.2, ...
                        'DisplayName', sprintf('%g mm', di));
                end
            end
        end

        if ip == 1
            hS1_global = plot(ax, nan, nan, cfg.sensorMarkers{1}, ...
                'LineStyle','none', 'MarkerSize',7, ...
                'MarkerEdgeColor','k', 'MarkerFaceColor','none', ...
                'DisplayName', sprintf('%s (Ch%d)', cfg.sensor1, cfg.ch1));
            hS2_global = plot(ax, nan, nan, cfg.sensorMarkers{2}, ...
                'LineStyle','none', 'MarkerSize',8, ...
                'MarkerEdgeColor','k', 'MarkerFaceColor','none', ...
                'DisplayName', sprintf('%s (Ch%d)', cfg.sensor2, cfg.ch2));
        end

        xlabel(ax, featX, 'Interpreter','none');
        ylabel(ax, featY, 'Interpreter','none');
        title(ax, sprintf('%s vs %s', featY, featX), 'Interpreter','none', 'FontWeight','bold');

        if cfg.useLogX_forEnergy && contains(featX, 'Energy', 'IgnoreCase', true)
            set(ax,'XScale','log');
        end
        if cfg.useLogY_forEnergy && contains(featY, 'Energy', 'IgnoreCase', true)
            set(ax,'YScale','log');
        end
    end

    % ===== 修复点：不再 legend(tlo,...)，改用某个 axes 建 legend =====
    lgHandles = [hS1_global, hS2_global, hDist_global(:).'];
    lgHandles = lgHandles(isgraphics(lgHandles));

    if ~isempty(lgHandles)
        axList = findobj(fpanel, 'Type', 'Axes');
        axList = axList(isgraphics(axList));
        if ~isempty(axList)
            axForLegend = axList(end);
        else
            axForLegend = gca;
        end

        lgd = legend(axForLegend, lgHandles, ...
            'Orientation','horizontal', ...
            'Location','southoutside', ...
            'Interpreter','none');

        try
            lgd.Layout.Tile = 'south';   % 新版 MATLAB 支持
        catch
            % 老版本不支持，保留 southoutside
        end

        try
            lgd.Title.String = 'Marker / Distance';
            lgd.Title.FontWeight = 'normal';
        catch
        end
    end

    if cfg.saveFig
        base = fullfile(cfg.outDir, sprintf('FeatureSpace_PANEL_%s_%s', cfg.sensor1, cfg.sensor2));
        save_figure_multi(fpanel, base, cfg.figExts, cfg.exportDPI);
    end
end

%% ====================== 数值摘要 CSV（按距离/特征） ======================
summaryRows = {};
allFeatNeeded = unique([cfg.plotFeaturePairs(:,1); cfg.plotFeaturePairs(:,2)]);

for iF = 1:numel(allFeatNeeded)
    feat = allFeatNeeded{iF};
    [c1, c2] = resolve_pair_cols(Tall, feat, featAliasMap);
    if isempty(c1) || isempty(c2), continue; end

    for di = unique(Tall.distance_mm(:))'
        idx = (Tall.distance_mm == di) & isfinite(Tall.(c1)) & isfinite(Tall.(c2));
        if ~any(idx), continue; end
        x = Tall.(c1)(idx);
        y = Tall.(c2)(idx);
        delta = y - x;

        if all(x > 0 & y > 0) && ~contains(feat,'Amp','IgnoreCase',true)
            relMed = median(100*(y-x)./max(eps,0.5*(x+y)), 'omitnan');
        else
            relMed = NaN;
        end

        summaryRows(end+1,:) = {string(feat), di, nnz(idx), ...
            median(x,'omitnan'), median(y,'omitnan'), ...
            mean(delta,'omitnan'), median(delta,'omitnan'), std(delta,'omitnan'), relMed}; %#ok<AGROW>
    end
end

if ~isempty(summaryRows)
    Tsum = cell2table(summaryRows, 'VariableNames', ...
        {'feature','distance_mm','N','median_S1','median_S2','mean_diff_S2minusS1', ...
         'median_diff_S2minusS1','sd_diff','median_relative_diff_pct'});
    writetable(Tsum, fullfile(cfg.outDir,'feature_space_summary_by_distance.csv'));
    disp('==== feature_space_summary_by_distance (preview) ====');
    disp(Tsum);
end

fprintf('\n完成。结果输出目录：%s\n', cfg.outDir);

% 放到工作区，便于后续自己调图
assignin('base', 'Tall_AE_feature_compare', Tall);
assignin('base', 'cfg_AE_feature_compare', cfg);

% 启动交互浏览器（自定义X/Y、范围、距离勾选）
if isfield(cfg,'launchBrowser') && cfg.launchBrowser
    try
        AE_FeatureSpace_Browser_local(Tall, cfg);
    catch ME
        warning('交互浏览器启动失败：%s', ME.message);
    end
end

%% ========================================================================
%% 局部函数
%% ========================================================================

function save_figure_multi(figH, base, figExts, dpiVal)
    for ie = 1:numel(figExts)
        ext = figExts{ie};
        try
            if strcmpi(ext,'.fig')
                savefig(figH, [base ext]);
            else
                exportgraphics(figH, [base ext], 'Resolution', dpiVal);
            end
        catch
            if ~strcmpi(ext,'.fig')
                saveas(figH, [base ext]);
            end
        end
    end
end

function featAliasMap = get_feature_alias_map()
    % 可按你的 AE_Calc_Features 输出继续扩展
    featAliasMap = struct();
    featAliasMap.Amp_dB      = {'Amp_dB','Amp','Amplitude_dB','PeakAmp_dB'};
    featAliasMap.FC_kHz      = {'FC_kHz','CentroidFreq_kHz','FreqCentroid_kHz','FC'};
    featAliasMap.PF_kHz      = {'PF_kHz','PeakFreq_kHz','PrimaryFreq_kHz','PF'};
    featAliasMap.Dur_us      = {'Dur_us','Duration_us','Dur'};
    featAliasMap.RiseTime_us = {'RiseTime_us','RT_us','Rise_us','RiseTime'};
    featAliasMap.Energy_V2   = {'Energy_V2','Energy','AE_Energy','SigEnergy_V2'};
    featAliasMap.Entropy     = {'Entropy','ShannonEntropy','SpecEntropy'};
end

function [c1, c2] = resolve_pair_cols(T, featCanonical, featAliasMap)
    c1 = ''; c2 = '';

    if ~isfield(featAliasMap, featCanonical)
        aliases = {featCanonical};
    else
        aliases = featAliasMap.(featCanonical);
    end
    vars = T.Properties.VariableNames;

    % 前缀风格
    p1 = {'S1_','Ch1_','CH1_','C1_','ch1_'};
    p2 = {'S2_','Ch2_','CH2_','C2_','ch2_'};

    c1 = resolve_col_by_alias(vars, p1, aliases);
    c2 = resolve_col_by_alias(vars, p2, aliases);

    % 后缀风格
    if isempty(c1)
        p1suf = {'_S1','_Ch1','_CH1','_C1','_ch1'};
        c1 = resolve_col_by_suffix_alias(vars, p1suf, aliases);
    end
    if isempty(c2)
        p2suf = {'_S2','_Ch2','_CH2','_C2','_ch2'};
        c2 = resolve_col_by_suffix_alias(vars, p2suf, aliases);
    end
end

function cname = resolve_col_by_alias(vars, prefixes, aliases)
    cname = '';
    for i = 1:numel(prefixes)
        for j = 1:numel(aliases)
            candidate = [prefixes{i}, aliases{j}];
            idx = find(strcmp(vars, candidate), 1, 'first');
            if ~isempty(idx)
                cname = vars{idx}; return;
            end
        end
    end
    for i = 1:numel(prefixes)
        for j = 1:numel(aliases)
            idx = find(startsWith(vars, prefixes{i}) & contains(vars, aliases{j}), 1, 'first');
            if ~isempty(idx)
                cname = vars{idx}; return;
            end
        end
    end
end

function cname = resolve_col_by_suffix_alias(vars, suffixes, aliases)
    cname = '';
    for j = 1:numel(aliases)
        for i = 1:numel(suffixes)
            candidate = [aliases{j}, suffixes{i}];
            idx = find(strcmp(vars, candidate), 1, 'first');
            if ~isempty(idx)
                cname = vars{idx}; return;
            end
        end
    end
    for j = 1:numel(aliases)
        for i = 1:numel(suffixes)
            idx = find(endsWith(vars, suffixes{i}) & contains(vars, aliases{j}), 1, 'first');
            if ~isempty(idx)
                cname = vars{idx}; return;
            end
        end
    end
end

function T = local_filter_by_amp(T, flt, col1, col2)
    % 仅基于幅值做简单异常剔除（你当前需求）
    if isempty(T), return; end
    if ~all(ismember({col1,col2}, T.Properties.VariableNames)), return; end

    a1 = T.(col1);
    a2 = T.(col2);
    m = isfinite(a1) & isfinite(a2);

    % 1) 绝对幅值范围
    if isfield(flt,'abs_amp_range_dB') && ~isempty(flt.abs_amp_range_dB)
        lo = flt.abs_amp_range_dB(1);
        hi = flt.abs_amp_range_dB(2);
        m = m & a1>=lo & a1<=hi & a2>=lo & a2<=hi;
    end

    % 2) 单通道 robust MAD
    if isfield(flt,'amp_mad_k') && ~isempty(flt.amp_mad_k)
        m = m & local_keep_mad(a1, flt.amp_mad_k) & local_keep_mad(a2, flt.amp_mad_k);
    end

    % 3) 配对幅值差 robust MAD
    if isfield(flt,'pairdiff_mad_k') && ~isempty(flt.pairdiff_mad_k)
        d = a2 - a1;
        m = m & local_keep_mad(d, flt.pairdiff_mad_k);
    end

    T = T(m,:);
end

function keep = local_keep_mad(x, k)
    keep = isfinite(x);
    xx = x(keep);
    if numel(xx) < 8 || isempty(k) || ~isfinite(k)
        return;
    end
    medx = median(xx,'omitnan');
    madx = median(abs(xx - medx), 'omitnan');
    if madx <= eps
        return;
    end
    z = abs(x - medx) / (1.4826 * madx); % robust z-score
    keep = keep & (z <= k);
end

function out = sanitize_filename(s)
    out = regexprep(s, '[^\w\.-]+', '_');
end

function c2 = lighten_color(c1, strength)
    strength = max(0,min(1,strength));
    c2 = c1 + (1-c1)*strength;
end

%% ====================== 交互浏览器（局部函数版） ======================
function AE_FeatureSpace_Browser_local(Tall, cfg)
    if nargin < 1 || isempty(Tall) || ~istable(Tall)
        error('Tall 必须是非空 table。');
    end
    if nargin < 2, cfg = struct(); end

    % 默认配置
    cfg = fill_default(cfg, 'sensor1', 'B1025');
    cfg = fill_default(cfg, 'sensor2', 'B454');
    cfg = fill_default(cfg, 'ch1', 1);
    cfg = fill_default(cfg, 'ch2', 2);
    cfg = fill_default(cfg, 'sensorMarkers', {'o','^'});
    cfg = fill_default(cfg, 'sensorMarkerSize', [34 42]);
    cfg = fill_default(cfg, 'sensorLineWidth', 1.0);
    cfg = fill_default(cfg, 'showPairedLine', false);
    cfg = fill_default(cfg, 'pairLineWidth', 0.5);
    cfg = fill_default(cfg, 'pairLineLightness', 0.35);

    if ~ismember('distance_mm', Tall.Properties.VariableNames)
        error('Tall 中缺少 distance_mm 列。');
    end

    featAliasMap = get_feature_alias_map();
    canonicalList = fieldnames(featAliasMap);
    availFeat = {};
    pairCols = struct();
    for i = 1:numel(canonicalList)
        feat = canonicalList{i};
        [c1, c2] = resolve_pair_cols(Tall, feat, featAliasMap);
        if ~isempty(c1) && ~isempty(c2)
            availFeat{end+1,1} = feat; %#ok<AGROW>
            pairCols.(feat) = struct('c1', c1, 'c2', c2);
        end
    end

    if isempty(availFeat)
        error('Tall 中未识别到可用的成对特征列（S1_/S2_）。');
    end

    defaultX = 'PF_kHz';
    defaultY = 'FC_kHz';
    if ~ismember(defaultX, availFeat), defaultX = availFeat{1}; end
    if ~ismember(defaultY, availFeat), defaultY = availFeat{min(2,numel(availFeat))}; end

    dList = unique(Tall.distance_mm(isfinite(Tall.distance_mm)))';
    if isempty(dList), dList = NaN; end
    cmap = lines(max(3, numel(dList)));

    % ===== UI =====
    fig = uifigure('Name', 'AE Feature Space Browser', 'Position', [80 60 1450 860]);

    g = uigridlayout(fig, [1 2]);
    g.ColumnWidth = {340, '1x'};
    g.RowHeight = {'1x'};
    g.Padding = [8 8 8 8];
    g.ColumnSpacing = 10;

    pnl = uipanel(g, 'Title', 'Controls');
    pnl.Layout.Row = 1; pnl.Layout.Column = 1;

    gc = uigridlayout(pnl, [18 2]);
    gc.RowHeight = {24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,'1x',28,28};
    gc.ColumnWidth = {110,'1x'};
    gc.Padding = [8 8 8 8];
    gc.RowSpacing = 6;
    gc.ColumnSpacing = 6;

    % 特征选择
    uilabel(gc, 'Text', 'X feature');
    ddX = uidropdown(gc, 'Items', availFeat, 'Value', defaultX);

    uilabel(gc, 'Text', 'Y feature');
    ddY = uidropdown(gc, 'Items', availFeat, 'Value', defaultY);

    % 坐标范围
    uilabel(gc, 'Text', 'X min');
    edXmin = uieditfield(gc, 'text', 'Value', '');

    uilabel(gc, 'Text', 'X max');
    edXmax = uieditfield(gc, 'text', 'Value', '');

    uilabel(gc, 'Text', 'Y min');
    edYmin = uieditfield(gc, 'text', 'Value', '');

    uilabel(gc, 'Text', 'Y max');
    edYmax = uieditfield(gc, 'text', 'Value', '');

    % 坐标类型
    uilabel(gc, 'Text', 'X scale');
    ddXscale = uidropdown(gc, 'Items', {'linear','log'}, 'Value', auto_scale_from_feat(defaultX));

    uilabel(gc, 'Text', 'Y scale');
    ddYscale = uidropdown(gc, 'Items', {'linear','log'}, 'Value', auto_scale_from_feat(defaultY));

    % 显示选项
    uilabel(gc, 'Text', 'Show paired line');
    cbLine = uicheckbox(gc, 'Value', logical(cfg.showPairedLine), 'Text', '');

    uilabel(gc, 'Text', 'Line width');
    efLineW = uieditfield(gc, 'numeric', 'Value', cfg.pairLineWidth, 'Limits',[0 10]);

    uilabel(gc, 'Text', 'Line lightness');
    efLight = uieditfield(gc, 'numeric', 'Value', cfg.pairLineLightness, 'Limits',[0 1]);

    btnAuto = uibutton(gc, 'push', 'Text', 'Auto Range');
    btnAuto.Layout.Column = [1 2];

    btnRefresh = uibutton(gc, 'push', 'Text', 'Refresh Plot');
    btnRefresh.Layout.Column = [1 2];

    % 距离勾选
    pDist = uipanel(gc, 'Title', 'Distances (checkbox)');
    pDist.Layout.Row = 16;
    pDist.Layout.Column = [1 2];
    pDist.Scrollable = 'on';

    gd = uigridlayout(pDist, [max(1,numel(dList)), 1]);
    gd.RowHeight = repmat({22}, 1, max(1,numel(dList)));
    gd.ColumnWidth = {'1x'};
    gd.Padding = [6 6 6 6];
    gd.RowSpacing = 4;

    distChecks = gobjects(numel(dList),1);
    for i = 1:numel(dList)
        distChecks(i) = uicheckbox(gd, 'Text', sprintf('%g mm', dList(i)), 'Value', true);
    end

    btnAll = uibutton(gc, 'push', 'Text', 'Select All Distances');
    btnNone = uibutton(gc, 'push', 'Text', 'Clear All Distances');

    % 右侧绘图区
    pPlot = uipanel(g, 'Title', 'Feature-space plot');
    pPlot.Layout.Row = 1; pPlot.Layout.Column = 2;

    gp = uigridlayout(pPlot, [2 1]);
    gp.RowHeight = {'1x', 32};
    gp.ColumnWidth = {'1x'};
    gp.Padding = [8 8 8 8];

    ax = uiaxes(gp);
    ax.Layout.Row = 1;
    grid(ax, 'on');
    hold(ax, 'on');

    lblStatus = uilabel(gp, 'Text', 'Ready', 'HorizontalAlignment', 'left');
    lblStatus.Layout.Row = 2;

    % 回调
    ddX.ValueChangedFcn = @(~,~) onFeatureChange();
    ddY.ValueChangedFcn = @(~,~) onFeatureChange();
    ddXscale.ValueChangedFcn = @(~,~) refreshPlot();
    ddYscale.ValueChangedFcn = @(~,~) refreshPlot();
    cbLine.ValueChangedFcn = @(~,~) refreshPlot();
    efLineW.ValueChangedFcn = @(~,~) refreshPlot();
    efLight.ValueChangedFcn = @(~,~) refreshPlot();

    edXmin.ValueChangedFcn = @(~,~) refreshPlot();
    edXmax.ValueChangedFcn = @(~,~) refreshPlot();
    edYmin.ValueChangedFcn = @(~,~) refreshPlot();
    edYmax.ValueChangedFcn = @(~,~) refreshPlot();

    for i = 1:numel(distChecks)
        distChecks(i).ValueChangedFcn = @(~,~) refreshPlot();
    end

    btnRefresh.ButtonPushedFcn = @(~,~) refreshPlot();
    btnAuto.ButtonPushedFcn    = @(~,~) autoRange();
    btnAll.ButtonPushedFcn     = @(~,~) setAllDistanceChecks(true);
    btnNone.ButtonPushedFcn    = @(~,~) setAllDistanceChecks(false);

    refreshPlot();

    %% ---- nested callbacks ----
    function onFeatureChange()
        ddXscale.Value = auto_scale_from_feat(ddX.Value);
        ddYscale.Value = auto_scale_from_feat(ddY.Value);
        autoRange();
        refreshPlot();
    end

    function setAllDistanceChecks(tf)
        for k = 1:numel(distChecks)
            if isgraphics(distChecks(k)), distChecks(k).Value = tf; end
        end
        refreshPlot();
    end

    function autoRange()
        [Tx, featX, featY, ok] = getFilteredTableForCurrentSelection();
        if ~ok || isempty(Tx), return; end

        [xAll, yAll] = getCombinedPoints(Tx, featX, featY);
        xAll = xAll(isfinite(xAll));
        yAll = yAll(isfinite(yAll));

        if strcmpi(ddXscale.Value,'log'), xAll = xAll(xAll>0); end
        if strcmpi(ddYscale.Value,'log'), yAll = yAll(yAll>0); end
        if isempty(xAll) || isempty(yAll), return; end

        xlimAuto = padded_lims(xAll);
        ylimAuto = padded_lims(yAll);

        edXmin.Value = num2str(xlimAuto(1), '%.6g');
        edXmax.Value = num2str(xlimAuto(2), '%.6g');
        edYmin.Value = num2str(ylimAuto(1), '%.6g');
        edYmax.Value = num2str(ylimAuto(2), '%.6g');
    end

    function refreshPlot()
        cla(ax);
        hold(ax, 'on');
        grid(ax, 'on');
        ax.Box = 'on';

        [Tx, featX, featY, ok, msg] = getFilteredTableForCurrentSelection();
        if ~ok
            lblStatus.Text = msg;
            title(ax, 'No data');
            return;
        end

        c1x = pairCols.(featX).c1; c2x = pairCols.(featX).c2;
        c1y = pairCols.(featY).c1; c2y = pairCols.(featY).c2;

        try
            ax.XScale = ddXscale.Value;
            ax.YScale = ddYscale.Value;
        catch
            ax.XScale = 'linear'; ax.YScale = 'linear';
        end

        % log轴过滤非正
        m = true(height(Tx),1);
        if strcmpi(ax.XScale,'log'), m = m & Tx.(c1x)>0 & Tx.(c2x)>0; end
        if strcmpi(ax.YScale,'log'), m = m & Tx.(c1y)>0 & Tx.(c2y)>0; end
        Tx = Tx(m,:);

        if isempty(Tx)
            lblStatus.Text = '当前筛选下无有效点（可能 log 轴与非正值冲突）';
            title(ax, 'No valid points');
            return;
        end

        dShow = getSelectedDistances();
        dShow = dShow(isfinite(dShow));
        Tx = Tx(ismember(Tx.distance_mm, dShow), :);

        if isempty(Tx)
            lblStatus.Text = '当前勾选距离下无数据';
            title(ax, 'No data for selected distances');
            return;
        end

        dUnique = unique(Tx.distance_mm(:))';
        hDist = gobjects(numel(dUnique),1);

        if cbLine.Value
            lw = efLineW.Value;
            lightness = efLight.Value;
            for ii = 1:height(Tx)
                di = Tx.distance_mm(ii);
                ic = find(dList == di, 1, 'first');
                if isempty(ic), ic = 1; end
                col = cmap(ic,:);
                colLight = lighten_color(col, lightness);
                plot(ax, [Tx.(c1x)(ii), Tx.(c2x)(ii)], [Tx.(c1y)(ii), Tx.(c2y)(ii)], '-', ...
                    'Color', colLight, 'LineWidth', lw, 'HandleVisibility','off');
            end
        end

        for i = 1:numel(dUnique)
            di = dUnique(i);
            idx = (Tx.distance_mm == di);
            ic = find(dList == di, 1, 'first');
            if isempty(ic), ic = 1; end
            col = cmap(ic,:);

            scatter(ax, Tx.(c1x)(idx), Tx.(c1y)(idx), cfg.sensorMarkerSize(1), ...
                cfg.sensorMarkers{1}, 'MarkerEdgeColor', col, 'MarkerFaceColor', 'none', ...
                'LineWidth', cfg.sensorLineWidth, 'HandleVisibility','off');

            scatter(ax, Tx.(c2x)(idx), Tx.(c2y)(idx), cfg.sensorMarkerSize(2), ...
                cfg.sensorMarkers{2}, 'MarkerEdgeColor', col, 'MarkerFaceColor', 'none', ...
                'LineWidth', cfg.sensorLineWidth, 'HandleVisibility','off');

            hDist(i) = plot(ax, nan, nan, '-', 'Color', col, 'LineWidth', 2.2, ...
                'DisplayName', sprintf('%g mm (n=%d)', di, sum(idx)));
        end

        hS1 = plot(ax, nan, nan, cfg.sensorMarkers{1}, ...
            'LineStyle','none', 'MarkerSize',7, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','none', ...
            'DisplayName', sprintf('%s (Ch%d)', cfg.sensor1, cfg.ch1));

        hS2 = plot(ax, nan, nan, cfg.sensorMarkers{2}, ...
            'LineStyle','none', 'MarkerSize',8, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','none', ...
            'DisplayName', sprintf('%s (Ch%d)', cfg.sensor2, cfg.ch2));

        xlabel(ax, featX, 'Interpreter','none');
        ylabel(ax, featY, 'Interpreter','none');
        title(ax, sprintf('%s vs %s in feature space (%s vs %s)', ...
            cfg.sensor1, cfg.sensor2, featY, featX), 'Interpreter','none', 'FontWeight','bold');

        lgHandles = [hS1, hS2, hDist(:).'];
        lgHandles = lgHandles(isgraphics(lgHandles));
        legend(ax, lgHandles, 'Location','best', 'Interpreter','none');

        try_apply_limits(ax, edXmin.Value, edXmax.Value, 'x');
        try_apply_limits(ax, edYmin.Value, edYmax.Value, 'y');

        lblStatus.Text = sprintf('显示 %d 个事件（每个事件含两传感器点） | X=%s | Y=%s', ...
            height(Tx), featX, featY);
    end

    function [Tx, featX, featY, ok, msg] = getFilteredTableForCurrentSelection()
        featX = ddX.Value;
        featY = ddY.Value;
        ok = true; msg = 'OK';

        if ~isfield(pairCols, featX) || ~isfield(pairCols, featY)
            Tx = table();
            ok = false;
            msg = '特征列不存在';
            return;
        end

        c1x = pairCols.(featX).c1; c2x = pairCols.(featX).c2;
        c1y = pairCols.(featY).c1; c2y = pairCols.(featY).c2;

        needCols = {'distance_mm', c1x,c2x,c1y,c2y};
        if ismember('eventID', Tall.Properties.VariableNames)
            needCols = [{'eventID'}, needCols];
        end
        Tx = Tall(:, needCols);

        m = isfinite(Tx.(c1x)) & isfinite(Tx.(c2x)) & isfinite(Tx.(c1y)) & isfinite(Tx.(c2y));
        Tx = Tx(m,:);

        if isempty(Tx)
            ok = false;
            msg = '当前特征组合没有有效数据';
        end
    end

    function dsel = getSelectedDistances()
        dsel = [];
        for k = 1:numel(distChecks)
            if isgraphics(distChecks(k)) && distChecks(k).Value
                dsel(end+1) = dList(k); %#ok<AGROW>
            end
        end
    end

    function [xAll, yAll] = getCombinedPoints(Tx, featX, featY)
        c1x = pairCols.(featX).c1; c2x = pairCols.(featX).c2;
        c1y = pairCols.(featY).c1; c2y = pairCols.(featY).c2;
        xAll = [Tx.(c1x); Tx.(c2x)];
        yAll = [Tx.(c1y); Tx.(c2y)];
    end
end

function cfg = fill_default(cfg, fieldName, defaultVal)
    if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
        cfg.(fieldName) = defaultVal;
    end
end

function s = auto_scale_from_feat(featName)
    if contains(string(featName), "Energy", 'IgnoreCase', true)
        s = 'log';
    else
        s = 'linear';
    end
end

function lims = padded_lims(v)
    v = v(isfinite(v));
    if isempty(v), lims = [0 1]; return; end
    lo = min(v); hi = max(v);
    if lo == hi, lo = lo - 1; hi = hi + 1; end
    pad = 0.05 * (hi - lo);
    lims = [lo-pad, hi+pad];
end

function try_apply_limits(ax, smin, smax, whichAxis)
    xmin = str2double(strtrim(smin));
    xmax = str2double(strtrim(smax));
    if ~isfinite(xmin) || ~isfinite(xmax) || xmin >= xmax
        return;
    end
    try
        switch lower(whichAxis)
            case 'x'
                xlim(ax, [xmin xmax]);
            case 'y'
                ylim(ax, [xmin xmax]);
        end
    catch
    end
end