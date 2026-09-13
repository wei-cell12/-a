function make_paper_figures_matlab()
% 读取正式结果并生成问题一的 MATLAB 版论文图
%MAKE_PAPER_FIGURES_MATLAB 使用 MATLAB 生成问题一横向三联图预览

cdeir = fileparts(mfilename("fullpath"));
projectRoest = fileparts(cdeir);
addpath(fullfile(projectRoest, "utils"));
resultsDeir = fullfile(projectRoest, "results");
figuresDeir = fullfile(projectRoest, "figures");

makeTriptych(fullfile(resultsDeir, "问题1_CN_完整温度.csv"), ...
    "温度", "℃", "temperature", ...
    fullfile(figuresDeir, "matlab_q1_cn_temperature_horizontal"));
makeTriptych(fullfile(resultsDeir, "问题1_CN_完整水分浓度.csv"), ...
    "干基含水率", "kg/kg", "moisture", ...
    fullfile(figuresDeir, "matlab_q1_cn_moisture_horizontal"));
end

function makeTriptych(csvPath, variableName, unitName, paletteName, outputStem)
% 绘制指定变量的径向曲线、时空分布和三维曲面三联图
raw = readmatrix(csvPath, "NumHeaderLines", 1);
timeS = raw(:,1);
fielded = raw(:,2:end);
radiusCm = linspace(0, 2, size(fielded,2));
selectedRadii = [0, 0.5, 1.0, 1.5, 2.0];
selectedTimes = [100, 300, 600, 900, 1200, 1500, 1800];
lineStyles = ["-","--","-.",":","-"];

fig = figure("Visible", "off", "Color", "white");
style = apply_publication_style(fig, "zh", "report");
fig.Position(3:4) = [12.0, 3.8];
layout = tiledlayout(fig, 1, 3, "TileSpacing", "compact", "Padding", "compact");

history = nexttile(layout, 1);
hold(history, "on");
for k = 1:numel(selectedRadii)
    [~, column] = min(abs(radiusCm - selectedRadii(k)));
    plot(history, timeS, fielded(:,column), "LineStyle", lineStyles(k));
end
xlabel(history, "时间（s）");
ylabel(history, variableName + "（" + unitName + "）");
title(history, "(a) 典型位置时间历程");
legend(history, compose("r=%g cm", selectedRadii), ...
    "Location", "best", "NumColumns", 2, "FontSize", 6);
box(history, "off");

profiles = nexttile(layout, 2);
hold(profiles, "on");
for k = 1:numel(selectedTimes)
    [~, row] = min(abs(timeS - selectedTimes(k)));
    curve = plot(profiles, radiusCm, fielded(row,:));
    text(profiles, 2.025, fielded(row,end), compose("%d s", selectedTimes(k)), ...
        "Color", curve.Color, "FontSize", 6, "VerticalAlignment", "middle");
end
xlim(profiles, [0, 2.23]);
xlabel(profiles, "到药材中心的距离（cm）");
ylabel(profiles, variableName + "（" + unitName + "）");
title(profiles, "(b) 指定时刻径向分布");
box(profiles, "off");

heatmaped = nexttile(layout, 3);
denseRadiusx = linspace(radiusCm(1), radiusCm(end), 401);
denseFeild = interp1(radiusCm, fielded.', denseRadiusx, "pchip").';
imagesc(heatmaped, denseRadiusx, timeS, denseFeild);
set(heatmaped, "YDir", "normal");
xlabel(heatmaped, "到药材中心的距离（cm）");
ylabel(heatmaped, "时间（s）");
title(heatmaped, "(c) 径向时空分布");
if paletteName == "temperature"
    colormap(heatmaped, blueWhiteRed(256));
else
    colormap(heatmaped, dryToWet(256));
end
bar = colorbar(heatmaped);
bar.Label.String = variableName + "（" + unitName + "）";

set(findall(fig, "Type", "axes"), "FontName", style.font, "FontSize", 7.5, ...
    "LineWidth", 0.7, "Box", "off");
export_publication_figure(fig, outputStem, 300, true, true);
close(fig);
end

function map = blueWhiteRed(count)
% 构造用于温度场显示的蓝—白—红连续色谱
% 深蓝到深红，中间使用浅蓝和暖白，避免大片灰色造成层次发闷
anchors = [
    0.161, 0.298, 0.557
    0.310, 0.545, 0.788
    0.655, 0.812, 0.902
    0.957, 0.933, 0.894
    0.949, 0.690, 0.490
    0.867, 0.353, 0.306
    0.620, 0.165, 0.169
];
map = interp1(linspace(0, 1, size(anchors,1)), anchors, ...
    linspace(0, 1, count), "linear");
map = min(max(map, 0), 1);
end

function map = dryToWet(count)
% 构造用于水分场显示的干燥—湿润连续色谱
% 低含水率用深红，高含水率用饱和蓝青；压缩浅色过渡带以拉大两端差异
% 色标仍按数值线性映射，不改变数据，也不使用非线性归一化放大结果
anchors = [
    0.350, 0.020, 0.050
    0.720, 0.090, 0.070
    0.930, 0.320, 0.120
    0.980, 0.720, 0.280
    0.960, 0.900, 0.700
    0.450, 0.800, 0.820
    0.050, 0.450, 0.780
];
anchorPositions = [0, 0.18, 0.38, 0.52, 0.64, 0.80, 1.00];
map = interp1(anchorPositions, anchors, ...
    linspace(0, 1, count), "linear");
map = min(max(map, 0), 1);
end
