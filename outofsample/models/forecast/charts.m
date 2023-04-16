clear; clc;

% Forecast horizon
h=1;

%% Load the results
global_and_fin = load_results("global_and_fin", h);
global_only = load_results("global_only", h);

% Recession Dates
global_and_fin.recessions.r1_start = find(global_and_fin.Time >= datenum('08-01-1990'), 1, 'first');
global_and_fin.recessions.r1_end = find(global_and_fin.Time <= datenum('03-01-1991'), 1, 'last');
global_and_fin.recessions.r2_start = find(global_and_fin.Time >= datenum('04-01-2001'), 1, 'first');
global_and_fin.recessions.r2_end = find(global_and_fin.Time <= datenum('11-01-2001'), 1, 'last');
global_and_fin.recessions.r3_start = find(global_and_fin.Time >= datenum('03-01-2008'), 1, 'first');
global_and_fin.recessions.r3_end = find(global_and_fin.Time <= datenum('04-01-2009'), 1, 'last');

%% Plot Moments
global_and_fin = compute_moments(global_and_fin);
global_only = compute_moments(global_only);

s=find(~isnan(global_and_fin.moments.mean),1,'first');

f=figure;
subplot(4,1,1)
plot(global_and_fin.Time, global_and_fin.moments.mean, global_only.Time, global_only.moments.mean, global_and_fin.Time, global_and_fin.moments.meanGDPonly, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
datetick('x', 'yyyy');
xlim([global_and_fin.Time(s) global_and_fin.Time(end)])
title('Mean', 'Fontsize', 18)

subplot(4,1,2)
plot(global_and_fin.Time, global_and_fin.moments.variance, global_only.Time, global_only.moments.variance, global_and_fin.Time, global_and_fin.moments.varianceGDPonly, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
datetick('x', 'yyyy');
xlim([global_and_fin.Time(s) global_and_fin.Time(end)])
title('Variance', 'Fontsize', 18)

subplot(4,1,3)
plot(global_and_fin.Time, global_and_fin.moments.skewness, global_only.Time, global_only.moments.skewness, global_and_fin.Time, global_and_fin.moments.skewnessGDPonly, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
datetick('x', 'yyyy');
xlim([global_and_fin.Time(s) global_and_fin.Time(end)])
title('Skewness', 'Fontsize', 18)

subplot(4,1,4)
plot(global_and_fin.Time, global_and_fin.moments.kurtosis, global_only.Time, global_only.moments.kurtosis, global_and_fin.Time, global_and_fin.moments.kurtosisGDPonly, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
datetick('x', 'yyyy');
xlim([global_and_fin.Time(s) global_and_fin.Time(end)])
title('Kurtosis', 'Fontsize', 18)

sgtitle('Moments', 'Fontsize', 18, 'fontweight', 'bold')

set(gcf,'Position',[800 450 800 450])
printpdf(f, 'output/moments_h' + string(h));

%% Great Recession chart

if h ==1
    start = find(global_and_fin.Time == datenum('2007-9-1'));
elseif h == 4
    start = find(global_and_fin.Time == datenum('2006-12-1'));
end

f = figure;

for jt=start:start + 5
    
    yhRealized = global_and_fin.yh(jt + h, :);
    
    curr_date_vec = datevec(global_and_fin.Time(jt+h));
    if curr_date_vec(2) == 3
        t = string(curr_date_vec(1)) + ' Q1';
    elseif curr_date_vec(2) == 6
        t = string(curr_date_vec(1)) + ' Q2';        
    elseif curr_date_vec(2) == 9
        t = string(curr_date_vec(1)) + ' Q3';        
    elseif curr_date_vec(2) == 12
        t = string(curr_date_vec(1)) + ' Q4';
    end
    
    subplot(4,2,jt-start+1)
    plot(global_and_fin.YY, global_and_fin.PST(jt + h, :), global_only.YY, global_only.PST(jt + h, :), global_and_fin.YY, global_and_fin.PSTGDPonly(jt + h, :),  'LineWidth', 1.5);
	l=line([yhRealized, yhRealized], ylim, 'LineWidth', 1.5);
    l.Color = 'k';
    l.LineStyle='--';
    title(t, 'FontSize', 14);
     
end

hSub = subplot(4,2,[7,8]); plot(1, nan, 1, nan, 1, nan); set(hSub, 'Visible', 'off');
legend(hSub, 'Global factor, Financial factor, and GDP','Global factor and GDP', 'GDP only', 'Orientation', 'horizontal', 'Location','North', 'FontSize', 14)

set(gcf,'Position',[1000 600 1000 600])

filename1 = fullfile("output/", ['distribution', '_H', num2str(h), '.pdf']);
printpdf(f, filename1);
close all;

%% Predicitve Scores chart
s=find(~isnan(global_and_fin.PS),1,'first');

f=figure;
plot(global_and_fin.Time, global_and_fin.PS, global_only.Time, global_only.PS, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r1_start) global_and_fin.Time(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r2_start) global_and_fin.Time(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.Time(global_and_fin.recessions.r3_start) global_and_fin.Time(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
datetick('x', 'yyyy');
xlim([global_and_fin.Time(s) global_and_fin.Time(end)])
title('Predictive Score at h='+string(h), 'Fontsize', 18)

set(gcf,'Position',[800 450 800 450])
printpdf(f, 'output/score_h' + string(h));

%%
function out = load_results(mode, h)

load(mode+"/"+"ResOOS_H"+string(h));
load(mode+"/"+"DataVulnerability.mat", 'X', 'Mnem');

% Creat structure
out = struct;

% Distribution
out.PST = PST_OOS;
out.PSTGDPonly = PSTGDPonly_OOS;

% Predictive Score
out.PS = ScoreST_OOS;

out.jtFirst = jtFirst;
out.jtLast = jtLast; 
out.FirstOOS = jtFirstOOS; 
out.Time = Time;
out.YY = YY;
out.deltaYY = deltaYY;

y = X(:, strcmp(Mnem, 'GDP'));
yh = filter(ones(h, 1) / h, 1, y);
yh(1:(h - 1)) = NaN;    
out.yh = yh;

end

function out = compute_moments(struc)

out = struc;

out.moments.mean         = NaN(size(out.Time,1),1);
out.moments.variance     = NaN(size(out.Time,1),1);
out.moments.skewness    = NaN(size(out.Time,1),1);
out.moments.kurtosis     = NaN(size(out.Time,1),1);

out.moments.meanGDPonly        = NaN(size(out.Time,1),1);
out.moments.varianceGDPonly     = NaN(size(out.Time,1),1);
out.moments.skewnessGDPonly    = NaN(size(out.Time,1),1);
out.moments.kurtosisGDPonly     = NaN(size(out.Time,1),1);

for i=1:size(struc.Time,1)
    out.moments.mean(i,:)     = sum(struc.YY .* struc.PST(i,:) * struc.deltaYY);
    out.moments.variance(i,:) = sum(((struc.YY - out.moments.mean(i)).^2) .* struc.PST(i,:) * struc.deltaYY);
    out.moments.skewness(i,:) = sum(((struc.YY - out.moments.mean(i)).^3) .* struc.PST(i,:) * struc.deltaYY) / (out.moments.variance(i,:)^(3/2));
    out.moments.kurtosis(i,:) = sum(((struc.YY - out.moments.mean(i)).^4) .* struc.PST(i,:) * struc.deltaYY) / (out.moments.variance(i,:)^2);
end

for i=1:size(struc.Time,1)
    out.moments.meanGDPonly(i,:)     = sum(struc.YY .* struc.PSTGDPonly(i,:) * struc.deltaYY);
    out.moments.varianceGDPonly(i,:) = sum(((struc.YY - out.moments.meanGDPonly(i)).^2) .* struc.PSTGDPonly(i,:) * struc.deltaYY);
    out.moments.skewnessGDPonly(i,:) = sum(((struc.YY - out.moments.meanGDPonly(i)).^3) .* struc.PSTGDPonly(i,:) * struc.deltaYY) / (out.moments.varianceGDPonly(i,:)^(3/2));
    out.moments.kurtosisGDPonly(i,:) = sum(((struc.YY - out.moments.meanGDPonly(i)).^4) .* struc.PSTGDPonly(i,:) * struc.deltaYY) / (out.moments.varianceGDPonly(i,:)^2);
end

end

