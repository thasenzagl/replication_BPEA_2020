clear; clc;

%% Load the results
global_and_fin = load_results("global_and_fin");
global_only = load_results("global_only");
nonfin_only = load_results("nonfin_only");

% Recession Dates
global_and_fin.recessions.r1_start = find(global_and_fin.name_vintage >= datenum('08-01-1990'), 1, 'first');
global_and_fin.recessions.r1_end = find(global_and_fin.name_vintage <= datenum('03-01-1991'), 1, 'last');
global_and_fin.recessions.r2_start = find(global_and_fin.name_vintage >= datenum('04-01-2001'), 1, 'first');
global_and_fin.recessions.r2_end = find(global_and_fin.name_vintage <= datenum('11-01-2001'), 1, 'last');
global_and_fin.recessions.r3_start = find(global_and_fin.name_vintage >= datenum('01-01-2008'), 1, 'first');
global_and_fin.recessions.r3_end = find(global_and_fin.name_vintage <= datenum('04-01-2009'), 1, 'last');

%% Plot Moments
global_and_fin = compute_moments(global_and_fin);
global_only = compute_moments(global_only);
nonfin_only = compute_moments(nonfin_only);

f=figure;
subplot(4,1,1)
plot(global_and_fin.name_vintage, global_and_fin.moments.mean.nc, global_only.name_vintage, global_only.moments.mean.nc, nonfin_only.name_vintage, nonfin_only.moments.mean.nc, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
datetick('x', 'yyyy');
title('Mean', 'Fontsize', 18)

subplot(4,1,2)
plot(global_and_fin.name_vintage, global_and_fin.moments.variance.nc, global_only.name_vintage, global_only.moments.variance.nc, nonfin_only.name_vintage, nonfin_only.moments.variance.nc, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
title('Variance', 'Fontsize', 18)

subplot(4,1,3)
plot(global_and_fin.name_vintage, global_and_fin.moments.skewness.nc, global_only.name_vintage, global_only.moments.skewness.nc, nonfin_only.name_vintage, nonfin_only.moments.skewness.nc, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
title('Skewness', 'Fontsize', 18)

subplot(4,1,4)
plot(global_and_fin.name_vintage, global_and_fin.moments.kurtosis.nc, global_only.name_vintage, global_only.moments.kurtosis.nc, nonfin_only.name_vintage, nonfin_only.moments.kurtosis.nc, 'Linewidth', 1.5)
hold on
yl=ylim;
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
title('Kurtosis', 'Fontsize', 18)

sgtitle('Moments of the Nowcast', 'Fontsize', 18, 'fontweight', 'bold')

set(gcf,'Position',[800 450 800 450])
printpdf(f, 'output/charts/moments');

%% Plot Factors
f=figure;
plot(datenum(global_and_fin.datesQ), global_and_fin.factors.global(:,end), 'Linewidth', 1.5)
hold on
plot(datenum(global_and_fin.datesQ), global_and_fin.factors.fin(:,end), 'Linewidth', 1.5)
hold on
%plot(datenum(nonfin_only.datesQ), nonfin_only.factors.global(:,end), 'Linewidth', 1.5)
%hold on
yl=ylim;
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
datetick('x', 'yyyy');

title('Global and Non-financial Factors', 'Fontsize', 18)

legend({'Global Factor','Financial Factor'}, 'Location', 'SouthWest', 'FontSize', 14)
set(gcf,'Position',[800 450 800 450])
printpdf(f, 'output/charts/factors');

%% Shortfall and Longrise
% Shortfall
f=figure;
subplot(2,1,1);
plot(global_and_fin.name_vintage, global_and_fin.SF.nc, global_only.name_vintage, global_only.SF.nc, nonfin_only.name_vintage, nonfin_only.SF.nc, 'Linewidth', 1.5);
hold on;
yl=ylim;
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
ylabel('Shortfall', 'FontSize', 14)
title('Shortfall', 'Fontsize', 18)

% Longrise
subplot(2,1,2);
plot(global_and_fin.name_vintage, global_and_fin.LR.nc, global_only.name_vintage, global_only.LR.nc, nonfin_only.name_vintage, nonfin_only.LR.nc, 'Linewidth', 1.5);
hold on;
yl=ylim;
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r1_start) global_and_fin.name_vintage(global_and_fin.recessions.r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r2_start) global_and_fin.name_vintage(global_and_fin.recessions.r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([global_and_fin.name_vintage(global_and_fin.recessions.r3_start) global_and_fin.name_vintage(global_and_fin.recessions.r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
ylabel('Longrise', 'FontSize', 14)
title('Longrise', 'Fontsize', 18)

legend({'Global Factor','Global and Financial Factors','Non-financial Factor'}, 'Location', 'NorthEast', 'FontSize', 14)
datetick('x', 'yyyy');
set(gcf,'Position',[800 450 800 450])
printpdf(f, 'output/charts/sf_lr');


%% Plot RMSFE and Predictive Score
global_and_fin = compute_fcst_errors(global_and_fin);
global_only = compute_fcst_errors(global_only);
nonfin_only = compute_fcst_errors(nonfin_only);

% RMSFE
f=figure;
h1=subplot(2,1,1);
plot(global_and_fin.distance(20:end), global_and_fin.meanSE(20:end), global_only.distance(20:end), global_only.meanSE(20:end), nonfin_only.distance(20:end), nonfin_only.meanSE(20:end), 'Linewidth', 1.5);
xlabel('Distance to Release (days)', 'FontSize', 14)
ylabel('RMSE', 'FontSize', 14)
title('RMSE', 'Fontsize', 18)

ylim([min(global_and_fin.meanSE(20:end))-1,max(global_and_fin.meanSE(20:end))+1])
set(h1, 'Xdir', 'reverse');

% Predictive Score
h2=subplot(2,1,2);
plot(global_and_fin.distance(20:end), global_and_fin.meanPS(20:end), global_only.distance(20:end), global_only.meanPS(20:end), nonfin_only.distance(20:end), nonfin_only.meanPS(20:end), 'Linewidth', 1.5);
xlabel('Distance to Release (days)', 'FontSize', 14)
ylabel('Predicitve Score', 'FontSize', 14)
title('Predicitve Score', 'Fontsize', 18)
legend({'Global Factor','Global and Financial Factors','Non-financial Factor'}, 'Location', 'NorthWest', 'FontSize', 14)

ylim([min(global_and_fin.meanPS(20:end))-0.01,max(global_and_fin.meanPS(20:end))+0.01])
set(h2, 'Xdir', 'reverse');
set(gcf,'Position',[800 450 800 450])
printpdf(f, 'output/charts/rmse_ps');

close all;

%%
function out = load_results(mode)

load("output/matlab_output/output_" + mode, "PST_bc", "PST_nc", "PST_fc", ...
    "PS_bc", "PS_nc", "PS_fc", "OT_bc", "OT_nc", "OT_fc","SF_bc", "SF_nc",...
    "SF_fc", "LR_bc", "LR_nc", "LR_fc", "YY", "deltaYY", "name_vintage", ...
    "datesQ", "factor_fin", "factor_global", "y", "ReleasesMat", "Period");

% Creat structure
out = struct;

% Distribution
out.PST.bc = PST_bc;
out.PST.nc = PST_nc;
out.PST.fc = PST_fc;

% Predictive Score
out.PS.bc = PS_bc;
out.PS.nc = PS_nc;
out.PS.fc = PS_fc;

% Outturn
out.OT.bc = OT_bc;
out.OT.nc = OT_nc;
out.OT.fc = OT_fc;

% Shortfall
out.SF.bc = SF_bc;
out.SF.nc = SF_nc;
out.SF.fc = SF_fc;

% Longrise
out.LR.bc = LR_bc;
out.LR.nc = LR_nc;
out.LR.fc = LR_fc;

out.YY = YY;
out.deltaYY = deltaYY;
out.name_vintage = name_vintage;
out.ReleasesMat = ReleasesMat;
out.Period = Period(:,1:2);

% Factors
out.datesQ = datenum([datesQ, ones(size(datesQ,1),1)]);
out.gdp = y;
out.factors.global = factor_global;
out.factors.fin = factor_fin;

end

function out = compute_moments(struc)

out = struc;

out.moments.mean.fc         = NaN(size(out.name_vintage,1),1);
out.moments.variance.fc     = NaN(size(out.name_vintage,1),1);
out.moments.skewness.fc     = NaN(size(out.name_vintage,1),1);
out.moments.kurtosis.fc     = NaN(size(out.name_vintage,1),1);

out.moments.mean.nc         = NaN(size(out.name_vintage,1),1);
out.moments.variance.nc     = NaN(size(out.name_vintage,1),1);
out.moments.skewness.nc     = NaN(size(out.name_vintage,1),1);
out.moments.kurtosis.nc     = NaN(size(out.name_vintage,1),1);

out.moments.mean.bc         = NaN(size(out.name_vintage,1),1);
out.moments.variance.bc     = NaN(size(out.name_vintage,1),1);
out.moments.skewness.bc     = NaN(size(out.name_vintage,1),1);
out.moments.kurtosis.bc     = NaN(size(out.name_vintage,1),1);

for i=1:size(struc.name_vintage,1)
    out.moments.mean.bc(i,:)     = sum(struc.YY .* struc.PST.bc(i,:) * struc.deltaYY);
    out.moments.variance.bc(i,:) = sum(((struc.YY - out.moments.mean.bc(i)).^2) .* struc.PST.bc(i,:) * struc.deltaYY);
    out.moments.skewness.bc(i,:) = sum(((struc.YY - out.moments.mean.bc(i)).^3) .* struc.PST.bc(i,:) * struc.deltaYY) / (out.moments.variance.bc(i,:)^(3/2));
    out.moments.kurtosis.bc(i,:) = sum(((struc.YY - out.moments.mean.bc(i)).^4) .* struc.PST.bc(i,:) * struc.deltaYY) / (out.moments.variance.bc(i,:)^2);
end

for i=1:size(struc.name_vintage,1)
    out.moments.mean.nc(i,:)     = sum(struc.YY .* struc.PST.nc(i,:) * struc.deltaYY);
    out.moments.variance.nc(i,:) = sum(((struc.YY - out.moments.mean.nc(i)).^2) .* struc.PST.nc(i,:) * struc.deltaYY);
    out.moments.skewness.nc(i,:) = sum(((struc.YY - out.moments.mean.nc(i)).^3) .* struc.PST.nc(i,:) * struc.deltaYY) / (out.moments.variance.nc(i,:)^(3/2));
    out.moments.kurtosis.nc(i,:) = sum(((struc.YY - out.moments.mean.nc(i)).^4) .* struc.PST.nc(i,:) * struc.deltaYY) / (out.moments.variance.nc(i,:)^2);
end

for i=1:size(struc.name_vintage,1)
    out.moments.mean.fc(i,:)     = sum(struc.YY .* struc.PST.fc(i,:) * struc.deltaYY);
    out.moments.variance.fc(i,:) = sum(((struc.YY - out.moments.mean.fc(i)).^2) .* struc.PST.fc(i,:) * struc.deltaYY);
    out.moments.skewness.fc(i,:) = sum(((struc.YY - out.moments.mean.fc(i)).^3) .* struc.PST.fc(i,:) * struc.deltaYY) / (out.moments.variance.fc(i,:)^(3/2));
    out.moments.kurtosis.fc(i,:) = sum(((struc.YY - out.moments.mean.fc(i)).^4) .* struc.PST.fc(i,:) * struc.deltaYY) / (out.moments.variance.fc(i,:)^2);
end

end

function out = compute_fcst_errors(struc)

out = struc;

target_releases = struc.ReleasesMat(end,:);
distance = 199:-1:0;
PS = nan(size(distance,2), size(struc.Period,1));
SE = nan(size(distance,2), size(struc.Period,1));

% Find distance for each release
for t = 1:size(struc.name_vintage,1)
    
    %% Find reference quarter
    temp = datevec(struc.name_vintage(t));
    temp = temp(1:2);
    
    % Turn months 1,2->3, 4,5->6, 7,8->9, 10,11->12
    if  mod(temp(2),3) == 1
        temp(2)= temp(2)+2;
    elseif mod(temp(2),3) == 2
        temp(2)= temp(2)+1;
    end
    
    % Index of DatesV where the date is the same as P.Qnews
    i_ref = find(struc.Period(:,1) == temp(1) & struc.Period(:,2) == temp(2));
    
    for i=-3:3:3
        if i_ref+i <= size(struc.ReleasesMat,2)
            dist = datenum(target_releases(i_ref+i)) - struc.name_vintage(t);
            
            idx = find(distance == dist);
            
            if i == -3
                PS(idx, i_ref+i) = struc.PS.bc(t);
                SE(idx, i_ref+i) = (struc.moments.mean.bc(t) - struc.OT.bc(t))^2;
            elseif i == 0
                PS(idx, i_ref+i) = struc.PS.nc(t);
                SE(idx, i_ref+i) = (struc.moments.mean.nc(t) - struc.OT.nc(t))^2;
            else
                PS(idx, i_ref+i) = struc.PS.fc(t);
                SE(idx, i_ref+i) = (struc.moments.mean.fc(t) - struc.OT.fc(t))^2;
            end
        end
    end

end
    %  Replace NaNs and average across quarters for SE & PS
    idx = sum(isnan(PS),1) ~= size(PS,1);
    PS = PS(:,idx);
    SE = SE(:,idx);
    
    for j = 1:size(PS,2)
        for i = 2:size(PS,1)
            if isnan(PS(i,j))
                PS(i,j) = PS(i-1,j);
                SE(i,j) = SE(i-1,j);
            end
        end
    end
    
    out.meanPS = nanmean(PS,2);
    out.meanSE = nanmean(SE,2);
    out.distance = distance;

end
