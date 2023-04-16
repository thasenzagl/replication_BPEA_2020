clear; clc; close all

mode = 'nonfin_only';
yhRealized = -4.8;

%% Load results
load("output/matlab_output/output_" + mode);

%% Plot

f=figure;
subplot(3,1,1)
plot(YY,PST_t0(1,:), YY,PST_t0(2,:), YY,PST_tm1(3,:), 'Linewidth', 1.5);
xlim([-10 10]);
l=line([yhRealized, yhRealized], ylim, 'LineWidth', 1.5);
l.Color = 'k';
l.LineStyle='--';
grid on;
title('Predictive Distributions for 2020 Q1', 'Fontsize', 18);
legend({datestr(name_vintage(1)),datestr(name_vintage(2)),datestr(name_vintage(3)), 'Advance Estimate'}, 'Location', 'Southwest', 'FontSize', 12);
    
subplot(3,1,2)
plot(YY,PST_t1(1,:), YY,PST_t1(2,:), YY,PST_t0(3,:), 'Linewidth', 1.5);
xlim([-10 10])
grid on;
title('Predictive Distributions for 2020 Q2', 'Fontsize', 18);

subplot(3,1,3)
plot(YY,PST_t4(1,:), YY,PST_t4(2,:), YY,PST_t3(3,:), 'Linewidth', 1.5);
xlim([-10 10])
grid on;
title('Predictive Distributions for 2021 Q1', 'Fontsize', 18);

if strcmp(mode,'global_only')
    sgtitle('Global Factor Only');
elseif strcmp(mode,'global_and_fin')
        sgtitle('Global and Financial Factor');
elseif strcmp(mode,'nonfin_only')
        sgtitle('Nonfinancial Factor Only');
end    

set(gcf,'Position',[800 450 800 450]);
printpdf(f, "output/charts/nowcast_covid_"+mode);

%% Factor Plot
clear, clc; 

% Load results
load("output/matlab_output/output_nonfin_only");
factor_nonfin = factor_global;
load("output/matlab_output/output_global_and_fin");

start = '01-Mar-2017';
last = '01-Sep-2020';
fcst_start = '01-Dec-2019';
dates = datenum([datesQ ones(size(datesQ,1),1)]);

idx_start = find(dates == datenum(start));
idx_last = find(dates == datenum(last));
idx_fcst = find(dates == datenum(fcst_start));

f=figure;
subplot(3,1,1)
hold all;
ax = gca;
plot(dates(idx_start:idx_fcst), factor_global(idx_start:idx_fcst,1), 'k', 'Linewidth', 1.5)
plot(dates(idx_start:idx_fcst), factor_global(idx_start:idx_fcst,2), 'b', 'Linewidth', 1.5)
plot(dates(idx_start:idx_fcst), factor_global(idx_start:idx_fcst,3), 'r', 'Linewidth', 1.5)

plot(dates(idx_fcst:idx_last), factor_global(idx_fcst:idx_last,1), 'k--', 'Linewidth', 1.5);
plot(dates(idx_fcst:idx_last), factor_global(idx_fcst:idx_last,2), 'b--', 'Linewidth', 1.5);
plot(dates(idx_fcst:idx_last), factor_global(idx_fcst:idx_last,3), 'r--', 'Linewidth', 1.5);

hold off;

legend(datestr(name_vintage), 'Location', 'Southwest', 'FontSize', 14);
title('Global factor', 'Fontsize', 18);
grid on;
set(gca, 'XTick', dates(idx_start:idx_last));
datetick('x', 'mmmyy', 'keepticks');

subplot(3,1,2)
hold all;
ax = gca;
plot(dates(idx_start:idx_fcst), factor_fin(idx_start:idx_fcst,1), 'k', 'Linewidth', 1.5)
plot(dates(idx_start:idx_fcst), factor_fin(idx_start:idx_fcst,2), 'b', 'Linewidth', 1.5)
plot(dates(idx_start:idx_fcst), factor_fin(idx_start:idx_fcst,3), 'r', 'Linewidth', 1.5)

plot(dates(idx_fcst:idx_last), factor_fin(idx_fcst:idx_last,1), 'k--', 'Linewidth', 1.5);
plot(dates(idx_fcst:idx_last), factor_fin(idx_fcst:idx_last,2), 'b--', 'Linewidth', 1.5);
plot(dates(idx_fcst:idx_last), factor_fin(idx_fcst:idx_last,3), 'r--', 'Linewidth', 1.5);

title('Financial factor', 'Fontsize', 18);
grid on;
set(gca, 'XTick', dates(idx_start:idx_last));
datetick('x', 'mmmyy', 'keepticks');

subplot(3,1,3)
hold all;
ax = gca;
plot(dates(idx_start:idx_fcst), factor_nonfin(idx_start:idx_fcst,1), 'k', 'Linewidth', 1.5)
plot(dates(idx_start:idx_fcst), factor_nonfin(idx_start:idx_fcst,2), 'b', 'Linewidth', 1.5)
plot(dates(idx_start:idx_fcst), factor_nonfin(idx_start:idx_fcst,3), 'r', 'Linewidth', 1.5)

plot(dates(idx_fcst:idx_last), factor_nonfin(idx_fcst:idx_last,1), 'k--', 'Linewidth', 1.5);
plot(dates(idx_fcst:idx_last), factor_nonfin(idx_fcst:idx_last,2), 'b--', 'Linewidth', 1.5);
plot(dates(idx_fcst:idx_last), factor_nonfin(idx_fcst:idx_last,3), 'r--', 'Linewidth', 1.5);

title('Nonfinancial factor', 'Fontsize', 18);
grid on;
set(gca, 'XTick', dates(idx_start:idx_last));
datetick('x', 'mmmyy', 'keepticks');

set(gcf,'Position',[800 450 800 450]);
printpdf(f, "output/charts/factors_covid");
