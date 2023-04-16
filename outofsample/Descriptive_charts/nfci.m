clear; clc;

%% Load the data
csv_in = strcat('data/nfci.csv');
data_struct = importdata(csv_in,',');
dates = data_struct.data(:,1)+693960;
data = data_struct.data(:,2:end);

% Recession Dates
r1_start = find(dates >= datenum('01-01-1974'), 1, 'first');
r1_end = find(dates <= datenum('01-01-1975'), 1, 'last');

r2_start = find(dates >= datenum('02-01-1980'), 1, 'first');
r2_end = find(dates <= datenum('07-01-1980'), 1, 'last');

r3_start = find(dates >= datenum('09-01-1981'), 1, 'first');
r3_end = find(dates <= datenum('11-01-1982'), 1, 'last');

r4_start = find(dates >= datenum('08-01-1990'), 1, 'first');
r4_end = find(dates <= datenum('03-01-1991'), 1, 'last');

r5_start = find(dates >= datenum('04-01-2001'), 1, 'first');
r5_end = find(dates <= datenum('11-01-2001'), 1, 'last');

r6_start = find(dates >= datenum('01-01-2008'), 1, 'first');
r6_end = find(dates <= datenum('06-01-2009'), 1, 'last');

%%
f=figure; 

plot(dates, data(:,1), dates, data(:,2), 'Linewidth', 1.5)
hold on
yl=ylim;

% Add Recession bars
area([dates(r1_start) dates(r1_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r2_start) dates(r2_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r3_start) dates(r3_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r4_start) dates(r4_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r5_start) dates(r5_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r6_start) dates(r6_end)], [yl(2) yl(2)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');

area([dates(r1_start) dates(r1_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r2_start) dates(r2_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r3_start) dates(r3_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r4_start) dates(r4_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r5_start) dates(r5_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');
area([dates(r6_start) dates(r6_end)], [yl(1) yl(1)], 'FaceAlpha', 0.1, 'EdgeColor','none', 'FaceColor', 'k');

datetick('x', 'yyyy');
title('NFCI and Real Activity Indicator', 'fontsize',14);
legend('NFCI','Real Activity Indicator')
xlim([dates(1) dates(end)]);
set(gcf,'Position',[1000 500 1000 500])

printpdf(f, 'output/nfci');

close all;