clear, clc

%% Load the data
addpath(genpath('../QR_library/'));

csv_in = strcat('data/histogram.csv');

% load csv
data_struct = importdata(csv_in,',');

% get series names
series = data_struct.textdata(1,2:end)';

%get reference dates
dates = datenum(data_struct.textdata(4:end,1));

%get data matrix
data=data_struct.data(3:end,:);

%find GDP
idx = find(strcmp(series,'GDPC1'));
gdp = 400*diff(log(data(:,idx)));
dates = dates(2:end);

%% Fit student t distribution
idx = find(dates == datenum('Mar-1-1984'));
QQ = 0.05:0.05:0.95; 
deltaYY = 0.1;        
YY = (-20):deltaYY:20;

qq=quantile(gdp,QQ);
[lc, sc, sh, df] = QuantilesInterpolation(qq, QQ);
dist1 = dskt(YY, lc, sc, sh, df);

qq=quantile(gdp(idx:end),QQ);
[lc, sc, sh, df] = QuantilesInterpolation(qq, QQ);
dist2 = dskt(YY, lc, sc, sh, df);

qq=quantile(gdp(1:idx),QQ);
[lc, sc, sh, df] = QuantilesInterpolation(qq, QQ);
dist3 = dskt(YY, lc, sc, sh, df);

loyolagreen = 1/255*[0,104,87];
loyolagray = 1/255*[200,200,200];

% Plot histogram
f=figure;
histogram(gdp, 'BinWidth', 1, 'FaceColor', 'b', 'Normalization', 'pdf', 'FaceAlpha', 0.5);
hold on;
histogram(gdp(idx:end), 'BinWidth', 1, 'FaceColor', 'r', 'Normalization', 'pdf', 'FaceAlpha', 0.5);
hold on;
plot(YY, dist1, 'b--', 'LineWidth', 1.5);
hold on;
plot(YY, dist2, 'r--', 'LineWidth', 1.5);

title('Histogram of Real GDP growth', 'FontSize', 16)
legend('Q2 1959 - Q3 2019','Q1 1984 - Q3 2019','FontSize', 14);
set(gcf,'Position',[800 450 800 450])
printpdf(f, 'output/density_t_distrib');


close all;

