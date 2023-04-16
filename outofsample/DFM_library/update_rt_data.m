function [P, iSer, iSer_idx]  = update_rt_data(P, DD_series, blocks, include)

iSer_idx = ~ismember(DD_series,P.SerNews);
iSer = find(ismember(DD_series,P.SerNews))';
series_idx=ismember(P.all_series,DD_series(iSer_idx));

% Blocks
P.include = logical(include(series_idx,:));
P.blocks = blocks(series_idx,:);
P.blocks = P.blocks(P.include,:);

% Series
series_idx_gdp=ismember(P.all_series,DD_series);
include_gdp = logical(include(series_idx_gdp,:));
P.Series = DD_series(include_gdp);

% Other parameters 
P.nM=length(P.blocks);
P.nQ = 0;
P.i_idio = logical([ones(P.nM,1);zeros(P.nQ,1)]);

end
