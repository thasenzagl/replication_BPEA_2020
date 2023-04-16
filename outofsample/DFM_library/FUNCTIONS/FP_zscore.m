function [ y, My, Wy ] = FP_zscore( y )
%FP_ZSCORE Zscore
%Author: Filippo Pellegrino

    My = nanmean(y, 1);
    Wy = nanstd(y, 0, 1);
    
    y  = (y-repmat(My, size(y, 1), 1))./repmat(Wy, size(y, 1), 1);
end