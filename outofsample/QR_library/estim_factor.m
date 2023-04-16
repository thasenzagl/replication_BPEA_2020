function [frQ, ffQ, y, datesQ, jt, first_not_na] = estim_factor(X, P, DatesV, iQ, target, annualize)
    
    T = size(DatesV,1);
    
    % Run EM algo
    R_new = EM_DFM_SS_block_idioQARMA_restrMQ_FP(X, P, []);
    frM = R_new.F(:,1);
    ffM = R_new.F(:,6);

    % Convert to Quarterly
    start = find(DatesV(:,2) == 1 | DatesV(:,2) == 4 | DatesV(:,2) == 7 | DatesV(:,2) == 10, 1, 'first');
    
    datesQ = [];
    frQ = [];
    ffQ = [];
    y = [];
    for i = start:3:floor(T/3)*3
        datesQ = [datesQ; DatesV(i+2,:)];
        frQ = [frQ; mean(frM(i:i+2))]; 
        ffQ = [ffQ; mean(ffM(i:i+2))];         
        y = [y; target(i+2)];
    end    
    
    jt = find(datesQ(:,1) == DatesV(iQ,1) & datesQ(:,2) == DatesV(iQ,2));
   
    first_not_na = find(~isnan(y), 1, 'first');
    last_not_na = find(~isnan(y), 1, 'last');

    % Flip factor if necessary
    sr = sign(corr(frQ(first_not_na:last_not_na), y(first_not_na:last_not_na,end)));
    frQ = sr*frQ;
    sf = (-1)*sign(corr(ffQ(first_not_na:last_not_na), y(first_not_na:last_not_na,end)));  
    ffQ = sf*ffQ;
    
    if annualize == 1
        y = 400 * y;
    end    
end