if EM_case == 0 % default
    R_new = EM_DFM_SS_block_idioQARMA_restrMQ_FP(X_new(1:iQ+9,:), P, []);
    
elseif EM_case == 1 % canada
    R_new = EM_DFM_SS_block_idioQARMA_restrMQ_VAR(X_new(1:iQ+9,:), P);
    
else % china
    R_new = EM_DFM_SS_idio_restrMQ_stat_out(X_new(1:iQ+9,:), P, Res_old);
    Res_old = R_new;
end

R_new.Groups = Group;
R_new.Series = SeriesID;

%% quarter on quarter
% Nowcast, Forecast and Backcast on release date t. iSer is index
% of variable we are estimating.

if QoQ == 1;
    
    ReleaseTime(t,:) = DDname{t} + datenum('12-30-1899');
    Nowcast(t,:)     = R_new.X_sm(iQ,iSer); % If error: check the name of the series
    Forecast(t,:)    = R_new.X_sm(iQ+3,iSer);
    Backcast(t,:)    = R_new.X_sm(iQ-3,iSer);
    
    if addFore == 1
        Forecast2S(t,:)=R_new.X_sm(iQ+6,iSer);
        Forecast3S(t,:)=R_new.X_sm(iQ+9,iSer);
    end
    
    % Nowcast the current quarter and forecast the next quarter using
    % the AR(1) process
    
    clear iSerInd
end


%% year on year
if YoY == 1
    
    XX = X_new(1:iQ+3, iSer);
    
    for iSerInd=1:length(iSer)
        if isnan(X_new(iQ-3, iSer(iSerInd)))
            XX(iQ-3, iSerInd) = R_new.X_sm(iQ-3,iSer(iSerInd));
        end
    end
    
    XX(iQ, :)   = R_new.X_sm(iQ, iSer);
    XX(iQ+3, :) = R_new.X_sm(iQ+3, iSer);
    
    XX(iQ+6, :) = R_new.X_sm(iQ+6, iSer);
    XX(iQ+9, :) = R_new.X_sm(iQ+9, iSer);
    
    Nowcast_yoy(t,:)  = sum(XX(iQ-9:3:iQ,1));
    Forecast_yoy(t,:) = sum(XX(iQ-6:3:iQ+3,1));
    Backcast_yoy(t,:) = sum(XX(iQ-12:3:iQ-3,1));
    Forecast2S_yoy(t,:) = sum(XX(iQ-3:3:iQ+6,1));
    Forecast3S_yoy(t,:) = sum(XX(iQ:3:iQ+9,1));
    
    ReleaseTime(t,:) = DDname{t};
    
end