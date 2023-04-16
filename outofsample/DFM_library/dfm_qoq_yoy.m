% quarter on quarter
if QoQ == 1
    Nowcast(t,:)  = y_newN;
    Forecast(t,:) = y_newF;
    Backcast(t,:) = y_newB;
    ReleaseTime(t,:) = DDname{t};
    
    if addFore == 1
        Forecast2S(t,:) = y_newF2S;
        Forecast3S(t,:) = y_newF3S;
    end
    
    for iSerInd=1:length(iSer)
        if isnan(X_new(iQ-3,iSer(iSerInd)))
            AR_Back(t, iSerInd) = [1 X_new(iQ-6,iSer(iSerInd)) ]*beta(:, iSerInd);
        else
            AR_Back(t, iSerInd) = X_new(iQ-3,iSer(iSerInd));
        end
        
        AR_Now(t, iSerInd) = [1 AR_Back(t, iSerInd)]*beta(:, iSerInd);
        AR_Fore(t, iSerInd) = [1 AR_Now(t, iSerInd) ]*beta(:, iSerInd);
        
        if addFore ==1
            AR_Fore2S(t, iSerInd) = [1 AR_Fore(t, iSerInd) ] * beta(:, iSerInd);
            AR_Fore3S(t, iSerInd) = [1 AR_Fore2S(t, iSerInd) ] * beta(:, iSerInd);
        end
        
        CurrentOutputN{iSerInd, t} = [];
        CurrentOutputB{iSerInd, t} = [];
        CurrentOutputF{iSerInd, t} = [];
        
        if addFore == 1
            CurrentOutputF2S{iSerInd, t} = [];
            CurrentOutputF3S{iSerInd, t} = [];
        end
        
    end
    clear iSerInd
    
end


%% year on year
if YoY == 1
    XX     = X_new(1:iQ+3, iSer);
    XX_old = X_old(1:iQ+3, iSer);
    
    for iSerInd=1:length(iSer)
        if isnan(X_new(iQ-3, iSer(iSerInd)))
            XX(iQ-3, iSerInd) = y_newB(:, iSerInd);
        end
        
        if isnan(X_old(iQ-3, iSer(iSerInd)))
            XX_old(iQ-3, iSerInd) = y_oldB(:, iSerInd);
        end
    end
    
    XX(iQ, :)     = y_newN;
    XX_old(iQ, :) = y_oldN;
    
    XX(iQ+3, :)     = y_newF;
    XX_old(iQ+3, :) = y_oldF;
    
    XX(iQ+6, :)     = y_newF2S;
    XX_old(iQ+6, :) = y_oldF2S;
    
    XX(iQ+9, :)     = y_newF3S;
    XX_old(iQ+9, :) = y_oldF3S;
    
    Nowcast_yoy(t, :)  = sum(XX(iQ-9:3:iQ, :));
    Forecast_yoy(t, :) = sum(XX(iQ-6:3:iQ+3, :));
    Backcast_yoy(t, :) = sum(XX(iQ-12:3:iQ-3, :));
    Forecast2S_yoy(t, :) = sum(XX(iQ-3:3:iQ+6, :));
    Forecast3S_yoy(t, :) = sum(XX(iQ:3:iQ+9, :));
    
    ReleaseTime(t, :)  = DDname{t};
    
    for iSerInd=1:length(iSer)
        if isnan(X_new(iQ-3, iSer(iSerInd)))
            AR_Back_yoy(t, iSerInd) = [1 sum(X_new(iQ-15:3:iQ-6, iSer(iSerInd))) ]*beta_yoy(:, iSerInd);
        else
            AR_Back_yoy(t, iSerInd) = sum(X_new(iQ-12:3:iQ-3, iSer(iSerInd)));
        end
        AR_Now_yoy(t, iSerInd)  = [1 AR_Back_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
        AR_Fore_yoy(t, iSerInd) = [1 AR_Now_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
        AR_Fore2S_yoy(t, iSerInd) = [1 AR_Fore_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
        AR_Fore3S_yoy(t, iSerInd) = [1 AR_Fore2S_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
        
        CurrentOutputN_yoy{iSerInd, t} = [];
        CurrentOutputB_yoy{iSerInd, t} = [];
        CurrentOutputF_yoy{iSerInd, t} = [];
        CurrentOutputF2S_yoy{iSerInd, t} = [];
        CurrentOutputF3S_yoy{iSerInd, t} = [];
        
        if addFore == 1
            CurrentOutputF2S{iSerInd, t} = [];
            CurrentOutputF3S{iSerInd, t} = [];
        end
        
    end
    clear iSerInd
    
end