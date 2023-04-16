for iSerInd=1:length(iSer)
    for jv = 1:length(v_miss)
        
        %% quarter on quarter
        
        if QoQ == 1
            temp = [ReleaseTime(t) y_oldN(:, iSerInd) y_newN(:, iSerInd) ...
                actual(v_miss(jv)) filtN(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                filtN(v_miss(jv),iSerInd)*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ,:)) X_Last(iQ,iSer(iSerInd)) AR_Now(t, iSerInd)];
            
            CurrentOutputN{iSerInd, t} = [CurrentOutputN{iSerInd, t}; temp];
            
            temp = [ReleaseTime(t) y_oldB(:, iSerInd) y_newB(:, iSerInd) ...
                actual(v_miss(jv)) filtB(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                filtB(v_miss(jv),iSerInd)*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ-3,:)) X_Last(iQ-3,iSer(iSerInd)) AR_Back(t, iSerInd)];
            
            CurrentOutputB{iSerInd, t} = [CurrentOutputB{iSerInd, t}; temp];
            
            temp = [ReleaseTime(t) y_oldF(:, iSerInd) y_newF(:, iSerInd) ...
                actual(v_miss(jv)) filtF(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                filtF(v_miss(jv),iSerInd)*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ+3,:)) X_Last(iQ+3,iSer(iSerInd)) AR_Fore(t, iSerInd)];
            
            CurrentOutputF{iSerInd, t} = [CurrentOutputF{iSerInd, t}; temp];
            
            if addFore == 1
                
                temp = [ReleaseTime(t) y_oldF2S(:, iSerInd) y_newF2S(:, iSerInd) ...
                    actual(v_miss(jv)) filtF2S(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                    filtF2S(v_miss(jv),iSerInd)*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ+6,:)) X_Last(iQ+6,iSer(iSerInd)) AR_Fore2S(t, iSerInd)];
                
                CurrentOutputF2S{iSerInd, t} = [CurrentOutputF2S{iSerInd, t}; temp];
                
                temp = [ReleaseTime(t) y_oldF3S(:, iSerInd) y_newF3S(:, iSerInd) ...
                    actual(v_miss(jv)) filtF3S(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                    filtF3S(v_miss(jv),iSerInd)*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ+9,:)) X_Last(iQ+9,iSer(iSerInd)) AR_Fore3S(t, iSerInd)];
                
                CurrentOutputF3S{iSerInd, t} = [CurrentOutputF3S{iSerInd, t}; temp];
                
                tMt=datevec(ReleaseTime(t));
            end
        end
        %% year on year
        if YoY == 1
            temp = [ReleaseTime(t) sum(XX_old(iQ-9:3:iQ, iSerInd)) sum(XX(iQ-9:3:iQ, iSerInd))...
                actual(v_miss(jv)) filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                (filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd))*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ,:))...
                sum(X_Last(iQ-9:3:iQ,iSer(iSerInd))) AR_Now_yoy(t, iSerInd)];
            
            CurrentOutputN_yoy{iSerInd, t} = [CurrentOutputN_yoy{iSerInd, t}; temp];
            
            temp = [ReleaseTime(t) sum(XX_old(iQ-12:3:iQ-3, iSerInd)) sum(XX(iQ-12:3:iQ-3, iSerInd))...
                actual(v_miss(jv)) filtB(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                filtB(v_miss(jv),iSerInd)*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ-3,:))...
                sum(X_Last(iQ-12:3:iQ-3,iSer(iSerInd))) AR_Back_yoy(t, iSerInd)];
            
            CurrentOutputB_yoy{iSerInd, t} = [CurrentOutputB_yoy{iSerInd, t}; temp];
            
            temp = [ReleaseTime(t) sum(XX_old(iQ-6:3:iQ+3, iSerInd))  sum(XX(iQ-6:3:iQ+3, iSerInd))...
                actual(v_miss(jv)) filtF(v_miss(jv),iSerInd)+filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                (filtF(v_miss(jv),iSerInd)+filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd))*(actual(v_miss(jv))-fore(v_miss(jv)))...
                v_miss(jv) datenum(Dates(iQ+3,:)) sum(X_Last(iQ-6:3:iQ+3,iSer(iSerInd))) AR_Fore_yoy(t, iSerInd)];
            
            CurrentOutputF_yoy{iSerInd, t} = [CurrentOutputF_yoy{iSerInd, t}; temp];
            
            temp = [ReleaseTime(t) sum(XX_old(iQ-3:3:iQ+6, iSerInd)) sum(XX(iQ-3:3:iQ+6, iSerInd))...
                actual(v_miss(jv)) filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                (filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd))*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ+6,:))...
                sum(X_Last(iQ-3:3:iQ+6,iSer(iSerInd))) AR_Fore2S_yoy(t, iSerInd)];
            
            CurrentOutputF2S_yoy{iSerInd, t} = [CurrentOutputF2S_yoy{iSerInd, t}; temp];
            
            temp = [ReleaseTime(t) sum(XX_old(iQ:3:iQ+9, iSerInd)) sum(XX(iQ:3:iQ+9, iSerInd))...
                actual(v_miss(jv)) filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd) actual(v_miss(jv))-fore(v_miss(jv))...
                (filtN(v_miss(jv),iSerInd)+filtB(v_miss(jv),iSerInd))*(actual(v_miss(jv))-fore(v_miss(jv))) v_miss(jv) datenum(Dates(iQ+9,:))...
                sum(X_Last(iQ:3:iQ+9,iSer(iSerInd))) AR_Fore3S_yoy(t, iSerInd)];
            
            CurrentOutputF3S_yoy{iSerInd, t} = [CurrentOutputF3S_yoy{iSerInd, t}; temp];
        end
    end
end

disp(' ')
disp('--------------------------')

disp(['Now running ', datestr(ReleaseTime(t),'dd-mm-yyyy HH:MM')])
disp(' ')

%%
for iSerInd=1:length(iSer)
    structName = strcat('var', num2str(iSerInd));
    if QoQ == 1
        disp(SeriesID{iSer(iSerInd)});
        disp('Released     Gain         News      Impact   Variable')
        for jv = 1:length(v_miss)
            disp([num2str(CurrentOutputN{iSerInd, t}(jv,4:7)),'     ', SeriesID{CurrentOutputN{1, t}(jv,8)},' ',Descr{CurrentOutputN{1, t}(jv,8)}])
        end;
        disp(' ')
        disp('Quarter on Quarter:')
        
        disp(['Previous nowcast: ', num2str(CurrentOutputN{iSerInd, t}(end,2))])
        disp(['New      nowcast: ', num2str(CurrentOutputN{iSerInd, t}(end,3))])
        disp(['Actual   GDP:     ', num2str(CurrentOutputN{iSerInd, t}(end, 10))])
        disp(' ')
        
        if isfield(OutputPrintN, structName)
            OutputPrintN.(structName)   = [OutputPrintN.(structName); CurrentOutputN{iSerInd, t}];
            OutputPrintB.(structName)   = [OutputPrintB.(structName); CurrentOutputB{iSerInd, t}];
            OutputPrintF.(structName)   = [OutputPrintF.(structName); CurrentOutputF{iSerInd, t}];
            
            if addFore == 1
                OutputPrintF2S.(structName)   = [OutputPrintF2S.(structName); CurrentOutputF2S{iSerInd, t}];
                OutputPrintF3S.(structName)   = [OutputPrintF3S.(structName); CurrentOutputF3S{iSerInd, t}];
            end
            
        else
            OutputPrintN.(structName) = CurrentOutputN{iSerInd, t};
            OutputPrintB.(structName) = CurrentOutputB{iSerInd, t};
            OutputPrintF.(structName) = CurrentOutputF{iSerInd, t};
            
            if addFore == 1
                OutputPrintF2S.(structName) = CurrentOutputF2S{iSerInd, t};
                OutputPrintF3S.(structName) = CurrentOutputF3S{iSerInd, t};
            end
        end
    end
    
    if YoY == 1
        if QoQ ~= 1
            disp('Released     Gain         News      Impact   Variable')
            for jv = 1:length(v_miss)
                disp([num2str(CurrentOutputN_yoy{iSerInd, t}(jv,4:7)),'     ',SeriesID{CurrentOutputN_yoy{iSerInd, t}(jv,8)},' ',Descr{CurrentOutputN_yoy{iSerInd, t}(jv,8)}])
            end;
        end
        
        disp(' ')
        disp('Year on Year:')
        
        disp(['Previous nowcast: ', num2str(CurrentOutputN_yoy{iSerInd, t}(end,2))])
        disp(['New      nowcast: ', num2str(CurrentOutputN_yoy{iSerInd, t}(end,3))])
        disp(['Actual   GDP:     ', num2str(CurrentOutputN_yoy{iSerInd, t}(end, 10))])
        disp(' ')
        
        if isfield(OutputPrintN_yoy, structName)
            OutputPrintN_yoy.(structName) = [OutputPrintN_yoy.(structName); CurrentOutputN_yoy{iSerInd, t}];
            OutputPrintB_yoy.(structName) = [OutputPrintB_yoy.(structName); CurrentOutputB_yoy{iSerInd, t}];
            OutputPrintF_yoy.(structName) = [OutputPrintF_yoy.(structName); CurrentOutputF_yoy{iSerInd, t}];
            OutputPrintF2S_yoy.(structName) = [OutputPrintF2S_yoy.(structName); CurrentOutputF2S_yoy{iSerInd, t}];
            OutputPrintF3S_yoy.(structName) = [OutputPrintF3S_yoy.(structName); CurrentOutputF3S_yoy{iSerInd, t}];
        else
            OutputPrintN_yoy.(structName) = CurrentOutputN_yoy{iSerInd, t};
            OutputPrintB_yoy.(structName) = CurrentOutputB_yoy{iSerInd, t};
            OutputPrintF_yoy.(structName) = CurrentOutputF_yoy{iSerInd, t};
            OutputPrintF2S_yoy.(structName) = CurrentOutputF2S_yoy{iSerInd, t};
            OutputPrintF3S_yoy.(structName) = CurrentOutputF3S_yoy{iSerInd, t};
        end
    end
end