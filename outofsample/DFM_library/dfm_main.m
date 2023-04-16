%% Dynamic Factor Models
%
% The stationary vector process x_t with dimentions n by 1 admits a factor
% model representation:
%
%  (1) x_t = Lambda * f_t + e_t
%
% where Lambda is the n by r matrix of factor loadings and f_t is a r by 1
% vector of factors. Here, lags of f_t ar not included into the model. e_t
% is a serially uncorrelated and follows an AR(1) process.
%
%  (2) e_t = rho * e_t-1 + v_t
%
% The factors follow a stationary VAR process of order p:
%
%  (3) f_t = sum_(p=0...P)(A_p * f_t-p) + u_t
%
% The following parameters will be chosen in the program:
% p: the number of lags in the VAR process (3). p takes on values 1,2.
% r: the number of factors

dfm_initial_settings;
for t=1:StartVintage:EndVintage
    
    if t > StartVintage
        CurrDates = datevec(name_vintage(t));
        PrevDates = datevec(name_vintage(t-1));
        
        if t+1 <= EndVintage
            NextDates = datevec(name_vintage(t+1));
        else
            NextDates = 0;
        end
        
         X_old = X_new;
    end
    
    % Matrix of monthly series up to and inluding release at date t
    X_new=DD{t};
    
    if t > 1
        X_last_iter = DD{t-1};
    elseif t==1
        X_last_iter = X_old;
    end
    
    %% Convert dates from month to quarter
    %
    % P.Qnews: For every month P.Qnews is the last month of the
    %          corresponding quarter
    
    temp = datevec(name_vintage(t));
    temp = temp(1:2);
    
    % Turn months 1,2->3, 4,5->6, 7,8->9, 10,11->12
    if  mod(temp(2),3) == 1
        temp(2)= temp(2)+2;
    elseif mod(temp(2),3) == 2
        temp(2)= temp(2)+1;
    end
    
    P.Qnews = temp;
    
    % Index of DatesV where the date is the same as P.Qnews
    iQ = find(DatesV(:,1) == P.Qnews(1) & DatesV(:,2) == P.Qnews(2));
    
    % 2/3 step forecast
    if addFore == 1
        if iQ+9>size(X_new,1);
            addFore = 0;
            disp('Two and three steps ahead forecast blocked');
        end
    end
    
    
    if strcmp(countryQ, 'CHN')
        if ~isempty(SubS)
            for ut=1:size(SubS,1)
                X_new(SubS(ut,1), SubS(ut,2))=nan;
            end
        end
    end    
    
    
    %%
    
    for iSerInd = 1:size(iSer, 2)
        if QoQ == 1;
            
            if t == 1
                iQ_last_iter = iQ;
            end
            
            % Regress GDP growth on lagged GDP growth to get AR(1) coefficient.
            if  t==1 || (t>2 && CurrDates(:,2)==1 && PrevDates(:,2)==12) || ...
                    (iQ == iQ_last_iter && ~isnan(X_new(iQ-3,iSer(iSerInd)))   && ...
                    isnan(X_last_iter(iQ-3,iSer(iSerInd))))
                
                temp=isfinite(X_new(:,iSer(iSerInd)));                                          % 1 if we there is a GDP release that month, 0 otherwise
                GdP=X_new(temp,iSer(iSerInd));                                                  % vector of quarterly GDP growth data
                intercept = ones(sum(temp)-1,1); % TO BE FIXED FP
                lagged = GdP(1:end-1,:);
                reg=[intercept lagged];
                beta(:, iSerInd) = (reg'*reg)\(reg'*GdP(2:end,:));
            end
            
            % If there is no observation for the previous quarter yet backcast
            % the quarter using estimated the AR(1) process:
            
            if isnan(X_new(iQ-3,iSer(iSerInd)))
                AR_Back(t, iSerInd) = [1 X_new(iQ-6,iSer(iSerInd))]*beta(:, iSerInd);
            else
                AR_Back(t, iSerInd) = X_new(iQ-3,iSer(iSerInd));
            end
        end
        
        AR_Now(t, iSerInd)  = [1 AR_Back(t, iSerInd) ] * beta(:, iSerInd);
        AR_Fore(t, iSerInd) = [1 AR_Now(t, iSerInd) ] * beta(:, iSerInd);
        
        if addFore == 1
            AR_Fore2S(t, iSerInd) = [1 AR_Fore(t, iSerInd) ] * beta(:, iSerInd);
            AR_Fore3S(t, iSerInd) = [1 AR_Fore2S(t, iSerInd) ] * beta(:, iSerInd);
        end
        
%         for iSerInd=1:length(iSer)
%             AR_Now(t, iSerInd)  = [1 AR_Back(t, iSerInd) ] * beta(:, iSerInd);
%             AR_Fore(t, iSerInd) = [1 AR_Now(t, iSerInd) ] * beta(:, iSerInd);
%             
%             if addFore == 1
%                 AR_Fore2S(t, iSerInd) = [1 AR_Fore(t, iSerInd) ] * beta(:, iSerInd);
%                 AR_Fore3S(t, iSerInd) = [1 AR_Fore2S(t, iSerInd) ] * beta(:, iSerInd);
%             end
%         end
        
        
        if YoY == 1;
            
            if t == 1
                iQ_last_iter = iQ;
            end
            
            % Regress GDP growth on lagged GDP growth to get AR(1) coefficient.
            if  t==1 || (t>2 && CurrDates(:,2)==1 && PrevDates(:,2)==12) || ...
                    (iQ == iQ_last_iter && ~isnan(X_new(iQ-3,iSer(iSerInd)))   && ...
                    isnan(X_last_iter(iQ-3,iSer(iSerInd))))
                
                temp=isfinite(X_new(:,iSer(iSerInd)));
                GdP=filter([1 1 1 1],1,X_new(temp,iSer(iSerInd)));
                intercept = ones(size(GdP(5:end-1,:),1),1);
                lagged = GdP(5:end-1,:);
                reg=[intercept lagged];
                beta_yoy(:, iSerInd)=(reg'*reg)\(reg'*GdP(6:end,:));
            end
            
            % If there is no observation for the previous quarter yet backcast
            % the quarter using estimated the AR(1) process:
            
            if isnan(X_new(iQ-3,iSer(iSerInd)))
                AR_Back_yoy(t, iSerInd) = [1 sum(X_new(iQ-15:3:iQ-6,iSer(iSerInd))) ]*beta_yoy(:, iSerInd);
                
            else
                AR_Back_yoy(t, iSerInd) = sum(X_new(iQ-12:3:iQ-3,iSer(iSerInd)));
            end
            
            AR_Now_yoy(t, iSerInd) = [1 AR_Back_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
            AR_Fore_yoy(t, iSerInd) = [1 AR_Now_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
            AR_Fore2S_yoy(t, iSerInd) = [1 AR_Fore_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
            AR_Fore3S_yoy(t, iSerInd) = [1 AR_Fore2S_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);

%             for iSerInd=1:length(iSer)
%                 AR_Now_yoy(t, iSerInd) = [1 AR_Back_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
%                 AR_Fore_yoy(t, iSerInd) = [1 AR_Now_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
%                 AR_Fore2S_yoy(t, iSerInd) = [1 AR_Fore_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
%                 AR_Fore3S_yoy(t, iSerInd) = [1 AR_Fore2S_yoy(t, iSerInd) ]*beta_yoy(:, iSerInd);
%             end
        end
    end
    
    
    %%
    
    if t==1
        dfm_parameters_estimation;
    end
    
    if em_frequency == 1 && (t>2 && CurrDates(:,2)==1 && PrevDates(:,2)==12)
        dfm_parameters_estimation;
    elseif em_frequency == 2 ...
            &&((t>2 && CurrDates(:,2)==1 && PrevDates(:,2)==12) ...
            || (t>2 && CurrDates(:,2)==4 && PrevDates(:,2)==3) ...
            || (t>2 && CurrDates(:,2)==7 && PrevDates(:,2)==6) ...
            || (t>2 && CurrDates(:,2)==10 && PrevDates(:,2)==9))
        dfm_parameters_estimation;
    elseif em_frequency == 3 && ((t>2 && CurrDates(:,2) ~= PrevDates(:,2)))
        dfm_parameters_estimation;
    end
    
    %%
    %if t > StartVintage
    
    dfm_news_estimation;
    % filtBB{1,t} = filtB;
    if length(filtN) > 0
        
        dfm_qoq_yoy;
        dfm_save_output;
        
        % Store the factors - FP
        factorCell{t} = R_new.F;
        
    else
        disp('Exception: variable released after the GDP');
    end
    %end
    
    iQ_last_iter = iQ;
end

if saveAddFore == 1
    addFore = 1;
end