if EM_case == 2 % CHINA
    [y_oldB,y_newB,groupnews,singlenewsB,gainB,gainSer,actual,fore,filtB,t_miss,v_miss,innov]=...
        News_DFM_MLdgN_OUT(X_old(1:iQ+3,:),X_new(1:iQ+3,:),R_new,iQ-3,iSer);
    
    [y_oldN,y_newN,groupnewsN,singlenewsN,gainN,gainSerN,actual,fore,filtN,t_miss,v_miss,innov,temp_P, factor]=...
        News_DFM_MLdgN_OUT(X_old(1:iQ+3,:),X_new(1:iQ+3,:),R_new,iQ,iSer);
    
    [y_oldF,y_newF,groupnews,singlenewsF,gainF,gainSer,actual,fore,filtF,t_miss,v_miss,innov]=...
        News_DFM_MLdgN_OUT(X_old(1:iQ+3,:),X_new(1:iQ+3,:),R_new,iQ+3,iSer);
    
    SubS=[SubS;temp_P];
    RealTimeFactor{t}=factor;
    
else % NOT CHINA
    [y_oldB,y_newB,groupnews,singlenewsB,gainB,gainSer,actual,fore,filtB,t_miss,v_miss,innov]=...
        News_DFM_ML_FP(X_old(1:iQ+3,:),X_new(1:iQ+3,:),R_new,iQ-3,iSer);
    
    [y_oldN,y_newN,groupnewsN,singlenewsN,gainN,gainSerN,actual,fore,filtN,t_miss,v_miss,innov,factor]=...
        News_DFM_ML_FP(X_old(1:iQ+3,:),X_new(1:iQ+3,:),R_new,iQ,iSer);
    
    [y_oldF,y_newF,groupnews,singlenewsF,gainF,gainSer,actual,fore,filtF,t_miss,v_miss,innov]=...
        News_DFM_ML_FP(X_old(1:iQ+3,:),X_new(1:iQ+3,:),R_new,iQ+3,iSer);
    
    RealTimeFactor{t}=factor; % save the factor
end

if length(filtN) > 0
    
    if length(filtB) == 0
        filtB = zeros(size(X_new, 2), 1);
    end
    
    if addFore == 1
        if EM_case == 2
            [y_oldF2S,y_newF2S,groupnews,singlenewsF2S,gainF2S,gainSer,actual,fore,filtF2S,t_miss,v_miss,innov]=...
                News_DFM_MLdgN_OUT(X_old(1:iQ+9,:),X_new(1:iQ+9,:),R_new,iQ+6,iSer);
            
            [y_oldF3S,y_newF3S,groupnewsN,singlenewsF3S,gainF3S,gainSerN,actual,fore,filtF3S,t_miss,v_miss,innov,~,factor3S]=...
                News_DFM_MLdgN_OUT(X_old(1:iQ+9,:),X_new(1:iQ+9,:),R_new,iQ+9,iSer);
            
            RealTimeFactor3S{t}=factor3S;
            
        else
            [y_oldF2S,y_newF2S,groupnews,singlenewsF2S,gainF2S,gainSer,actual,fore,filtF2S,t_miss,v_miss,innov]=...
                News_DFM_ML_FP(X_old(1:iQ+9,:),X_new(1:iQ+9,:),R_new,iQ+6,iSer);
            
            [y_oldF3S,y_newF3S,groupnewsN,singlenewsF3S,gainF3S,gainSerN,actual,fore,filtF3S,t_miss,v_miss,innov,factor3S]=...
                News_DFM_ML_FP(X_old(1:iQ+9,:),X_new(1:iQ+9,:),R_new,iQ+9,iSer);
            
            RealTimeFactor3S{t}=factor3S;
        end
    end
end