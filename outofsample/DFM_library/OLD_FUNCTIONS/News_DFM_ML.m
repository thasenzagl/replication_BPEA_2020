function [y_old,y_new,groupnews,singlenews,gain,gainSer,actual,fore,filt,t_miss,v_miss,innov] = ...
    News_DFM_ML(X_old,X_new,Q,t_fcst,v_news)


r = size(Q.C,2);
[T,N] = size(X_new);
gList = unique(Q.Groups);
groupnews = zeros(1,length(gList));
%singlenews = NaN(1,N); %zeros
singlenews = zeros(1,N); %zeros
gain = [];
gainSer = {};

if ~isnan(X_new(t_fcst,v_news))
    Res_old = para_const(X_old, Q, 0);
    temp = X_new(t_fcst,v_news) - Res_old.X_sm(t_fcst,v_news);
    %singlenews(:,v_news) = NaN;
    singlenews(:,v_news) = 0;
    groupnews(:,ismember(gList,Q.Groups(v_news))) = temp;
    y_old = Res_old.X_sm(t_fcst,v_news);
    %y_new = num2cell(X_new(t_fcst,v_news));
    %y_new = single(X_new(t_fcst,v_news));
    %y_new =single(NaN);
    y_new =single(0);

       
    %%to check and decide what to send!!!
%   %-----------------
% 
    miss_old=isnan(X_old);
    miss_new=isnan(X_new);
    temp=miss_old-miss_new;
    [t_miss,v_miss]=find(temp==1);
    
    %gain=cell(1,size(v_miss,1));
    gain=NaN(1,size(v_miss,1));
    %-----------------
    gainSer=[];%gain=[];
    fore=[];filt=[];t_miss=[];v_miss=[];
    innov=[];
    actual=[];
    return 

else
    
    Mx = Q.Mx;
    Wx = Q.Wx;
    
    miss_old=isnan(X_old);
    miss_new=isnan(X_new);
    temp=miss_old-miss_new;
    [t_miss,v_miss]=find(temp==1);
    if isempty(v_miss)
    y_old=[]; y_new=[];groupnews=[];singlenews=[];gain=[];gainSer=[];
    actual=[];fore=[];filt=[];t_miss=[];v_miss=[];innov=[];
    return 
    
    end    %----------------------------------------------------------------------
%     v_miss=[1:size(X_new,2)]';
%     t_miss=t_miss(1)*ones(size(X_new,2),1);
    %----------------------------------------------------------------------
    lag = t_fcst-t_miss;
        
    k = max([abs(lag); max(lag)-min(lag)]);
    
    C = Q.C;
    R = Q.R';
    
    n_news = size(lag,1);
    
   
        
        Res_old = para_const(X_old, Q, k);
        Res_new = para_const(X_new, Q, 0);
        
        y_old = Res_old.X_sm(t_fcst,v_news);
        y_new = Res_new.X_sm(t_fcst,v_news);
        
        if isempty(t_miss)
        actual=[];fore=[];filt=[];t_miss=[];v_miss=[];innov=NaN;
        return
        else
        P = Res_old.P(:,:,2:end);
        P1=[];
        
        for i=1:size(lag,1)
            h = abs(t_fcst-t_miss(i));
            m = max([t_miss(i) t_fcst]);
            if t_miss(i)>t_fcst
                Pp=P(1:r,h*r+1:h*r+r,m)';
            else
                Pp=P(1:r,h*r+1:h*r+r,m);
            end
            P1=[P1 Pp*C(v_miss(i),1:r)'];
        end
        
        for i=1:size(t_miss,1)
            X_new_norm = (X_new(t_miss(i),v_miss(i)) - Mx(v_miss(i)))./Wx(v_miss(i));
            X_sm_norm = (Res_old.X_sm(t_miss(i),v_miss(i))- Mx(v_miss(i)))./Wx(v_miss(i));
            innov(i)= X_new_norm-X_sm_norm;
        end
        
        ins=size(innov,2);
        P2=[];
        p2=[];
        for i=1:size(lag,1)
            for j=1:size(lag,1)
                h=abs(lag(i)-lag(j));
                m=max([t_miss(i),t_miss(j)]);
                if t_miss(j)>t_miss(i)
                    Pp=P(1:r,h*r+1:(h+1)*r,m)';
                else
                    Pp=P(1:r,h*r+1:(h+1)*r,m);
                end
                if v_miss(i)==v_miss(j) & t_miss(i)~=t_miss(j)
                    WW(v_miss(i),v_miss(j))=0;
                else
                    WW(v_miss(i),v_miss(j))=R(v_miss(i),v_miss(j));
                end
                p2=[p2 C(v_miss(i),1:r)*Pp*C(v_miss(j),1:r)'+WW(v_miss(i),v_miss(j))];
            end
            P2=[P2;p2];
            p2=[];
        end
        
        
        
        totnews = Wx(v_news)*C(v_news,1:r)*P1*inv(P2)*innov';
        temp = Wx(v_news)*C(v_news,1:r)*P1*inv(P2).*innov;
        gain = Wx(v_news)*C(v_news,1:r)*P1*inv(P2);
        gainSer = Q.Series(v_miss);
        
        singlenews = zeros(max(t_miss)-min(t_miss)+1,N);
        actual = zeros(N,1);
        fore = zeros(N,1);
        filt = zeros(N,1);
        
        for i=1:size(innov,2)
            singlenews(t_miss(i)-min(t_miss)+1,v_miss(i)) = temp(i);
            actual(v_miss(i),1) = X_new(t_miss(i),v_miss(i));
            fore(v_miss(i),1) = Res_old.X_sm(t_miss(i),v_miss(i));
            filt(v_miss(i),:) = gain(i)/Wx(v_miss(i));
        end
        
        singlenews = sum(singlenews,1);
        
        %         for i=1:length(gList)
        %             groupnews(i)=gain(:,ismember(Q.Groups(v_miss),gList(i)))*innov(ismember(Q.Groups(v_miss),gList(i)))';
        %         end
        
        [v_miss, idx, j] = unique(v_miss);
        gain = gain(idx);
        gainSer = gainSer(idx);
        
    end
end

