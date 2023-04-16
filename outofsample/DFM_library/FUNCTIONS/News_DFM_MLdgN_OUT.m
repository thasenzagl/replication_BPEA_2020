function [y_old,y_new,groupnews,singlenews,gain,gainSer,actual,fore,filt,t_miss,v_miss,innov,SubS, factor] = ...
    News_DFM_MLdgN_OUT(X_old,X_new,Q,t_fcst,v_news,nM,X_Last)


r = size(Q.C,2);
[T,N] = size(X_new);
gList = unique(Q.Groups);
groupnews = zeros(1,length(gList));
singlenews = zeros(1,N);
gain = [];
gainSer = {};

if ~isnan(X_new(t_fcst,v_news)) % ALL should be not NAN to be TRUE
    Res_old = para_const_FP(X_old, Q, 0);

    for i=1:size(v_news,2)
        temp = X_new(t_fcst,v_news(i)) - Res_old.X_sm(t_fcst,v_news(i));
        singlenews(:,v_news(i)) = temp;
        groupnews(:,ismember(gList,Q.Groups(v_news(i)))) = temp;
        y_old(1,i) = Res_old.X_sm(t_fcst,v_news(i));
        y_new(1,i) = X_new(t_fcst,v_news(i));
    end
    groupnews=[];singlenews=[];gain=[];gainSer=[];
    actual=[];fore=[];filt=[];t_miss=[];v_miss=[];innov=[];
    SubS=[];factor=[];is_fcst=[];is_Qam=[];
else

    Mx = Q.Mx;
    Wx = Q.Wx;

    miss_old=isnan(X_old);
    miss_new=isnan(X_new);
    temp=miss_old-miss_new;
    [t_miss,v_miss]=find(temp==1); %find release (time & variable)
    if isempty(v_miss)
        y_old=[]; y_new=[];groupnews=[];singlenews=[];gain=[];gainSer=[];
        actual=[];fore=[];filt=[];t_miss=[];v_miss=[];innov=[];
        SubS=[];factor=[];is_fcst=[];is_Qam=[];
        return
    end
    %----------------------------------------------------------------------
    lag = t_fcst-t_miss; %distance of release from forecast horizon

    k = max([abs(lag); max(lag)-min(lag)]);

    C = Q.C;
    R = Q.R';

    n_news = size(lag,1);


    Res_old = para_const_FP(X_old, Q, k);
    Plag = Res_old.Plag; %covariance of states up to lag k


    Res_new = para_const_FP(X_new, Q, 0);
    
    factor = Res_new.F; % save the factor

    y_old = Res_old.X_sm(t_fcst,v_news);
    y_new = Res_new.X_sm(t_fcst,v_news);
    
    %P1 is E[(X_{t}-X_{t|t-k})(S_{t}-S_{t|t-k})']
    P = Res_old.P(:,:,2:end); %variance of the states
    P1=[];

    for i=1:size(lag,1)
        h = abs(t_fcst-t_miss(i)); %fcst horizon
        m = max([t_miss(i) t_fcst]); 
        if t_miss(i)>t_fcst
            Pp=Plag{h+1}(:,:,m);
        else
            Pp=Plag{h+1}(:,:,m)'; %covariance of states t&(t-k) @time max(t_miss,t_fcst)
        end
        P1=[P1 Pp*C(v_miss(i),1:r)']; %cov states with X (released variable)
    end

    for i=1:size(t_miss,1)
        X_new_norm = (X_new(t_miss(i),v_miss(i)) - Mx(v_miss(i)))./Wx(v_miss(i));
        X_sm_norm = (Res_old.X_sm(t_miss(i),v_miss(i))- Mx(v_miss(i)))./Wx(v_miss(i));
        innov(i)= X_new_norm-X_sm_norm;

    end
    
    %P2 is E[(X_{t}-X_{t|t-k})(X_{t}-X_{t|t-k})']
    ins=size(innov,2);
    P2=[];
    p2=[];
    for i=1:size(lag,1)
        for j=1:size(lag,1)
            h=abs(lag(i)-lag(j));
            m=max([t_miss(i),t_miss(j)]);

            if t_miss(j)>t_miss(i)
                Pp=Plag{h+1}(:,:,m); %P(1:r,h*r+1:(h+1)*r,m)';
            else
                Pp=Plag{h+1}(:,:,m)'; %P(1:r,h*r+1:(h+1)*r,m);
            end
            if v_miss(i)==v_miss(j) & t_miss(i)~=t_miss(j) %multiple release for same variable
                WW(v_miss(i),v_miss(j))=0;
            else
                WW(v_miss(i),v_miss(j))=R(v_miss(i),v_miss(j)); %idio variance
            end
            p2=[p2 C(v_miss(i),1:r)*Pp*C(v_miss(j),1:r)'+WW(v_miss(i),v_miss(j))]; %cov of X up to max distance lag of two releases
        end
        P2=[P2;p2];
        p2=[];
    end

    clear temp
    % loop on v_news

    SubS=[];
    for jr=1:size(innov,1)
        if P2(jr,jr)*(4^2)<(innov(:,jr))^2
            SubS=[SubS; t_miss(jr) v_miss(jr)] ;
            [P2(jr,jr) (innov(:,jr))^2 innov(:,jr)]
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(SubS)
        for ut=1:size(SubS,1)
            X_new(SubS(ut,1),SubS(ut,2))=Res_old.X_sm(SubS(ut,1),SubS(ut,2));
        end
        [y_old,y_new,groupnews,singlenews,gain,gainSer,actual,fore,filt,t_miss,v_miss,innov] = ...
            News_DFM_MLdg(X_old,X_new,Q,t_fcst,v_news);
        %[y_old,y_new,groupnews,singlenews,gain,gainSer,actual,fore,filt,t_miss,v_miss,innov] = ...
        %    News_DFM_MLdgN(X_old,X_new,Q,t_fcst,v_news);
    else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:size(v_news,2)
            totnews(1,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1*inv(P2)*innov';
            temp(1,:,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1*inv(P2).*innov;
            gain(:,:,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1*inv(P2);
        end
        gainSer = Q.Series(v_miss);

        singlenews = zeros(max(t_miss)-min(t_miss)+1,N); %,size(v_news,2)
        actual = zeros(N,1);
        fore = zeros(N,1);
        filt = zeros(N,1,size(v_news,2));

        for i=1:size(innov,2)

            actual(v_miss(i),1) = X_new(t_miss(i),v_miss(i));
            fore(v_miss(i),1) = Res_old.X_sm(t_miss(i),v_miss(i));

            for j=1:size(v_news,2)
                singlenews(t_miss(i)-min(t_miss)+1,v_miss(i),j) = temp(1,i,j);
                filt(v_miss(i),:,j) = gain(:,i,j)/Wx(v_miss(i));
            end
        end

        singlenews = sum(singlenews,1);

        [v_miss, idx, j] = unique(v_miss);
        gain = gain(:,idx,:);
        gainSer = gainSer(idx);

    end
end