function [y_old,y_new,groupnews,singlenews,gain,gainSer,actual,fore,filt,t_miss,v_miss,innov, factor] = ...
    News_DFM_ML_FP(X_old,X_new,Q,t_fcst,v_news)


    r = size(Q.C,2);
    [T,N] = size(X_new);
    gList = unique(Q.Groups);
    groupnews = zeros(3,length(gList));
    singlenews = zeros(1,N);
    gain = [];
    gainSer = {};
    factor = {};
    
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
    else
        
        % comment TH: Why do we use the mean and std from the EM?
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
        end
        %----------------------------------------------------------------------
        %     v_miss=[1:size(X_new,2)]';
        %     t_miss=t_miss(1)*ones(size(X_new,2),1);
        %----------------------------------------------------------------------
        lag = t_fcst-t_miss; 
        % t_fcst is the last month of the quarter we are forecasting
        % t_miss is the reference month of the release

        k = max([abs(lag); max(lag)-min(lag)]); 
        % max distance between the variables that are released and the final month of ref. quarter

        C = Q.C;
        R = Q.R'; % cov matrix of idio (ALB)

        n_news = size(lag,1);   % can be deleted


        %         Res_old = para_const_FP(X_old, Q, k);
        %         P = Res_old.P(:,:,2:end);
        %         for jt = 1:size(P,3)
        %             for jk = 0:k
        %                 Plag{jk+1}(:,:,jt) = P(1:r,jk*r+1:(jk+1)*r,jt)';
        %             end
        %         end

        Res_old = para_const_FP(X_old, Q, k);
        Plag = Res_old.Plag;


        Res_new = para_const_FP(X_new, Q, 0);

        factor = Res_new.F; % save the factor
        
        y_old = Res_old.X_sm(t_fcst,v_news);
        y_new = Res_new.X_sm(t_fcst,v_news);

        %         if isempty(t_miss)
        %         actual=[];fore=[];filt=[];t_miss=[];v_miss=[];innov=NaN;
        %         return
        %         else
        P = Res_old.P(:,:,2:end);       % it is unused!!! it can be deleted!!
        P1=[];

        for i=1:size(lag,1)
            h = abs(t_fcst-t_miss(i));
            m = max([t_miss(i) t_fcst]);
            if t_miss(i)>t_fcst % "backcast" period (ALB)
                Pp=Plag{h+1}(:,:,m);     %P(1:r,h*r+1:h*r+r,m)';
            else
                Pp=Plag{h+1}(:,:,m)';    %P(1:r,h*r+1:h*r+r,m);
            end
            P1=[P1 Pp*C(v_miss(i),1:r)'];
        end
    
        % Calculate the news
        for i=1:size(t_miss,1)
            X_new_norm = (X_new(t_miss(i),v_miss(i)) - Mx(v_miss(i)))./Wx(v_miss(i));
            X_sm_norm = (Res_old.X_sm(t_miss(i),v_miss(i))- Mx(v_miss(i)))./Wx(v_miss(i));
            innov(i)= X_new_norm-X_sm_norm; % NEWS (ALB)
        end

        ins = size(innov,2);    % unused!! we cna delete it
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

        clear temp
        % loop on v_news
        for i=1:size(v_news,2)
            totnews(1,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1/P2*innov';
            temp(1,:,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1/P2.*innov;
            gain(:,:,i) = Wx(v_news(i))*C(v_news(i),1:r)*P1/P2;
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

        %         for i=1:length(gList)
        %             groupnews(i)=gain(:,ismember(Q.Groups(v_miss),gList(i)))*innov(ismember(Q.Groups(v_miss),gList(i)))';
        %         end

        [v_miss, idx, j] = unique(v_miss);
        gain = gain(:,idx,:);
        gainSer = gainSer(idx);

    end
end

