function funNewsQK(P,DD,DatesV)

Qnews = P.Qnews;
SerNews = P.SerNews;
Series = P.Series;
Group = P.Group;
P.max_iter = 500;
nM = P.nM;
nQ = P.nQ;
%--------------------------------------------------------------------------
% restrictions
%--------------------------------------------------------------------------
P.Rconstr = [2 -1 0 0 0;...
        3 0 -1 0 0;...
        2 0 0 -1 0;...
        1 0 0 0 -1];
P.q = zeros(4,1);
P.restr = '_restrMQ';
%--------------------------------------------------------------------------
% out-of-sample evaluation
%--------------------------------------------------------------------------
iQ = find(DatesV(:,1) == Qnews(1) & DatesV(:,2) == Qnews(2)); 
iSer = find(ismember(Series,SerNews));  
P.i_idio = logical([ones(nM,1);zeros(nQ,1)]);

%---------------------

for i= 9: size(DD,2)
    X_old=DD{i-1};
    X_new=DD{i};
    if i == 9
        eval(['R_new = EM_DFM_SS',P.method,P.idio,P.restr,'(X_new,P,Res);'])
        R_new.Groups = Group;
        R_new.Series = Series;
    end
    [y_old,y_new,groupnews,singlenews,gain,gainSer,actual,fore,filt]=...
        News_DFM_ML(X_old(1:iQ,:),X_new(1:iQ,:),R_new,iQ,iSer); %!! 190
    popup.name=gainSer;
    pupup.actual=
    pupup.unit=
    pupup.date =
    pupup.period =
    pupup.weight =
    pupup.news = 
    pupup.impact =
    pupup.yold = y_old;
    pupup.ynew = y_new;
end


% DatesNews = nan(length(OldFcst),2);
% DatesNews(1:4:end,:) = DatesV(iS:iE,:);
% DatesNews(2:4:end,:) = DatesV(iS:iE,:);
% 
% GroupNames = unique(Group)';
% TrueSer = Data(iQ,iSer);

% check whether the new forecats is equal to the old forecast plus the news 
check = y_new-y_old-sum(singlenews,2);

%save the results in structures

