%% Preallocate to store results
YQ_bc = nan(length(StartVintage:EndVintage), length(QQ));
YQ_nc = nan(length(StartVintage:EndVintage), length(QQ));
YQ_fc = nan(length(StartVintage:EndVintage), length(QQ));

YQunc = nan(length(StartVintage:EndVintage), length(QQ));
PSTunc = nan(length(StartVintage:EndVintage), length(YY));

PST_bc = nan(length(StartVintage:EndVintage), length(YY));
PST_nc = nan(length(StartVintage:EndVintage), length(YY));
PST_fc = nan(length(StartVintage:EndVintage), length(YY));

QST_bc = nan(length(StartVintage:EndVintage), length(QQ));
QST_nc = nan(length(StartVintage:EndVintage), length(QQ));
QST_fc = nan(length(StartVintage:EndVintage), length(QQ));

Par_bc = nan(length(StartVintage:EndVintage), 4);
Par_nc = nan(length(StartVintage:EndVintage), 4);
Par_fc = nan(length(StartVintage:EndVintage), 4);

OT_bc = nan(length(StartVintage:EndVintage), 1);
OT_nc = nan(length(StartVintage:EndVintage), 1);
OT_fc = nan(length(StartVintage:EndVintage), 1);

PS_bc = nan(length(StartVintage:EndVintage), 1);
PS_nc = nan(length(StartVintage:EndVintage), 1);
PS_fc = nan(length(StartVintage:EndVintage), 1);

LR_bc = nan(length(StartVintage:EndVintage), 1);
LR_nc = nan(length(StartVintage:EndVintage), 1);
LR_fc = nan(length(StartVintage:EndVintage), 1);

SF_bc = nan(length(StartVintage:EndVintage), 1);
SF_nc = nan(length(StartVintage:EndVintage), 1);
SF_fc = nan(length(StartVintage:EndVintage), 1);

factor_global = nan(floor(length(DatesV)/3), length(name_vintage));
factor_fin = nan(floor(length(DatesV)/3), length(name_vintage));