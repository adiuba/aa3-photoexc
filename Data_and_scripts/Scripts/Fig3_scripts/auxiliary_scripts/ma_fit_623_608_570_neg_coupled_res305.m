% clear global
% close all
% clear
global X drw Npos Nneg W err_all sp dit nmdit
err_all=[];
%---constants-----------
load ..\mat\param              %load parameters
load ..\mat\aa3570MA.mat
load ..\mat\aa3608MA.mat
load ..\mat\aa3623MA.mat
load ..\mat\dit_data
Npos=6; % number of positive gaussians
W=[1 1 2 2]; % number of half-width
Nneg=2; % number of negative gaussians
%-------------------------------
nm=aa3560MA(:,1);
nmlmt=[409 485];
sp=[aa3560MA(nm>nmlmt(1)&nm<nmlmt(2),2) aa3600MA(nm>nmlmt(1)&nm<nmlmt(2),2) aa3620MA(nm>nmlmt(1)&nm<nmlmt(2),2)];
nm=nm(nm>nmlmt(1)&nm<nmlmt(2)); X=1e7./nm;
msp=max(sp);
msp1=ones(numel(X),1)*msp;
sp=sp./msp1;

%structure of l0:
%1. parameters for positive gaussians
%   a. Positions (Npos)
%   b. Widths (numel(W))
%   c. Amplitudes (Npos*3)
%
%2. parameters for negative gaussians
%   a. Positions (Nneg*3)
%   b. Widths (Nneg)
%   c. Amplitudes (Nneg)

drw=0;
lapos=param.LB.lapos; lpneg=param.LB.lpneg; lwneg=param.LB.lwneg; laneg=param.LB.laneg;
uapos=param.UB.uapos; upneg=param.UB.upneg; uwneg=param.UB.uwneg; uaneg=param.UB.uaneg;
lppos=param.LB.lppos; lwpos=param.LB.lwpos;
uppos=param.UB.uppos; uwpos=param.UB.uwpos;
lb=[lppos lwpos lapos lpneg lwneg laneg];
ub=[uppos uwpos uapos upneg uwneg uaneg];
l0=(lb+ub)/2.*(1+(rand(1,numel(lb))-0.5)*0.5);

OPTIONS=param.options;

[lres,errl,exitflag]=fmincon(@aa3_anygauss_fit_ma_623_608_570_neg_coupled,l0,[],[],[],[],lb,ub,@fixed_ratio_aa3_ma_623_608_570,OPTIONS);

drw=0;
aa3_anygauss_fit_ma_623_608_570_neg_coupled(lres)

l=lres;

ppos=l(1:Npos); l(1:Npos)=[];
wpos=l(1:numel(W)); l(1:numel(W))=[];
apos=l(1:Npos*3); l(1:Npos*3)=[];
apos=reshape(apos,Npos,3);

pneg=l(1:sum(Nneg)); l(1:sum(Nneg))=[];
wneg=l(1:sum(Nneg)); l(1:sum(Nneg))=[];
aneg=l(1:sum(Nneg)*3);
aneg=reshape(aneg,sum(Nneg),3);

ws=[];
for k=1:numel(W)
    ws=[ws wpos(k)*ones(1,W(k))];
end

res.LB=param.LB;
res.UB=param.UB;
res.options=param.options;
res.l0=l0;
res.errl=errl;
res.res=lres;
res.exitflag=exitflag;
res.Npos=Npos;
res.W=W;
save ..\mat\res res