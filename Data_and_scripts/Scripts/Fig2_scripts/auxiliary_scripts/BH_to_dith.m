% clear
% close all
% clc

global nm_y dit_y BH_y ox_y
load ..\mat\BHdata
load ..\mat\dith_ruth_ready

% nmlow=400; nmhigh=500; %Soret
% nmlow=500; nmhigh=700; %Vis
nmlow=400; nmhigh=700; %all
dit1=dit_ready(end:-1:1);
nmdit1=nmdit(end:-1:1);
nm_y=nm(nm>nmlow&nm<nmhigh);
BH_y=BH(nm>nmlow&nm<nmhigh);
ox_y=ox(nm>nmlow&nm<nmhigh);

dit_y=new_grid(nmdit1,nm_y,dit1);
% dit_y=dit_y-diff([dit1(end) dit_y(end)]);
% BH1=BH1+max(BH)-max(BH1);
dit_y=dit_y+max(dit1(nmdit1>nmlow&nmdit1<nmhigh))-max(dit_y);

% plot(nmdit1,dit1), hold on, plot(nm_y,dit_y)
% break

%l0 : linepar1 linepar2 c kox

lb=[-0.5 -0.5 2 0.1];
ub=[0.5 0.5 5 0.4];
l0=lb+rand(1,numel(lb)).*(ub-lb);
[lres,errl,exitflag]=fmincon(@BH_to_dith_fit,l0,[],[],[],[],lb,ub);

x00 = [350 700];
y00 = [lres(1) lres(2)];
c00 = [[1; 1]  x00(:)]\y00(:);                        % Calculate Parameter Vector
slope_m = c00(2);
intercept_b = c00(1);
c=lres(3)*1e-3;
kox=lres(4);
c=c/(1-kox);
BH_ready=(BH-(nm*slope_m + intercept_b)-kox*ox)/(1-kox);
BH_ready=BH_ready/c;

nmdit_ready2=nm(nm>=370&nm<700);
dit_ready2=new_grid(nmdit1,nm(nm>=370&nm<700),dit1);
dit_ready2=dit_ready2+max(dit1)-max(dit_ready2);

plot(nm,BH_ready,'r','linewidth',1.5)
hold on
plot(nmdit,dit_ready,'b','linewidth',1.5)
% plot(nmdit_ready2,dit_ready2,'g','linewidth',1.5)

BH.nm=nm;
BH.nmdit=nmdit;
BH.BH_ready=BH_ready;
BH.dit_ready=dit_ready;
BH.lres=lres;
BH.slope_m=slope_m;
BH.intercept_b=intercept_b;
BH.c=c;
BH.kox=kox;
BH.nmdit1=nmdit1;

% save BH_dith_ready nm nmdit BH_ready dit_ready lres slope_m intercept_b c kox nmdit1