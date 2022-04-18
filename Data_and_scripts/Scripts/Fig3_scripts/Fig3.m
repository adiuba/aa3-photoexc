clc
clear
close all

run('auxiliary_scripts\ma_fit_623_608_570_neg_coupled_res305') % run fitting
load mat\res                                                   % load fitting result
run('auxiliary_scripts\selectivity_estimation')                % run selectivity estimation
load mat\selectivity                                           % load selectivity estimation workspace
load mat\aa3570MA.mat
load mat\aa3608MA.mat
load mat\aa3623MA.mat
Cashift=load('..\Fig2_scripts\mat\MV18all_res_v2_02042021');

figure('position', [10 50 1300 1300], 'paperpositionmode', 'auto');
axpos{1}=[0.08 0.74 0.27 0.24];
axpos{2}=[0.43 0.74 0.5 0.24];
axpos{3}=[0.08 0.455 0.27 0.24];
axpos{4}=[0.395 0.455 0.27 0.24];
axpos{5}=[0.71 0.455 0.27 0.24];
axpos{6}=[0.08 0.07 0.27 0.34];
axpos{7}=[0.395 0.07 0.27 0.34];
axpos{8}=[0.71 0.07 0.27 0.34];

olive=[51 153 51]/255;
ca=[204 51 0]/255;
ca3=[0 102 0]/255;
darkred=[153 0 51]/255;

ylabels={'-\DeltaA', 'A/[aa_3] (mM^{-1}cm^{-1})',...
    '-\DeltaA','','','Absorption (r.u.)','',''};
xlabels='Wavelength (nm)';
xl=[400 470];
letter='abcdefgh';

police=16; lw=1.5;
ax=cell(length(axpos),1);
for k=1:length(axpos)
    ax{k}=axes('position',axpos{k});
    set(gca,'fontsize',police)
    set(gca,'linewidth',lw)
    set(gca,'ticklength',[0.015 0.015])
    box on
    ylabel(ylabels{k},'fontsize',police)
    xlim(xl)
    pos = get(gca,'Position');
    xpos=0.02; xposnorm=1-xpos/pos(3);
    ypos=0.01; yposnorm=1-ypos/pos(4);
    text(xposnorm,yposnorm,['\it',letter(k)],'fontsize',police,'Units','normalized')
end


%% EASs
axes(ax{1})
line(xl,[0 0],'linewidth',lw,'color','k')
hold on
plot(aa3560MA(:,1),aa3560MA(:,2),'linewidth',lw,'color','b')
plot(aa3600MA(:,1),aa3600MA(:,2),'linewidth',lw,'color',olive)
plot(aa3620MA(:,1),aa3620MA(:,2),'linewidth',lw,'color','r')

yl1=[-0.03 0.08];
ylim(yl1)


text(xl(1)+0.1*diff(xl),yl1(1)+0.9*diff(yl1),'570exc','fontsize',police,'color','b')
text(xl(1)+0.1*diff(xl),yl1(1)+0.75*diff(yl1),'608exc','fontsize',police,'color',olive)
text(xl(1)+0.1*diff(xl),yl1(1)+0.6*diff(yl1),'623exc','fontsize',police,'color','r')

mxps560=sortrows(aa3560MA,2); max560=mean(mxps560(end-1:end,1));
mxps600=sortrows(aa3600MA,2); max600=mean(mxps600(end-1:end,1))-0.3;
mxps620=sortrows(aa3620MA,2); max620=mean(mxps620(end-1:end,1))+0.3;

line([max560 max560],[yl1(1)+0.03*diff(yl1) yl1(1)+0.1*diff(yl1)],'color','b','linewidth',lw)
line([max600 max600],[yl1(1)+0.03*diff(yl1) yl1(1)+0.1*diff(yl1)],'color',olive,'linewidth',lw)
line([max620 max620],[yl1(1)+0.03*diff(yl1) yl1(1)+0.1*diff(yl1)],'color','r','linewidth',lw)

text(0.03,0.08,'max positions\rightarrow','fontsize',police,'Units','normalized')

%% selectivity estimation
axes(ax{2})
hold on
fill([nm; nm(end:-1:1)],[zeros(length(a_Ca),1); a_Ca(end:-1:1)]/max(a_Ca)*a_Ca_abs_max,ca,'EdgeColor','none')
fill([nm; nm(end:-1:1)],[zeros(length(a3_Ca),1); a3_Ca(end:-1:1)],ca3,'EdgeColor','none')

alpha(0.5)
plot(nm,aa3(:,2),'color',0.5*ones(1,3),'linewidth',lw)
xl2=[500 660];
xlim(xl2)
ylim([0 50])
set(gca,'xtick',500:40:700)

text(0.02,0.78,'570exc','fontsize',police,'color','b','Units','normalized')
text(0.02,0.63,'608exc','fontsize',police,'color',olive,'Units','normalized')
text(0.02,0.48,'623exc','fontsize',police,'color','r','Units','normalized')

text(0.02,0.91,'selectivity {\ita}/{\ita}_3','fontsize',police,'color','k','Units','normalized')
text(0.21,0.78,num2str(sel_570_aa3_Ca2,'%.2f'),'fontsize',police,'color','b','Units','normalized')
text(0.21,0.63,num2str(sel_608_aa3_Ca2,'%.2f'),'fontsize',police,'color',olive,'Units','normalized')
text(0.21,0.48,num2str(sel_623_aa3_Ca2,'%.2f'),'fontsize',police,'color','r','Units','normalized')

text(600,6,'{\ita}_3','fontsize',police,'HorizontalAlignment','center')
text(604,20,'{\ita}','fontsize',police,'HorizontalAlignment','center')
text(525,15,'{\itaa}_3','fontsize',police,'HorizontalAlignment','center','color',0.6*ones(1,3))

set(gca,'YTick',0:10:50)

ax_right = axes('Position', ax{2}.Position, 'Color', 'none', 'YAxisLocation', 'Right');
hold on
set(gca,'XTick',[],'ycolor','k')
plot(pulse623_aa3(:,1),pulse623_aa3(:,2),'r','linewidth',lw)
plot(pulse608_aa3(:,1),pulse608_aa3(:,2),'color',olive,'linewidth',lw)
plot(pulse570_aa3_ba3(:,1),pulse570_aa3_ba3(:,2),'b','linewidth',lw)
xlim(xl2)
ylim([0 1])
set(gca,'YTick',0:0.2:1)
set(gca,'YTickLabels',{'0',[],[],[],[],'1',[]})
set(gca,'linewidth',lw)
set(gca,'fontsize',police)
ylabel('Light power (r.u.)','fontsize',police,'color','k','rotation',270,'position',[675 0.5])
nmsel=nm;

%% fitting results

titles{1}='570 exc';
titles{2}='608 exc';
titles{3}='623 exc';


nm=aa3560MA(:,1);
nmlmt=[nm(1) nm(end)];
sp=[aa3560MA(nm>nmlmt(1)&nm<nmlmt(2),2) aa3600MA(nm>nmlmt(1)&nm<nmlmt(2),2) aa3620MA(nm>nmlmt(1)&nm<nmlmt(2),2)];
nm=nm(nm>nmlmt(1)&nm<nmlmt(2)); X=1e7./nm;
msp=max(sp);
msp1=ones(numel(X),1)*msp;
sp=sp./msp1;
X1=X; 
X=25000:-1:20000;
nm=1e7./X;

Npos=res.Npos;
W=res.W;
l=res.res;
Nneg=2;

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

for m=1:3
    
    s1=zeros(numel(X),Npos+Nneg);
    %----collecting positive gaussians
    for k=1:Npos
        s1(:,k)=apos(k,m)*exp(-(X-1e7/ppos(k)).^2/(2*ws(k)^2));
    end
    
    %----collecting negative gaussians
    for k=Npos+1:Npos+Nneg
        f=k-Npos;
        s1(:,k)=aneg(f,m)*exp(-(X-1e7/pneg(f)).^2/(2*wneg(f)^2));
    end
    S1(:,m)=sum(s1,2);
    S{m}=s1;
end


for m=1:3
    axes(ax{m+2})
    hold on
    plot(1e7./X, S{m}(:,[1 3 4]),'--','color',ca3,'linewidth',lw)
    plot(1e7./X, S{m}(:,[2 5 6]),'--','color',ca,'linewidth',lw)
    plot(1e7./X, S{m}(:,7:end),':','color',0.5*ones(1,3),'linewidth',lw)
    plot(1e7./X1, sp(:,m),'color',[0.8 0.8 0.8],'linewidth',6)
    plot(1e7./X, S1(:,m),'color','r','linewidth',3)
    xl3=[410 480]; yl3=[-0.2 1.1];
    set(gca,'Xlim',xl3,'Ylim',yl3)
    box on
    text(470,0.4,'{\ita}','fontsize',16,'color',ca)
    text(470,0.2,'{\ita}_3','color',ca3,'fontsize',16)
    text(xl3(1)+0.07*diff(xl3),yl3(1)+0.9*diff(yl3),titles{m},'fontsize',police)
end

%% "pure" spectral lineshape of the hemes

lb=0.04;

a3=sum(S{1}(:,[1 3 4]),2); ma3=max(a3); m1a3=max(S{1}(:,[1 3 4]),[],1); ksa3=m1a3/ma3;
a3_norm=a3/ma3; a3_norm_comp=S{1}(:,[1 3 4]).*(ones(length(a3),1)*ksa3)./(ones(length(a3),1)*m1a3);

a=sum(S{1}(:,[2 5 6]),2); ma=max(a); m1a=max(S{1}(:,[2 5 6]),[],1); ksa=m1a/ma;
a_norm=a/ma; a_norm_comp=S{1}(:,[2 5 6]).*(ones(length(a),1)*ksa)./(ones(length(a),1)*m1a);

ma=max(a_norm); man=nm(a_norm==ma);
ma3=max(a3_norm); ma3n=nm(a3_norm==ma3);

axes(ax{6})
hold on
plot(Cashift.nmMV18V(1:end-1),Cashift.yfull,'color',0.5*ones(1,3),'linewidth',lw)
plot(1e7./X, a_norm,'color',ca,'linewidth',lw)
plot(1e7./X, a_norm_comp,'--','color',ca,'linewidth',lw)

xl4=[400 480]; yl4=[0 1.2];
set(gca,'Xlim',xl4,'Ylim',yl4)
box on
text(0.05,0.9,'heme {\ita}^{2+}','fontsize',police,'Units','normalized')
text(0.05,0.75,'bleaching','fontsize',police,'color',ca,'Units','normalized')
text(0.05,0.65,'Ca^{2+} shift','fontsize',police,'color',0.5*ones(1,3),'Units','normalized')

line([man man],[ma-lb ma+lb],'color','k','linewidth',1)
line(ppos(2)*ones(1,2),ksa(1)+[-lb +lb],'color','k','linewidth',1)
line(ppos(5)*ones(1,2),ksa(2)+[-lb +lb],'color','k','linewidth',1)
line([ppos(5) 430],ksa(2)+[lb 2*lb],'color','k','linewidth',1)
line(ppos(6)*ones(1,2),ksa(3)+[-lb +lb],'color','k','linewidth',1)
line([ppos(6) 460],ksa(3)+[lb 2*lb],'color','k','linewidth',1)

text(man,ma+2*lb,num2str(round(man)),'fontsize',police,'HorizontalAlignment','center')
text(ppos(2)+3,ksa(1)+2*lb,num2str(round(ppos(2))),'fontsize',police,'HorizontalAlignment','center')
text(429,ksa(2)+3*lb-0.02,num2str(round(ppos(5))),'fontsize',police,'HorizontalAlignment','right')
text(461,ksa(3)+3*lb-0.02,num2str(round(ppos(6))),'fontsize',police,'HorizontalAlignment','left')

[posnm_a, poscm_a, wcm_a, wnm_a, n1nm_a, n2nm_a, hm_a] = fwhm_gauss(ppos([2 5 6])',wpos([2 4 4])',ksa','nm');
probe_a=[posnm_a poscm_a wcm_a wnm_a];
center_cm_a=mean(poscm_a(2:3));
center_nm_a=mean(posnm_a(2:3));

ngauss=1;
line([n1nm_a(ngauss) n2nm_a(ngauss)],hm_a(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
line([419 417],[hm_a(ngauss)+0.01 hm_a(ngauss)+0.07],'color','k','linewidth',1)
text(417,hm_a(ngauss)+0.1,[(num2str(wnm_a(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','right')

ngauss=3;
line([n1nm_a(ngauss) n2nm_a(ngauss)],hm_a(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
line([452 460],[hm_a(ngauss)+0.01 hm_a(ngauss)+0.07],'color','k','linewidth',1)
text(460,hm_a(ngauss)+0.1,[(num2str(wnm_a(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','left')

set(gca,'YTick',0:0.2:1)

axes(ax{7})
hold on
plot(1e7./X, a3_norm,'color',ca3,'linewidth',lw)
plot(1e7./X, a3_norm_comp,'--','color',ca3,'linewidth',lw)

xl4=[400 480]; yl4=[0 1.2];
set(gca,'Xlim',xl4,'Ylim',yl4)
box on
text(0.05,0.9,'heme {\ita}_3^{2+}','fontsize',police,'Units','normalized')

line([ma3n ma3n],[ma3-lb ma3+lb],'color','k','linewidth',1)
line(ppos(1)*ones(1,2),ksa3(1)+[-lb +lb],'color','k','linewidth',1)
line(ppos(3)*ones(1,2),ksa3(2)+[-lb +lb],'color','k','linewidth',1)
line(ppos(4)*ones(1,2),ksa3(3)+[-lb +lb],'color','k','linewidth',1)
line([ppos(4) 454],ksa3(3)+[lb 2*lb],'color','k','linewidth',1)

text(ma3n,ma3+2*lb,num2str(round(ma3n)),'fontsize',police,'HorizontalAlignment','center')
text(ppos(1)-1,ksa3(1)+2*lb,num2str(round(ppos(1))),'fontsize',police,'HorizontalAlignment','center')
text(ppos(3)+1,ksa3(2)+2*lb,num2str(round(ppos(3))),'fontsize',police,'HorizontalAlignment','center')
text(455,ksa3(3)+3*lb-0.02,num2str(round(ppos(4))),'fontsize',police,'HorizontalAlignment','left')

[posnm_a3, poscm_a3, wcm_a3, wnm_a3, n1nm_a3, n2nm_a3, hm_a3] = fwhm_gauss(ppos([1 3 4])',wpos([1 3 3])',ksa3','nm');
probe_a3=[posnm_a3 poscm_a3 wcm_a3 wnm_a3];
center_cm_a3=mean(poscm_a3(2:3));
center_nm_a3=mean(posnm_a3(2:3));

ngauss=1;
line([n1nm_a3(ngauss) n2nm_a3(ngauss)],hm_a3(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
text(423,hm_a3(ngauss)+0.04,[(num2str(wnm_a3(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','center')

ngauss=3;
line([n1nm_a3(ngauss) n2nm_a3(ngauss)],hm_a3(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
line([444 455],[hm_a3(ngauss)+0.01 hm_a3(ngauss)+0.07],'color','k','linewidth',1)
text(455,hm_a3(ngauss)+0.1,[(num2str(wnm_a3(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','left')
set(gca,'YTick',0:0.2:1)

xlabel(xlabels,'fontsize',police)

%% Soret band reconstruction
load ..\Fig2_scripts\mat\dith_ruth_ready

% %estimated selectivity at 570exc as a reference
% a3=sum(S{1}(:,[1 3 4]),2); a=sum(S{1}(:,[2 5 6]),2); Soret=a/sel_570_aa3_Ca2+a3;
% ma=max(a)/sel_570_aa3_Ca2/max(Soret);

%estimated selectivity at 608exc as a reference
a3=sum(S{2}(:,[1 3 4]),2); a=sum(S{2}(:,[2 5 6]),2); Soret=a/sel_608_aa3_Ca2+a3;
ma=max(a)/sel_608_aa3_Ca2/max(Soret);

% %estimated selectivity at 623exc as a reference
% a3=sum(S{3}(:,[1 3 4]),2); a=sum(S{3}(:,[2 5 6]),2); Soret=a/sel_623_aa3_Ca2+a3;
% ma=max(a)/sel_623_aa3_Ca2/max(Soret);

aSoret=a/max(a)*ma;
ma3=max(a3)/max(Soret); a3Soret=a3/max(a3)*ma3;
Soret=Soret/max(Soret);
ditnorm=dit_ready/max(dit_ready);

axes(ax{8})
hold on

plot(nm,a3Soret,'color',ca3,'linewidth',lw)
plot(nm,aSoret,'color',ca,'linewidth',lw)
plot(nm,Soret,'r','linewidth',lw)
plot(nmdit,ditnorm,'color',0.5*ones(1,3),'linewidth',lw)
set(gca,'Xlim',xl4,'Ylim',yl4)
set(gca,'YTick',0:0.2:1)

text(0.05,0.9,'{\ita}^{2+}{\ita}_3^{2+}','fontsize',police,'color',0.5*ones(1,3),'Units','normalized')
text(0.05,0.77,'heme {\ita}^{2+}','fontsize',police,'color',ca,'Units','normalized')
text(0.05,0.64,'heme {\ita}_3^{2+}','fontsize',police,'color',ca3,'Units','normalized')
text(0.05,0.5,'{\ita}^{2+}+{\ita}_3^{2+}','fontsize',police,'color','r','Units','normalized')

a3_area=-trapz(X,a3Soret); a_area=-trapz(X,aSoret);
a_to_a3=a_area/a3_area;

realsel560=trapz(X,sum(S{1}(:,[2 5 6]),2))/trapz(X,sum(S{1}(:,[1 3 4]),2))/a_to_a3
realsel600=trapz(X,sum(S{2}(:,[2 5 6]),2))/trapz(X,sum(S{2}(:,[1 3 4]),2))/a_to_a3
realsel620=trapz(X,sum(S{3}(:,[2 5 6]),2))/trapz(X,sum(S{3}(:,[1 3 4]),2))/a_to_a3

print('Fig3','-dpng')

%---------- save data in csv format ----------------------------
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3a\1_a2+a32+_570exc.txt',aa3560MA)
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3a\2_a2+a32+_608exc.txt',aa3600MA)
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3a\3_a2+a32+_623exc.txt',aa3620MA)

csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3b\1_a2+a32+_alpha_beta_region.txt',[nmsel aa3(:,2)])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3b\2_a2+.txt',[nmsel a_Ca/max(a_Ca)*a_Ca_abs_max])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3b\3_a32+.txt',[nmsel a3_Ca])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3b\4_pulse_570exc.txt',pulse570_aa3_ba3)
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3b\5_pulse_608exc.txt',pulse608_aa3)
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3b\6_pulse_623exc.txt',pulse623_aa3)

csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3c\1_a2+a32+_exc570_norm.txt',[1e7./X1 sp(:,1)])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3c\2_a2+a32+_exc570_fit.txt',[1e7./X', S1(:,1)])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3c\3_a2+_exc570_gaussians.txt',[1e7./X' S{1}(:,[2 5 6])])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3c\4_a32+_exc570_gaussians.txt',[1e7./X' S{1}(:,[1 3 4])])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3c\5_a2+a32+_exc570_negative_gaussians.txt',[1e7./X' S{1}(:,7:end)])

csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3d\1_a2+a32+_exc608_norm.txt',[1e7./X1 sp(:,2)])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3d\2_a2+a32+_exc608_fit.txt',[1e7./X', S1(:,2)])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3d\3_a2+_exc608_gaussians.txt',[1e7./X' S{2}(:,[2 5 6])])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3d\4_a32+_exc608_gaussians.txt',[1e7./X' S{2}(:,[1 3 4])])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3d\5_a2+a32+_exc608_negative_gaussians.txt',[1e7./X' S{2}(:,7:end)])

csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3e\1_a2+a32+_exc623_norm.txt',[1e7./X1 sp(:,3)])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3e\2_a2+a32+_exc623_fit.txt',[1e7./X', S1(:,3)])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3e\3_a2+_exc623_gaussians.txt',[1e7./X' S{3}(:,[2 5 6])])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3e\4_a32+_exc623_gaussians.txt',[1e7./X' S{3}(:,[1 3 4])])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3e\5_a2+a32+_exc623_negative_gaussians.txt',[1e7./X' S{3}(:,7:end)])

csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3f\1_a2+_bleaching_sum.txt',[1e7./X' a_norm])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3f\2_a2+_bleaching_gaussians.txt',[1e7./X' a_norm_comp])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3f\3_a2+_Ca2+_reconstruction.txt',[Cashift.nmMV18V(1:end-1) Cashift.yfull])

csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3g\1_a32+_bleaching_sum.txt',[1e7./X' a3_norm])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3g\2_a32+_bleaching_gaussians.txt',[1e7./X' a3_norm_comp])

csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3h\1_a2+_bleaching_sum_scaled.txt',[nm' aSoret])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3h\2_a32+_bleaching_sum_scaled.txt',[nm' a3Soret])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3h\3_a2+_scaled+a32+_scaled.txt',[nm' Soret])
csvwrite('..\..\Data\Figures\Fig3_data\Panel_Fig3h\4_a2+a32+_norm.txt',[nmdit ditnorm])
