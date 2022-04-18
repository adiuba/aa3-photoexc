clc
clear
close all

load mat\FR_MV2018
run('auxiliary_scripts\antiderMV18all_v2') % take integral over a2+a33+-CN +/- Ca2+ difference and decompose the integral
run('auxiliary_scripts\BH_to_dith')
figure('position', [10 50 1100 1300], 'paperpositionmode', 'auto');
axpos{1}=[0.08 0.74 0.4 0.25];
axpos{2}=[0.58 0.74 0.40 0.25];
axpos{3}=[0.08 0.39 0.90 0.30];
axpos{4}=[0.08 0.09 0.40 0.25];
axpos{5}=[0.26 0.18 0.13 0.13];
axpos{6}=[0.58 0.09 0.40 0.25];

olive=[51 153 51]/255;
ca=[204 51 0]/255;
ca3=[0 102 0]/255;

ylabels={'A/[aa_3] (mM^{-1}cm^{-1})', '\DeltaA/[aa_3] (mM^{-1}cm^{-1})',...
    'Absorption (r.u.)','A/[aa_3] (mM^{-1}cm^{-1})','A/[aa_3]','normalized MCD (r.u.)'};
xlabels='Wavelength (nm)';
xl=[400 630];
letter='abcdefgh';

police=16; lw=1.5;
ax=cell(length(axpos),1);
k_let=[1 2 3 4 4 5];
for k=1:length(axpos)
    ax{k}=axes('position',axpos{k});
    set(gca,'fontsize',police)
    set(gca,'linewidth',lw)
    set(gca,'ticklength',[0.015 0.015])
    box on

%     xlim(xl)
    pos = get(gca,'Position');
    xpos=0.02; xposnorm=1-xpos/pos(3);
    ypos=0.01; yposnorm=1-ypos/pos(4);
    if k~=5
            ylabel(ylabels{k},'fontsize',police)
        text(xposnorm,yposnorm,['\it',letter(k_let(k))],'fontsize',police,'Units','normalized')
    end
    
    pos = get(gca,'Position');
    desired_length = 0.01; %cm
    normalized_length_h = desired_length/sqrt(sum([pos(3) pos(4)]).^2);
    normalized_length_w = desired_length/pos(3);
    set(gca,'TickLength', [normalized_length_h, normalized_length_w])
end

%% Parent spectra
a=y;

axes(ax{1})
hold on
d=10;
plot(nm_MV18,(absMV18-absMV18(1))/cMV18/1e3+absFR18(1)/cFR18/1e3,'k','linewidth',lw)
plot(nm_FR18,absFR18/cFR18/1e3,'r','linewidth',lw)
xlim([400 650])
ylim([0 270])
set(gca,'xtick',400:50:650)
line([430.3 430.3],[-d d]+195,'color','k','linewidth',1)
line([441.9 441.9],[-d d]+200.8,'color','k','linewidth',1)
line([441.9 455],[d+200.8 d+200.8+10],'color','k','linewidth',1)
line([443.9 443.9],[-d d]+246.7,'color','k','linewidth',1)
line([604.8 604.8],[-d d]+60.6,'color','k','linewidth',1)
text(408,2*d+195,'430','color','k','fontsize',police)
text(456,2*d+201,'442','color','k','fontsize',police)
text(448,d+247,'444','color','k','fontsize',police)
text(602,2*d+60.6,'605','color','k','fontsize',police)

text(530,200,'{\ita}^{2+}{\ita}_3^{3+}-CN','color','k','fontsize',police)
text(530,160,'{\ita}^{2+}{\ita}_3^{2+}','color','r','fontsize',police)

axes(ax{2})
hold on
d=0.5;
line([400 630],[0 0],'color','k')
plot(nm_FR18,CashiftFR18,'r','linewidth',lw)
plot(nm_MV18,CashiftMV18,'k','linewidth',lw)
% plot(nm_MV18,baseline,'color',0.5*ones(1,3),'linewidth',lw+1,'linestyle',':')
plot(nm_MV18,baseline,'color',olive,'linewidth',lw+1,'linestyle',':')
set(gca,'ytick',-5:2.5:5,'xtick',400:50:650)
box on
xlim([400 650])
ylim([-7 7])
line([438.3 438.3],[-d d]-4.2,'color','k','linewidth',1)
line([460 451],[-1 -.1],'color','k','linewidth',1)
line([454.6 454.6],[-d d]+4.6,'color','k','linewidth',1)
line([598 598],[-d d]-2.1,'color','k','linewidth',1)
line([616 607],[-1 -.1],'color','k','linewidth',1)
line([611.4 611.4],[-d d]+2.8,'color','k','linewidth',1)
text(425,-5.2,'438','color','k','fontsize',police)
text(446,5.5,'460','color','k','fontsize',police)
text(448,-1.5,'449','color','k','fontsize',police)
text(585,-3,'598','color','k','fontsize',police)
text(598,3.7,'611','color','k','fontsize',police)
text(615,-1.5,'605','color','k','fontsize',police)

%% Ca shift reconstruction
nmCa=nmMV18V(1:end-1);
axes(ax{3})
hold on
plot(nmCa,yfull,'k','linewidth',lw)
plot(nmCa,[ss1 ss2 ss3 ss4 ss5 ss6 ss7],'r--','linewidth',lw)
plot(nmCa,SS,'r','linewidth',lw)
xlim([400 630])
ylim([0 1.2])
set(gca,'xtick',400:50:630)
lb=0.035;
line(1e7./[n1 n1],[h1-lb h1+lb],'color','k','linewidth',1)
line(1e7./[n2 n2],[h2-lb h2+lb],'color','k','linewidth',1)
line(1e7./[n3 n3],[h3-lb h3+lb],'color','k','linewidth',1)
line(1e7./[n4 n4],[h4-lb h4+lb],'color','k','linewidth',1)
line(1e7./[n5 n5],[h5-lb h5+lb],'color','k','linewidth',1)
line(1e7./[n6 n6],[h6-lb h6+lb],'color','k','linewidth',1)
line(1e7./[n7 n7],[h7-lb h7+lb],'color','k','linewidth',1)

line([435 1e7/n2],[h2+2*lb h2+lb],'color','k','linewidth',1)
line([1e7/n3 460],[h3+lb h3+2*lb],'color','k','linewidth',1)
line([586 1e7/n6],[h6+3*lb h6+lb],'color','k','linewidth',1)
line([1e7/n7 612],[h7+lb h7+2*lb],'color','k','linewidth',1)

line([605 605],[-lb lb]+0.41,'color','k','linewidth',1)
line([447.6 447.6],[-lb lb]+1,'color','k','linewidth',1)

text(1e7/n1-5,h1+0.08,num2str(round(1e7/n1)),'fontsize',police)
text(1e7/n2-17.5,h2+0.1,num2str(round(1e7/n2)),'fontsize',police)
text(1e7/n3+10,h3+0.08,num2str(round(1e7/n3)),'fontsize',police)
text(1e7/n4-5,h4+0.08,num2str(round(1e7/n4)),'fontsize',police)
text(1e7/n5-5,h5+0.08,num2str(round(1e7/n5)),'fontsize',police)
text(1e7/n6-17.5,h6+0.12,num2str(round(1e7/n6)),'fontsize',police)
text(1e7/n7+7,h7+0.08,num2str(round(1e7/n7)),'fontsize',police)

text(443,1.08,'447.5','fontsize',police)
text(595,0.47,'605','fontsize',police)

line([480 490],[1.01 1.01],'linewidth',lw,'Color','k')
line([480 490],[0.866 0.866],'linewidth',lw,'Color','r','linestyle','--')
line([480 490],[0.734 0.734],'linewidth',lw,'Color','r')
text(0.4,0.85,'Integral over {\ita}^{2+}{\ita}_3^{3+}-CN +/- Ca^{2+} difference spectrum','fontsize',police,'Units','normalized','Color','k')
text(0.4,0.73,'gaussians','fontsize',police,'Units','normalized','Color','r')
text(0.4,0.62,'sum','fontsize',police,'Units','normalized','Color','r')

[posnm, poscm, wcm, wnm, n1nm, n2nm, hm] = fwhm_gauss(res([1 4 7 9 12 15 18])',res([2 5 5 10 13 16 16])',res([3 6 8 11 14 17 19])','nm');
probe=[posnm poscm wcm wnm];

ngauss=3;
line([n1nm(ngauss) n2nm(ngauss)]+0.6,hm(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
text(460,hm(ngauss)+0.1,[(num2str(wnm(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','left')
line([452 460],[hm(ngauss)+0.01 hm(ngauss)+0.07],'color','k','linewidth',1)

ngauss=7;
line([n1nm(ngauss) n2nm(ngauss)],hm(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
line([605 614],[hm(ngauss)+0.01 hm(ngauss)+0.1],'color','k','linewidth',1)
text(614,hm(ngauss)+0.135,[(num2str(wnm(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','left')

%% BH spectrum

axes(ax{4})
hold on
plot(BH.nm,BH.BH_ready,'r','linewidth',lw)
xlim([350 650])
ylim([0 270])

lb=10;
line([444 444],218+[-lb lb],'linewidth',1,'Color','k')
line([420 420],112.6+[-lb lb],'linewidth',1,'Color','k')
line([605 605],49.5+[-lb lb],'linewidth',1,'Color','k')
text(444,lb+232,'444','color','k','fontsize',police,'HorizontalAlignment','center')
text(410,lb+126,'420','color','k','fontsize',police,'HorizontalAlignment','center')
text(615,lb+63,'605','color','k','fontsize',police,'HorizontalAlignment','center')

axes(ax{5})
hold on
plot(BH.nm,BH.BH_ready,'r','linewidth',lw)
xlim([360 420])
ylim([70 115])
set(gca,'fontsize',police-3)
set(gca,'ytick',70:20:110)

%% MCD

axes(ax{6})
hold on
run('auxiliary_scripts\MCDfit_manual')
nm=MCDdata(:,1);
S=MCDdata(:,2);
s1=MCDdata(:,3);
s2=MCDdata(:,4);
s3=MCDdata(:,5);
MCD=MCDdata(:,6);

plot(nm,[s1 s2 s3],'--r','linewidth',lw)
plot(nm,MCD,'-k','linewidth',lw)
plot(nm,S,'-r','linewidth',lw)

lb=0.06;

box on
xlim([550 650])
ylim([-0.9 0.9])
set(gca,'xtick',400:25:650,'ytick',-0.8:0.4:0.8)

line(xl,[0 0],'linewidth',lw,'Color','k')

line([574 574],[-lb lb]-0.25,'linewidth',1,'Color','k')
line([594 594],[-lb lb]-0.66,'linewidth',1,'Color','k')
line([605.5 605.5],[-lb lb]+0.6589,'linewidth',1,'Color','k')
line([609 609],[-lb lb]+0.4618,'linewidth',1,'Color','k')
line([609 612],[lb 0.1]+0.4618,'linewidth',1,'Color','k')
line([590 590],[-lb lb]-0.48,'linewidth',1,'Color','k')
line([590 587],[-lb -0.1]-0.48,'linewidth',1,'Color','k')

text(574,-0.42,'574','HorizontalAlignment','center','fontsize',police)
text(596,-0.80,'594','HorizontalAlignment','right','fontsize',police)
text(606,0.78,'605.5','HorizontalAlignment','left','fontsize',police)
text(613,0.56,'609','HorizontalAlignment','left','fontsize',police)
text(586.5,-0.6,'590','HorizontalAlignment','right','fontsize',police)


line([555 565],[0.75 0.75],'linewidth',lw,'Color','k')
line([555 565],[0.48 0.48],'linewidth',lw,'Color','r','linestyle','--')
line([555 565],[0.25 0.25],'linewidth',lw,'Color','r')

text(0.17,0.9,'{\ita}^{2+}{\ita}_3^{2+}','fontsize',police,'Units','normalized','Color','k')
text(0.17,0.77,'gaussians','fontsize',police,'Units','normalized','Color','r')
text(0.17,0.64,'sum','fontsize',police,'Units','normalized','Color','r')

text(-0.3,-0.25,xlabels,'fontsize',police,'Units','normalized')


print('Fig2','-dpng')

%---------- save data in csv format ----------------------------
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2a\1_a2+a32+_ascorbate_RuAm.txt',[nm_FR18 absFR18/cFR18/1e3])
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2a\2_a2+a33+_CN_ascorbate_TMPD.txt',[nm_MV18 (absMV18-absMV18(1))/cMV18/1e3+absFR18(1)/cFR18/1e3])

csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2b\1_a2+a32+_ascorbate_RuAm_diff.txt',[nm_FR18 CashiftFR18])
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2b\2_a2+a33+_CN_ascorbate_TMPD_diff.txt',[nm_MV18 CashiftMV18])
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2b\3_baseline.txt',[nm_MV18 baseline])

csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2c\1_a2+_Ca2+_recontruction_integral.txt',[nmCa yfull])
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2c\2_a2+_Ca2+_recontruction_fit.txt',[nmCa SS])
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2c\3_a2+_Ca2+_recontruction_gaussians.txt',[nmCa ss1 ss2 ss3 ss4 ss5 ss6 ss7])

csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2d\1_a2+a32+_borohydride.txt',[BH.nm BH.BH_ready])

csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2e\1_a2+a32+_ascorbate_RuAm_MCD.txt',[nm MCD])
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2e\2_a2+a32+_ascorbate_RuAm_MCD_fit.txt',[nm S])
csvwrite('..\..\Data\Figures\Fig2_data\Panel_Fig2e\3_a2+a32+_ascorbate_RuAm_MCD_gaussians.txt',[nm s1 s2 s3])