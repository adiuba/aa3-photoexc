clear
close all
clc

load mat\Liao

police=16; lw=1.5;

figure('position', [10 50 800 800], 'paperpositionmode', 'auto');

plot(monohisA(:,1),monohisA(:,2),'color','k','linewidth',lw)
hold on
plot(a3redox(:,1),a3redox(:,2),'color','r','linewidth',lw)

set(gca,'fontsize',police)
set(gca,'linewidth',lw)
set(gca,'ticklength',[0.015 0.015])
box on
line([500 700],[0 0],'linestyle','--','linewidth',lw,'color','k')

yl=get(gca,'ylim'); d=0.02*diff(yl); dd=0.02*diff(yl);
text(0.65,0.85,'1,2 diMe-Im heme A','fontsize',police,'units','normalized','color','k')
text(0.65,0.75,'bovine CcO heme {\ita}_3','fontsize',police,'units','normalized','color','r')

[m1,i1]=max(monohisA(:,2));
[m2,i2]=max(a3redox(:,2));

line(monohisA(i1,1)*ones(1,2),m1+[-d d],'linewidth',lw,'color','k')
line(a3redox(i2,1)*ones(1,2),m2+[-d d],'linewidth',lw,'color','k')

text(monohisA(i1,1),m1+d+dd,num2str(round(monohisA(i1,1))),'fontsize',police,'HorizontalAlignment','center')
text(a3redox(i2,1)+5,m2+d+dd,num2str(round(a3redox(i2,1))),'fontsize',police,'HorizontalAlignment','center')

xlabel('Wavelength (nm)')
ylabel('\DeltaA (mM^{-1}cm^{-1})')
% ,'position',[365 -40]
% ,'position',[242 560]

print('-dpng','-r400','FigS6')

%---------- save data in csv format ----------------------------
csvwrite('..\..\Data\Figures\FigS6_data\1_bov_a3_redox.txt',a3redox)
csvwrite('..\..\Data\Figures\FigS6_data\2_1_2_diMeIm_heme_A.txt',monohisA)
