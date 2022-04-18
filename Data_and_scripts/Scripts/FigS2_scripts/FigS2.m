clc
clear
close all

global ba3x ba3y

ba3=dlmread('..\Fig3_scripts\mat\ba3_metallomics_static_abs.txt');
load mat\ba3MA560
load ..\Fig3_scripts\mat\pulses
load mat\SQR_c
load mat\ba3620MA
load mat\SMR

figure('position', [10 50 1100 900], 'paperpositionmode', 'auto');
axpos{1}=[0.08 0.56 0.38 0.42];
axpos{2}=[0.53 0.56 0.38 0.42];
axpos{3}=[0.08 0.09 0.38 0.42];
axpos{4}=[0.53 0.09 0.38 0.42];

olive=[51 153 51]/255;
ca=[204 51 0]/255;
ca3=[0 102 0]/255;
csqr=[153 102 0]/255;

ylabels={'A/[ba_3] or A/[SMR] (mM^{-1}cm^{-1})', 'A/[ba_3] or A/[SMR] (mM^{-1}cm^{-1})',...
    '-\DeltaA',''};
xlabels='Wavelength (nm)';
xl=[400 630];
letter='abcdefgh';

police=16; lw=1.5;
ax=cell(length(axpos),1);
for k=1:length(axpos)
    ax{k}=axes('position',axpos{k});
    set(gca,'fontsize',police)
    set(gca,'linewidth',lw)
    set(gca,'ticklength',[0.015 0.015])
    box on

    pos = get(gca,'Position');
    xpos=0.02; xposnorm=1-xpos/pos(3);
    ypos=0.01; yposnorm=1-ypos/pos(4);
    if k~=6
            ylabel(ylabels{k},'fontsize',police)
        text(xposnorm,yposnorm,['\it',letter(k)],'fontsize',police,'Units','normalized')
    end
    
    pos = get(gca,'Position');
    desired_length = 0.01; %cm
    normalized_length_h = desired_length/sqrt(sum([pos(3) pos(4)]).^2);
    normalized_length_w = desired_length/pos(3);
    set(gca,'TickLength', [normalized_length_h, normalized_length_w])
    hold on
end
text(-0.3,-0.14,xlabels,'fontsize',police,'Units','normalized')

%% ba3 absorption spectrum

axes(ax{1})
plot(ba3(:,1),ba3(:,2),'k','linewidth',lw)
plot(SQR_c_nm,SQR_c,'k--','linewidth',lw)
xlim([400 650])
ylim([0 320])
lb=8;
line([425.75 425.75],218+[-lb lb],'color','k','linewidth',1)
line([425.75 441],218+[lb lb+10],'color','k','linewidth',1)
line([443 443],109.3+[-lb lb],'color','k','linewidth',1)
line([558.5 558.5],38.07+[-lb lb]-2,'color','k','linewidth',1)
line([558.5 545],38.07+[lb lb+10]-2,'color','k','linewidth',1)
line([613 613],15.76+[-lb lb],'color','k','linewidth',1)

text(444,218+lb+15,'426','fontsize',police,'HorizontalAlignment','left')
text(452,109.3+lb+7,'443','fontsize',police,'HorizontalAlignment','center')
text(545,38.07+lb+15,'559','fontsize',police,'HorizontalAlignment','right')
text(613,15.76+lb+9,'613','fontsize',police,'HorizontalAlignment','center')


line([520 544],280*ones(1,2),'color','k','linewidth',lw,'linestyle','-')
line([520 544],250*ones(1,2),'color','k','linewidth',lw,'linestyle','--')
text(550,280,'{\itb}^{2+}{\ita}_3^{2+}','fontsize',police,'HorizontalAlignment','left')
text(550,250,'SMR','fontsize',police,'HorizontalAlignment','left')

%% ba3 selectivity illustration

kpulse=1;
ksqr=282;
ba3vis=ba3(ba3(:,1)>500&ba3(:,1)<680,:);
% bamp=25; x=1e7./ba3vis(:,1);
% b=bamp*exp(-(1e7/558.5-x).^2/(2*300^2)); % Gauss in nm
% a3_ba3=ba3vis(:,2)-b;

axes(ax{2})
plot(ax{2},ba3vis(:,1),ba3vis(:,2),'k','linewidth',lw)
plot(SQR_c_nm,SQR_c,'k--','linewidth',lw)
set(gca,'YTick',0:20:80)
xlim([535 660])
ylim([0 60])

ax_right = axes('Position', ax{2}.Position, 'Color', 'none', 'YAxisLocation', 'Right');
hold on
set(gca,'XTick',[],'ycolor','k')
plot(ax_right,pulse570_aa3_ba3(:,1),pulse570_aa3_ba3(:,2)*kpulse,'b','linewidth',lw)
plot(ax_right,pulse623_ba3(:,1),pulse623_ba3(:,2)*kpulse,'r','linewidth',lw)

xlim([535 660])
ylim([0 1.2])
set(gca,'linewidth',lw)
set(gca,'fontsize',police)
set(gca,'YTick',0:0.2:1)
set(gca,'YTickLabels',{'0',[],[],[],[],'1',[]})
ylabel('Light power (r.u.)','fontsize',police,'color','k','rotation',270,'position',[675 0.6])
text(575,1.05,'570exc','fontsize',police,'color','b','HorizontalAlignment','center')
text(625,1.05,'623exc','fontsize',police,'color','r','HorizontalAlignment','center')
text(600,0.4,'{\itb}{\ita}_3','fontsize',police,'HorizontalAlignment','center')
text(630,0.05,'SMR','fontsize',police,'HorizontalAlignment','center')
%% 570exc bleaching spectrum

axes(ax{3})
plot(nm_ba3MA560,ba3MA560,'k','linewidth',lw)
line([405 470],[0 0],'color','k','linewidth',1)
line([443 443],[0.008 0.016],'color','k','linewidth',1)
line([426 426],[0.082 0.09],'color','k','linewidth',1)

xlim([405 470])
ylim([-0.03 0.1])

text(443.5,0.0186,'443','fontsize',police)
text(426.5,0.092,'426','fontsize',police)
text(0.05,0.9,'570exc','fontsize',police,'Units','normalized')


%% 623exc bleaching spectrum
axes(ax{4})
hold on
nm_ba3MA620=ba3620MA(:,1);
ba3x=1e7./ba3620MA(:,1); ba3y=ba3620MA(:,2)/max(ba3620MA(:,2));
lppos=[410 434 438];
lwpos=[350 250];
lapos=ones(1,3)*0;

lpneg=[370 460];
lwneg=ones(1,2)*200;
laneg=[-0.5 -0.5];

uppos=[427 438 446];
uwpos=[700 400];
uapos=ones(1,3)*0.9;

upneg=[410 470];
uwneg=[600 400];
uaneg=[0 0];

LB=[lppos lwpos lapos lpneg lwneg laneg];
UB=[uppos uwpos uapos upneg uwneg uaneg];
l0=(LB+UB)/2.*(1+(rand(1,numel(LB))-0.5)*0.5);
nelpos=numel(lapos);

OPTIONS =  optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',1500000,'TolFun',1e-6,'MaxIter',1500000);

[lres,errl,exitflag]=fmincon(@ba3620MA_fit,l0,[],[],[],[],LB,UB,[],OPTIONS);

l=lres;
ppos=l(1:3); l(1:3)=[];
wpos=l(1:2); l(1:2)=[];
apos=l(1:3)*max(ba3620MA(:,2)); l(1:3)=[];

pneg=l(1:2); l(1:2)=[];
wneg=l(1:2); l(1:2)=[];
aneg=l(1:2)*max(ba3620MA(:,2));

ws=[wpos(1) wpos(2) wpos(2)];

s1=zeros(numel(ba3x),5);

%----collecting positive gaussians
for k=1:3
    s1(:,k)=apos(k)*exp(-(ba3x-1e7/ppos(k)).^2/(2*ws(k)^2));
end

%----collecting negative gaussians
for k=4:5
    f=k-3;
    s1(:,k)=aneg(f)*exp(-(ba3x-1e7/pneg(f)).^2/(2*wneg(f)^2));
end

s_ba3MA620=s1;
sum_ba3MA620=sum(s1,2);

line([405 470],[0 0],'color','k')
plot(nm_ba3MA620,ba3620MA(1:length(nm_ba3MA620),2),'k','linewidth',lw)
plot(nm_ba3MA620,s_ba3MA620(:,1:3),'r','linewidth',lw,'linestyle','--')
plot(nm_ba3MA620,s_ba3MA620(:,4:5),'color',0.5*ones(1,3),'linewidth',lw,'linestyle','--')
plot(nm_ba3MA620,sum_ba3MA620,'r','linewidth',lw)
plot(nm_SMR,SMR,'k--','linewidth',lw)

sum_ba3MA620_pos=sum(s1(:,1:3),2);

save mat\ba3decomposition ba3620MA nm_ba3MA620 ba3x ba3y s1 sum_ba3MA620 sum_ba3MA620_pos

lb=0.0005;
line([405 470],[0 0],'color','k')
line(ppos(1)*ones(1,2),apos(1)+[-lb lb],'color','k')
line([422 ppos(1)],[0.005 apos(1)+lb],'color','k')
line(ppos(2)*ones(1,2),apos(2)+[-lb lb],'color','k')
line(ppos(3)*ones(1,2),apos(3)+[-lb lb],'color','k')
line([ppos(3) 450],[apos(3)+lb 0.0136],'color','k')
line([443 443],0.0131+[-lb lb],'color','k')
line([463 466],[-0.0027 -0.0018],'color','k')
line([429 427],[0.006 0.0075],'color','r')

[posnm_a3, poscm_a3, wcm_a3, wnm_a3, n1nm_a3, n2nm_a3, hm_a3] = fwhm_gauss(ppos',ws',apos','nm');
probe_a3=[posnm_a3 poscm_a3 wcm_a3 wnm_a3];
center_cm_a3=mean(poscm_a3(2:3));
center_nm_a3=mean(posnm_a3(2:3));

ngauss=1;
line([n1nm_a3(ngauss) n2nm_a3(ngauss)],hm_a3(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
line([425 417],[hm_a3(ngauss)+0.0002 hm_a3(ngauss)+0.0015],'color','k','linewidth',1)
text(417,hm_a3(ngauss)+0.002,[(num2str(wnm_a3(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','right')

ngauss=3;
line([n1nm_a3(ngauss) n2nm_a3(ngauss)],hm_a3(ngauss)*ones(1,2),'color','k','linestyle','--','linewidth',1)
line([444 451],[hm_a3(ngauss)+0.0003 hm_a3(ngauss)+0.003],'color','k','linewidth',1)
text(451,hm_a3(ngauss)+0.0035,[(num2str(wnm_a3(ngauss),'%0.0f')),' nm'],'fontsize',police,'HorizontalAlignment','left')

text(416,0.0052,num2str(ppos(1),'%0.0f'),'fontsize',police)
text(ppos(2)-0.5,apos(2)+lb+0.0006,num2str(ppos(2),'%0.0f'),'fontsize',police,'HorizontalAlignment','center')
text(450,0.014,num2str(ppos(3),'%0.0f'),'fontsize',police)
text(440,0.0141,'442.5','fontsize',police)
text(466,-0.0016,'\it1','fontsize',police)
text(425,0.008,'\it2','fontsize',police,'color','r')

line([436 439.7],[-0.0009 -0.0017],'color','k')
text(440,-0.002,'\it3','fontsize',police,'color','k')

xlim([405 470])
ylim([-0.004 0.015])

text(0.05,0.9,'623exc','fontsize',police,'Units','normalized')

x=[0 5 10 15]*1e-3; % define the x values where you want to have a tick
set(gca,'YTick',x);  % Apply the ticks to the current axes
set(gca,'YTickLabel', arrayfun(@(v) sprintf('%0.3f',v), x, 'UniformOutput', false) );

print('FigS2','-dpng')


%---------- save data in csv format ----------------------------
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2a\1_b2+a32+_dithionite.txt',ba3)
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2a\2_SMR_dithionite.txt',[SQR_c_nm SQR_c])

csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2b\1_b2+a32+_dithionite.txt',ba3)
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2b\2_SMR_dithionite.txt',[SQR_c_nm SQR_c])
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2b\3_pulse_570exc.txt',pulse570_aa3_ba3)
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2b\4_pulse_623exc.txt',pulse623_ba3)

csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2c\1_b2+a32+_570exc.txt',[nm_ba3MA560 ba3MA560])

csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2d\1_b2+a32+_623exc.txt',[nm_ba3MA620 ba3620MA(1:length(nm_ba3MA620),2)])
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2d\2_b2+a32+_623exc_fit.txt',[nm_ba3MA620 sum_ba3MA620])
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2d\3_b2+a32+_623exc_gaussians.txt',[nm_ba3MA620 s_ba3MA620(:,1:3)])
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2d\4_b2+a32+_623exc_negative_gaussians.txt',[nm_ba3MA620 s_ba3MA620(:,4:5)])
csvwrite('..\..\Data\Figures\FigS2_data\Panel_FigS2d\5_SMR_623exc.txt',[nm_SMR SMR])