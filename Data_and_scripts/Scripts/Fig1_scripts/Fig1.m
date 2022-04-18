clc
clear
close all

load mat\dit_data;

set(0,'defaulttextinterpreter','tex')

figure('position', [10 50 900 900], 'paperpositionmode', 'auto');
axpos{1}=[0.10 0.55 0.87 0.44];
axpos{2}=[0.10 0.07 0.38 0.4];
axpos{3}=[0.57 -0.01 0.45 0.45];

police=16; lw=1.5;

ca=[204 51 0]/255;
ca3=[0 102 0]/255;

ax=cell(length(axpos),1);
letter='abcdefgh';

for k=1:length(axpos)
    ax{k}=axes('position',axpos{k});
    hold on
    set(gca,'fontsize',police)
    set(gca,'linewidth',lw)
    set(gca,'ticklength',[0.015 0.015])
    box on
    pos = get(gca,'Position');
    if k==1 || k==2
        xpos=0.02; xposnorm=1-xpos/pos(3);
        ypos=0.01; yposnorm=1-ypos/pos(4);
        text(xposnorm,yposnorm,['\it',letter(k)],'fontsize',police,'Units','normalized')
    elseif k==3
        xpos=0.05; xposnorm=1-xpos/pos(3);
        ypos=-0.02; yposnorm=1-ypos/pos(4);
        text(xposnorm,yposnorm,['\it',letter(k)],'fontsize',police,'Units','normalized')
    end
    set(gca,'ytick',0:100:200)
    if k==2 || k==3, axis off, end
end

axes(ax{1})
hold on
% x00 = [500 700];
% y00 = [10.5 0];
% c00 = [[1; 1]  x00(:)]\y00(:);
% slope_m = c00(2);
% intercept_b = c00(1);
% aa3=dit-(nmdit*slope_m + intercept_b);


plot(nmdit,dit,'k','linewidth',lw)
xlim([390 700])
ylim([0 320])
xlabel('Wavelength (nm)','fontsize',police)
ylabel('A/[{\itaa}_3] (mM^{-1}cm^{-1})','fontsize',police)



yl=get(gca,'ylim');
df=0.05*diff(yl);
line([420 420],120+[-df df],'linewidth',1,'color','k')
line([443.5 443.5],253.6+[-df df],'linewidth',1,'color','k')
line([604.5 604.5],50.75+[-df df],'linewidth',1,'color','k')
line([400 460],286*ones(2,1),'linewidth',lw,'color','k')
line([500 570],45*ones(2,1),'linewidth',lw,'color','k')
line([590 630],100*ones(2,1),'linewidth',lw,'color','k')

text(630,250,'{\ita}^{2+}{\ita}_3^{2+}','fontsize',police,'HorizontalAlignment','left')

text(422,120+df+10,'420','fontsize',police,'HorizontalAlignment','right')
text(445,253+df+3,'444','fontsize',police,'HorizontalAlignment','left')
text(605,51+df+10,'605','fontsize',police,'HorizontalAlignment','left')

text(430,300,'\gamma = B = Soret','fontsize',police,'HorizontalAlignment','center')
text(535,59,'\beta = Q_v','fontsize',police,'HorizontalAlignment','center')
text(610,114,'\alpha = Q_{0}','fontsize',police,'HorizontalAlignment','center')


axes(ax{2})
hold on
kxshift=-0.02; kyshift=0.08; lw2=1;
annotation('arrow',0.12*ones(1,2)+kxshift,[0 0.33]+kyshift,'linewidth',lw2)
annotation('arrow',0.14*ones(1,2)+kxshift,[0 0.3]+kyshift,'linewidth',lw2)
annotation('arrow',0.22*ones(1,2)+kxshift,[0 0.32]+kyshift,'linewidth',lw2)
annotation('arrow',0.24*ones(1,2)+kxshift,[0 0.29]+kyshift,'linewidth',lw2)

annotation('arrow',0.37*ones(1,2)+kxshift,[0 0.18]+kyshift,'linewidth',lw2)
annotation('arrow',0.39*ones(1,2)+kxshift,[0 0.15]+kyshift,'linewidth',lw2)
annotation('arrow',0.47*ones(1,2)+kxshift,[0 0.17]+kyshift,'linewidth',lw2)
annotation('arrow',0.49*ones(1,2)+kxshift,[0 0.14]+kyshift,'linewidth',lw2)

annotation('line',[-0.00 0.41]+0.1+kxshift,0*ones(1,2)+kyshift-0.005,'linewidth',lw2)

annotation('line',[-0.02 0.06]+0.12+kxshift,0*ones(1,2)+kyshift+0.305,'linewidth',lw2)
annotation('line',[-0.02 0.04]+0.12+kxshift,0*ones(1,2)+kyshift+0.335,'linewidth',lw2)

annotation('line',[-0.02 0.06]+0.22+kxshift,0*ones(1,2)+kyshift+0.295,'linewidth',lw2)
annotation('line',[-0.02 0.04]+0.22+kxshift,0*ones(1,2)+kyshift+0.325,'linewidth',lw2)

annotation('line',[-0.02 0.06]+0.37+kxshift,0*ones(1,2)+kyshift+0.155,'linewidth',lw2)
annotation('line',[-0.02 0.04]+0.37+kxshift,0*ones(1,2)+kyshift+0.185,'linewidth',lw2)

annotation('line',[-0.02 0.06]+0.47+kxshift,0*ones(1,2)+kyshift+0.145,'linewidth',lw2)
annotation('line',[-0.02 0.04]+0.47+kxshift,0*ones(1,2)+kyshift+0.175,'linewidth',lw2)

annotation('textbox',[0.04+kxshift kyshift-0.035 0.05 0.05],'String','S_0','fontsize',police,'EdgeColor','none')

annotation('textbox',[0.105+kxshift kyshift+0.33 0.05 0.05],'String','B_y','fontsize',police,'EdgeColor','none')
annotation('textbox',[0.064+kxshift kyshift+0.12 0.05 0.05],'String',{'B_{1y}','(\gamma)'},'fontsize',police,'EdgeColor','none')
annotation('textbox',[0.136+kxshift kyshift+0.15 0.05 0.05],'String',{'B_{0y}','(\gamma)'},'fontsize',police,'EdgeColor','none')

annotation('textbox',[0.205+kxshift kyshift+0.32 0.05 0.05],'String','B_x','fontsize',police,'EdgeColor','none')
annotation('textbox',[0.164+kxshift kyshift+0.05 0.05 0.05],'String',{'B_{1x}','(\gamma)'},'fontsize',police,'EdgeColor','none')
annotation('textbox',[0.236+kxshift kyshift+0.13 0.05 0.05],'String',{'B_{0x}','(\gamma)'},'fontsize',police,'EdgeColor','none')


annotation('textbox',[0.355+kxshift kyshift+0.18 0.05 0.05],'String','Q_y','fontsize',police,'EdgeColor','none')
annotation('textbox',[0.313+kxshift kyshift+0.05 0.05 0.05],'String',{'Q_{1y}','(\beta)'},'fontsize',police,'EdgeColor','none')
annotation('textbox',[0.386+kxshift kyshift+0.09 0.05 0.05],'String',{'Q_{0y}','(\alpha)'},'fontsize',police,'EdgeColor','none')

annotation('textbox',[0.455+kxshift kyshift+0.17 0.05 0.05],'String','Q_x','fontsize',police,'EdgeColor','none')
annotation('textbox',[0.413+kxshift kyshift+0.02 0.05 0.05],'String',{'Q_{1x}','(\beta)'},'fontsize',police,'EdgeColor','none')
annotation('textbox',[0.486+kxshift kyshift+0.05 0.05 0.05],'String',{'Q_{0x}','(\alpha)'},'fontsize',police,'EdgeColor','none')


axes(ax{3})
hold on
hemeaxy = imread('mat\Heme_a_2.tif');
imshow(hemeaxy)
text(0.58,0.21,'x','fontsize',police,'Units','normalized')
text(0.58,0.76,'y','fontsize',police,'Units','normalized')

print('-dpng','-r400','Fig1')

%---------- save data in csv format
csvwrite('..\..\Data\Figures\Fig1_data\Panel_Fig1a\1_aa3_dithionite.txt',[nmdit dit])