clear
close all
clc

load mat\Fig_S1_data

police=16; lw=2; xl=[400 470];
figure('position', [100 50 1200 600], 'paperpositionmode', 'auto');
axpos{1}=[0.08 0.12 0.4 0.8];
axpos{2}=[0.56 0.121     0.4 0.8];

letter='abcdefgh';

% heme a
axes('position',axpos{1});
hold on
plot(aVanneste(:,1),aVanneste(:,2),'k','linewidth',lw)
plot(aOrii(:,1),aOrii(:,2),'r','linewidth',lw)
ylabel('Absorption (r.u.)','fontsize',police)
xlabel('Wavelength (nm)','fontsize',police)
box on
ylim([0 1.2])
xlim(xl)
set(gca,'fontsize',police)
pos = get(gca,'Position');
xpos=0.02; xposnorm=1-xpos/pos(3);
ypos=0.02; yposnorm=1-ypos/pos(4);
text(xposnorm,yposnorm,['\it',letter(1)],'fontsize',police,'Units','normalized','VerticalAlignment','middle')

line([444 444],[-0.05 0.05]+1,'color','k')
line([445 445],[-0.05 0.05]+1,'color','k')
line([441 444],[0.1 0.05]+1,'color','k')
line([448 445],[0.1 0.05]+1,'color','k')
text(440.5,1.12,'444','fontsize',police,'HorizontalAlignment','right')
text(448.5,1.12,'445','fontsize',police,'HorizontalAlignment','left')

line([425 425],[-0.05 0.05]+0.6297,'color','k')
text(425,0.705,'425','fontsize',police,'HorizontalAlignment','center')

text(420,1.15,'heme \ita^{2+}','fontsize',police,'HorizontalAlignment','center')
text(405,0.95,'\itVanneste (1966)','fontsize',police,'HorizontalAlignment','left')
text(405,0.85,'\itOrii (1998)','fontsize',police,'HorizontalAlignment','left','color','r')

% heme a3
axes('position',axpos{2});
hold on
plot(a3Vanneste(:,1),a3Vanneste(:,2),'k','linewidth',lw)
plot(a3Orii(:,1),a3Orii(:,2),'r','linewidth',lw)
xlabel('Wavelength (nm)','fontsize',police)
box on
ylim([0 1.2])
xlim(xl)
set(gca,'fontsize',police)

pos = get(gca,'Position');
xposnorm=1-xpos/pos(3);
yposnorm=1-ypos/pos(4);
text(xposnorm,yposnorm,['\it',letter(2)],'fontsize',police,'Units','normalized','VerticalAlignment','middle')

line([443 443],[-0.05 0.05]+1,'color','k')
text(443,1.07,'443','fontsize',police,'HorizontalAlignment','center')

line([426.1 426.1],[-0.05 0.05]+0.8766,'color','k')
text(426.1,0.95,'426','fontsize',police,'HorizontalAlignment','center')

line([415 415],[-0.05 0.05]+0.2548,'color','k')
text(415,0.33,'415','fontsize',police,'HorizontalAlignment','center')
text(420,1.135,'heme {\ita}_3^{2+}','fontsize',police,'HorizontalAlignment','center')

print('FigS1','-dpng')

%---------- save data in csv format ----------------------------
csvwrite('..\..\Data\Figures\FigS1_data\Panel_FigS1a\1_a2+_Vanneste.txt',aVanneste)
csvwrite('..\..\Data\Figures\FigS1_data\Panel_FigS1a\2_a2+_Orii.txt',aOrii)

csvwrite('..\..\Data\Figures\FigS1_data\Panel_FigS1b\1_a32+_Vanneste.txt',a3Vanneste)
csvwrite('..\..\Data\Figures\FigS1_data\Panel_FigS1b\2_a32+_Orii.txt',a3Orii)