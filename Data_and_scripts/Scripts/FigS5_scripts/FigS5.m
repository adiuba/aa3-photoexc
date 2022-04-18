close all
clear

police=16;
figure('position', [300 50 800 1000], 'paperpositionmode', 'auto');

ax{1}=axes('position', [0.17 0.57 0.8 0.41]); hold on
ax{2}=axes('position', [0.17 0.08 0.8 0.41]); hold on
ax{3}=axes('position', [0.74 0.15 0.2 0.13]); hold on
% letter='ABCDEFGH';
letter='abcdefgh';

ylabels={'Simulated absorption (r.u.)','Second derivative','Number of minima'};

x=linspace(1e7./500,1e7./400,10000)';
x1=443;
x2=433;
A0=[130 130];

Y=[]; D2Y=[]; N=[]; D3Y=[]; SN=[];

color=[1   0   0;
    1 0.2   0.2;
    1 0.4   0.4;
    1 0.6   0.6;
    0.7 0.7   1;
    0.6   0.6 1;
    0.4   0.4 1;
    0.2   0.2 1;
    0     0   1];
S=390:-10:310;
k=0;
for s1=S
    k=k+1;
    
    s=[s1 s1 600];
    
    A=A0*max(S)/s1;
    
    y0(:,1)=A(1)*exp(-(x-1e7./x1).^2/(2*s(1)^2));
    y0(:,2)=A(2)*exp(-(x-1e7./x2).^2/(2*s(2)^2));
    
    i11=find(y0(:,1) > max(y0(:,1))/2,1,'first');
    i12=find(y0(i11+1:end,1) < max(y0(:,1))/2,1,'first');
    i12=i11+i12;
    sn1=abs(1e7/x(i11)-1e7/x(i12));
    
    i21=find(y0(:,2) > max(y0(:,2))/2,1,'first');
    i22=find(y0(i11+1:end,2) < max(y0(:,2))/2,1,'first');
    i22=i21+i22;
    sn2=abs(1e7/x(i21)-1e7/x(i22));
    
    sn=mean([sn1 sn2]);
    SN=[SN; sn];
    
    y=sum(y0,2);
    
    Y=[Y y];
    dy=diff(y)./diff(x);
    d2y=diff(dy)./diff(x(2:end));
    d3y=diff(d2y);
    
    d3ya=d3y./abs(d3y);
    dd3ya=diff(d3ya);
    dd3yas=numel(find(dd3ya==2));
    n=dd3yas;
    N=[N; n];
    d2ys=sort(d2y);
%     d2y(1:2)
    
    D2Y=[D2Y d2y];
    D3Y=[D3Y d3y];
    axes(ax{1}), plot(1e7./x,y,'color',color(k,:),'linewidth',1.2)
    if s1==S(1) || s1==S(end)
        plot(1e7./x,y0,'color',color(k,:),'linewidth',1.2,'linestyle',':')
        if s1==S(1), Y1=y0; else YF=y0; end
    end
    axes(ax{2}), plot(1e7./x(3:end),d2y,'color',color(k,:),'linewidth',1.2)
    axes(ax{3}), plot(sn,n,'o','markerfacecolor',color(k,:),'markersize',5,'markeredgecolor',color(k,:))
        
end

for k=1:2
    axes(ax{k})
    xlim([410 480])
    set(gca,'fontsize',police,'linewidth',1)
    xl=get(gca,'xlim'); dxl=diff(xl);
    yl=get(gca,'ylim'); dyl=diff(yl);
    ylabel(ylabels{k},'fontsize',police+5,'position',[xl(1)-0.08*dxl mean(yl)])
    box on
    if k==2, xlabel('Wavelength (nm)','fontsize',police+5), end
    pos = get(gca,'Position');
    xpos=0.04; xposnorm=1-xpos/pos(3);
    ypos=0.02; yposnorm=1-ypos/pos(4);
    text(xposnorm,yposnorm,['\it',letter(k)],'fontsize',police,'Units','normalized')
end

axes(ax{3})
ylim([0.5 2.5])
xlim([6.5 9.5])
ylabel(ylabels{3},'fontsize',police)
xlabel('FWHH, nm','fontsize',police)
set(gca,'fontsize',police,'ytick',[1 2],'ticklength',[0.035 0.035])
box off

print('-r300', '-dpng', 'FigS5');

%---------- save data in csv format ----------------------------
csvwrite('..\..\Data\Figures\FigS5_data\Panel_FigS5a\1_narrowest_gaussians.txt',[1e7./x YF])
csvwrite('..\..\Data\Figures\FigS5_data\Panel_FigS5a\2_widest_gaussians.txt',[1e7./x Y1])
csvwrite('..\..\Data\Figures\FigS5_data\Panel_FigS5a\3_sum_of_two_gaussians.txt',[1e7./x Y])

csvwrite('..\..\Data\Figures\FigS5_data\Panel_FigS5b\1_second_derivatives.txt',[1e7./x(3:end) D2Y])

csvwrite('..\..\Data\Figures\FigS5_data\Panel_FigS5b_insert\1_nb_of_maxima_vs_FWHH.txt',[SN N])