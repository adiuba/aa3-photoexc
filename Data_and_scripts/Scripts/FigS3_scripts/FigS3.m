clear
close all
clc

load mat\Orii

koef_c=2;

Orii_x{1}=x1_Orii;
Orii_x{2}=x2_Orii;

Orii{1}=Orii_var1;
Orii{2}=Orii_var2;
Orii{3}=Orii2_var1*koef_c;
Orii{4}=Orii2_var2*koef_c;

Orii_hemes{1}=Orii_var1_hemes;
Orii_hemes{2}=Orii_var2_hemes;
Orii_hemes{3}=Orii2_var1_hemes*koef_c;
Orii_hemes{4}=Orii2_var2_hemes*koef_c;

set(0,'defaulttextinterpreter','tex')

figure('position', [10 50 800 1200], 'paperpositionmode', 'auto');
axpos{1}=[0.12 0.77 0.38 0.22];
axpos{2}=[0.58 0.77 0.38 0.22];
axpos{3}=[0.12 0.54 0.38 0.22];
axpos{4}=[0.58 0.54 0.38 0.22];
axpos{5}=[0.12 0.30 0.38 0.22];
axpos{6}=[0.58 0.30 0.38 0.22];
axpos{7}=[0.12 0.07 0.38 0.22];
axpos{8}=[0.58 0.07 0.38 0.22];

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
    xpos=0.02; xposnorm=1-xpos/pos(3);
    ypos=0.01; yposnorm=1-ypos/pos(4);
    text(xposnorm,yposnorm,['\it',letter(k)],'fontsize',police,'Units','normalized')
    xlim([380 480])
    ylim([0 250])
    set(gca,'ytick',0:100:200)
    if k==1 || k==2 || k==3 || k==4 || k==5 || k==6
        set(gca,'xticklabel',[])
    end
    
    if k==1
        
        line([385 398],[225 225],'color','k','linewidth',lw)
        line([385 398],[180 180],'color','r','linewidth',lw)
        line([430 443],[225 225],'color',0.5*ones(1,3),'linewidth',lw,'linestyle','--')
        line([430 443],[180 180],'color',ca,'linewidth',lw,'linestyle','--')
        
        text(402,225,'{\ita}^{2+}{\ita}_3^{3+}','HorizontalAlignment','left','fontsize',police)
        text(402,180,'sum','HorizontalAlignment','left','fontsize',police)
        text(447,225,'{\ita}_3^{3+}','HorizontalAlignment','left','fontsize',police)
        text(447,180,'{\ita}^{2+}','HorizontalAlignment','left','fontsize',police)
        
    elseif k==3
        
        line([385 398],[220 220],'color',0.5*ones(1,3),'linewidth',lw)
        line([385 398],[175 175],'color',ca,'linewidth',lw)
        line([385 398],[130 130],'color',ca3,'linewidth',lw)
        
        text(402,220,'{\ita}^{2+}{\ita}_3^{2+}','HorizontalAlignment','left','fontsize',police,'color',0.5*ones(1,3))
        text(402,175,'{\ita}^{2+}','HorizontalAlignment','left','fontsize',police,'color',ca)
        text(402,130,'{\ita}_3^{2+}','HorizontalAlignment','left','fontsize',police,'color',ca3)
    end
    
    if k==1 || k==5
        text(0.8,0.6,'Var. 1','fontsize',police,'Units','normalized')
    elseif k==2 || k==6
        text(0.8,0.6,'Var. 2','fontsize',police,'Units','normalized')
    end
end

xlabel('Wavelength (nm)','position',[365 -40])
ylabel('A/[{\itaa}_3] (mM^{-1}cm^{-1})','position',[242 560])

l=[1 2 5 6]; xl=[2 2 1 1]; varl=[3 4 1 2];
for k=1:4
    axes(ax{l(k)})
    hold on
    
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,1),'k','linewidth',lw)
    
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,2),'--','color',0.5*ones(1,3),'linewidth',lw)
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,3),'--','color',0.5*ones(1,3),'linewidth',lw)
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,4),'--','color',ca,'linewidth',lw)
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,5),'--','color',ca,'linewidth',lw)
    
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,6),'r','linewidth',lw)
    
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,8),'--','color',ca,'linewidth',lw)
    plot(Orii_x{xl(k)},Orii{varl(k)}(:,9),'--','color',0.5*ones(1,3),'linewidth',lw)
    
    panel=['..\..\Data\Figures\FigS3_data\Panel_FigS3',letter(l(k))];
    
    csvwrite([panel,'\1_mixed_valence.txt'],[Orii_x{xl(k)} Orii{varl(k)}(:,1)])
    csvwrite([panel,'\2_mixed_valence_fit.txt'],[Orii_x{xl(k)} Orii{varl(k)}(:,6)])
    csvwrite([panel,'\3_a2+_gaussians.txt'],[Orii_x{xl(k)} Orii{varl(k)}(:,4) Orii{varl(k)}(:,5)])
    csvwrite([panel,'\4_a2+_baseline.txt'],[Orii_x{xl(k)} Orii{varl(k)}(:,8)])
    csvwrite([panel,'\5_a33+_gaussians.txt'],[Orii_x{xl(k)} Orii{varl(k)}(:,2) Orii{varl(k)}(:,3)])
    csvwrite([panel,'\6_a33+_baseline.txt'],[Orii_x{xl(k)} Orii{varl(k)}(:,9)])
    
    axes(ax{l(k)+2})
    hold on
    
    plot(Orii_x{xl(k)},Orii_hemes{varl(k)}(:,1),'color',ca,'linewidth',lw)
    plot(Orii_x{xl(k)},Orii_hemes{varl(k)}(:,2),'color',ca3,'linewidth',lw)
    plot(Orii_x{xl(k)},Orii_hemes{varl(k)}(:,4),'color',0.5*ones(1,3),'linewidth',lw)
    
    panel=['..\..\Data\Figures\FigS3_data\Panel_FigS3',letter(l(k)+2)];
    csvwrite([panel,'\1_a2+a32+.txt'],[Orii_x{xl(k)} Orii_hemes{varl(k)}(:,4)])
    csvwrite([panel,'\2_a2+.txt'],[Orii_x{xl(k)} Orii_hemes{varl(k)}(:,1)])
    csvwrite([panel,'\3_a32+.txt'],[Orii_x{xl(k)} Orii_hemes{varl(k)}(:,2)])
end

annotation('arrow',0.32*ones(1,2),[0.05 0]+0.715)
annotation('arrow',0.32*ones(1,2),[0.05 0]+0.245)
annotation('arrow',0.78*ones(1,2),[0.05 0]+0.715)
annotation('arrow',0.78*ones(1,2),[0.05 0]+0.245)

print('-dpng','-r400','FigS3')
