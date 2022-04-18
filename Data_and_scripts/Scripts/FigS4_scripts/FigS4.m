close all
clear

lwY=5; lws=1; lwS=2;
police=18;
% letter='ABCDEFGH';
letter='abcdefgh';

ld{1}='Soret_formyl_polar_2g';
ld{2}='Soret_formyl_polar_3g';
ld{3}='Soret_formyl_nonpolar_2g';
ld{4}='Soret_formyl_nonpolar_3g';
ld{5}='Soret_schiff_polar_2g';
ld{6}='Soret_schiff_polar_3g';
ld{7}='Soret_schiff_nonpolar_2g';
ld{8}='Soret_schiff_nonpolar_3g';

txt1{1}='Heme \ita^{2+}'; txt2{1}='aqueous';
txt1{2}='Heme \ita^{2+}'; txt2{2}='aqueous';
txt1{5}='Schiff base'; txt2{5}='aqueous';
txt1{6}='Schiff base'; txt2{6}='aqueous';
txt1{3}='Heme \ita^{2+}'; txt2{3}='CH_2Cl_2';
txt1{4}='Heme \ita^{2+}'; txt2{4}='CH_2Cl_2';
txt1{7}='Schiff base'; txt2{7}='CH_2Cl_2';
txt1{8}='Schiff base'; txt2{8}='CH_2Cl_2';

xtxtpos=[0.1 0.58];
ixtxtpos=[1 1 1 1 2 2 2 2];

ytxtpos=[0.9 0.83];
iytxtpos=[1 1 1 1 2 2 2 2];

xl=[400 460]; yl=[0 1.1];
d=0.05*diff(yl);
dt=0.5*d;


figure('position', [5 50 1900 1000], 'paperpositionmode', 'auto');
ax{1}=axes('position', [0.05 0.58 0.21 0.4]);
ax{2}=axes('position', [0.29 0.58 0.21 0.4]);
ax{3}=axes('position', [0.53 0.58 0.21 0.4]);
ax{4}=axes('position', [0.77 0.58 0.21 0.4]);
ax{5}=axes('position', [0.05 0.08 0.21 0.4]);
ax{6}=axes('position', [0.29 0.08 0.21 0.4]);
ax{7}=axes('position', [0.53 0.08 0.21 0.4]);
ax{8}=axes('position', [0.77 0.08 0.21 0.4]);

for k=1:length(ax)
    axes(ax{k}), hold on, box on
    if exist(['mat\',ld{k},'.mat'],'file')
        load(['mat\',ld{k},'.mat'])
        plot(1e7./Xf,Yf,'color',0.5*ones(1,3),'linewidth',lwY)
        plot(1e7./Xf,sf,'--r','linewidth',lws)
        plot(1e7./Xf,Sf,'r','linewidth',lwS)
        
        MS=[]; MSI=[];
        for si=1:size(sf,2)
            [ms,msi]=max(sf(:,si));
            hsL{k}{si}=line(1e7./([1 1]*Xf(msi)),[-d d]+ms,'color','k');
            hsT{k}{si}=text(1e7./Xf(msi),ms+d+dt,num2str(round(1e7./Xf(msi))),'fontsize',police-5,...
                'HorizontalAlignment','center');
            MS=[MS ms]; MSI=[MSI msi];
        end
        [mS,mSi]=max(Sf);
        hSL{k}=line(1e7./([1 1]*Xf(mSi)),[-d d]+mS,'color','k');
        hST{k}=text(1e7./Xf(mSi),mS+d+dt,num2str(round(1e7./Xf(mSi))),'fontsize',police-5,...
            'HorizontalAlignment','center');
    end
    if k==1 || k==5, ylabel('Absorption (r.u.)','fontsize',police+5), end
    set(gca,'xlim',xl,'ylim',yl,'linewidth',1,'fontsize',police)
    text(xl(1)+diff(xl)*xtxtpos(ixtxtpos(k)),yl(1)+diff(yl)*ytxtpos(iytxtpos(k)),txt1{k},'fontsize',police+5)
    text(xl(1)+diff(xl)*xtxtpos(ixtxtpos(k)),yl(1)+diff(yl)*(ytxtpos(iytxtpos(k))-0.12),txt2{k},'fontsize',police+5)
    xlabel('Wavelength (nm)','fontsize',police+5)
    
    %additional annotations
    if k==1
        line(1e7./Xf(MSI(2)) + [0 0.08*diff(xl)],[d+MS(2) d+MS(2)+0.08*diff(yl)],'color','k');
        text(1e7./Xf(MSI(2))+0.09*diff(xl),d+MS(2)+0.09*diff(yl),num2str(round(1e7./Xf(MSI(2)))),'fontsize',police-5,...
                'HorizontalAlignment','left')
        line([421 453],max(Yf)/2*[1 1],'linestyle',':','color','k','linewidth',1.5);
        line([1 1]*440,[-d d]+mS,'color','k');
        text(440,max(Yf)/2+0.03*diff(yl),'35 nm','fontsize',police-5,...
                'HorizontalAlignment','center')
        text(440,mS+d+dt,'440','fontsize',police-5,...
                'HorizontalAlignment','center');
    elseif k==3
        line(1e7./Xf(MSI(2)) + [0 0.08*diff(xl)],[d+MS(2) d+MS(2)+0.08*diff(yl)],'color','k');
        text(1e7./Xf(MSI(2))+0.09*diff(xl),d+MS(2)+0.09*diff(yl),num2str(round(1e7./Xf(MSI(2)))),'fontsize',police-5,...
                'HorizontalAlignment','left')
        line(1e7./Xf(mSi) + [0 -0.08*diff(xl)],[d+mS d+mS+0.08*diff(yl)],'color','k');
        text(1e7./Xf(mSi)-0.09*diff(xl),d+mS+0.09*diff(yl),num2str(round(1e7./Xf(mSi))),'fontsize',police-5,...
                'HorizontalAlignment','right')
        line([424 450],max(Yf)/2*[1 1],'linestyle',':','color','k','linewidth',1.5);
        text(436,max(Yf)/2+0.03*diff(yl),'27 nm','fontsize',police-5,...
                'HorizontalAlignment','center')
    elseif k==5
        line([418 435],max(Yf)/2*[1 1],'linestyle',':','color','k','linewidth',1.5);
        text(426.5,max(Yf)/2+0.03*diff(yl),'19 nm','fontsize',police-5,...
                'HorizontalAlignment','center')
    elseif k==6
        line(1e7./Xf(MSI(3)) + [0 0.12*diff(xl)],[d+MS(3) d+MS(3)+0.12*diff(yl)],'color','k');
        text(1e7./Xf(MSI(3))+0.13*diff(xl),d+MS(3)+0.13*diff(yl),num2str(round(1e7./Xf(MSI(3)))),'fontsize',police-5,...
                'HorizontalAlignment','left')
    elseif k==7
        line(1e7./Xf(MSI(1)) + [0 0],MS(1)+0.03*diff(yl)+[-d d]*1.5,'color','k');
        text(1e7./Xf(MSI(1)),d*1.5+MS(1)+0.05*diff(yl),num2str(round(1e7./Xf(MSI(1)))),'fontsize',police-5,...
                'HorizontalAlignment','right')
        line([419 435.5],max(Yf)/2*[1 1],'linestyle',':','color','k','linewidth',1.5);
        text(428,max(Yf)/2+0.03*diff(yl),'18 nm','fontsize',police-5,...
                'HorizontalAlignment','center')
    elseif k==8
        line(1e7./Xf(MSI(3)) + [0 0.12*diff(xl)],[d+MS(3) d+MS(3)+0.12*diff(yl)],'color','k');
        text(1e7./Xf(MSI(3))+0.13*diff(xl),d+MS(3)+0.13*diff(yl),num2str(round(1e7./Xf(MSI(3)))),'fontsize',police-5,...
                'HorizontalAlignment','left')
    end
    pos = get(gca,'Position');
    xpos=0.02; xposnorm=1-xpos/pos(3);
    ypos=0.02; yposnorm=1-ypos/pos(4);
    text(xposnorm,yposnorm,['\it',letter(k)],'fontsize',police,'Units','normalized')
    
    panel=['..\..\Data\Figures\FigS4_data\Panel_FigS4',letter(k)];
    csvwrite([panel,'\1_digitized_spectrum.txt'],[1e7./Xf Yf])
    csvwrite([panel,'\2_fit.txt'],[1e7./Xf Sf])
    csvwrite([panel,'\3_gaussians.txt'],[1e7./Xf sf])
end

delete(hSL{1})
delete(hsL{7}{1})
delete(hsT{1}{2})
delete(hsT{3}{2})
delete(hsT{5}{2})
delete(hsT{6}{3})
delete(hsT{7}{1})
delete(hsT{8}{3})
delete(hST{1})
delete(hST{3})

print('-r300', '-dpng', 'FigS4');