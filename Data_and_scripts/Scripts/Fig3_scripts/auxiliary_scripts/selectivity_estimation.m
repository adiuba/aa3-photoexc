% clear
% clc
% close all
% break

global a3_Ca x

% load aa3
load ..\mat\dit_data
load ..\mat\pulses
ba3=dlmread('..\mat\ba3_metallomics_static_abs.txt'); ba3=ba3(ba3(:,1)>500&ba3(:,1)<680,:);
Ca=load('..\..\Fig2_scripts\mat\MV18all_res_v2_02042021');
yCa=Ca.yfull;

police=18;

aa3(:,1)=nmdit(nmdit>500 & nmdit<650);
aa3(:,2)=dit(nmdit>500 & nmdit<650);
aa3=aa3(end:-1:1,:); aa3=aa3(aa3(:,1)<=650,:);
aa3(:,2)=smooth(aa3(:,2),10);

x00 = [500 700];
y00 = [10.5 0];
c00 = [[1; 1]  x00(:)]\y00(:);
slope_m = c00(2);
intercept_b = c00(1);
aa3(:,2)=aa3(:,2)-(aa3(:,1)*slope_m + intercept_b);

nm_aa3=aa3(:,1); nm=nm_aa3; x=1e7./nm;

ba3=ba3(ba3(:,1)>500&ba3(:,1)<680,:);
x00 = [500 674];
y00 = [10 1.8];
c00 = [[1; 1]  x00(:)]\y00(:);
slope_m = c00(2);
intercept_b = c00(1);
ba3(:,2)=ba3(:,2)-(ba3(:,1)*slope_m + intercept_b);

% %-------heme b----------
% y2=-imag(5.714e6*A2*(x2^2-x_2.^2-1i*x_2*s2)./((x2.^2-x_2.^2).^2+x_2.^2.*s2^2)); %Lorentz
bamp=13;
% y2_2=bamp*exp(-(1e7/558-x_2).^2/(2*200^2)); % Gauss in nm2
b=bamp*exp(-(1e7/558-x).^2/(2*200^2)); % Gauss in nm

pulse623_aa3_new=new_grid(pulse623_aa3(:,1),nm,pulse623_aa3(:,2));
pulse623_aa3_new(2:end)=pulse623_aa3_new(2:end)-diff(pulse623_aa3_new(1:2));

pulse623_ba3_new=new_grid(pulse623_ba3(:,1),nm,pulse623_ba3(:,2));
pulse623_ba3_new(2:end)=pulse623_ba3_new(2:end)-diff(pulse623_ba3_new(1:2));

pulse608_aa3_new=new_grid(pulse608_aa3(:,1),nm,pulse608_aa3(:,2));
pulse608_aa3_new(2:end)=pulse608_aa3_new(2:end)-diff(pulse608_aa3_new(1:2));

pulse600_ba3_new=new_grid(pulse600_ba3(:,1),nm,pulse600_ba3(:,2));
pulse600_ba3_new(2:end)=pulse600_ba3_new(2:end)-diff(pulse600_ba3_new(1:2));

pulse570_aa3_ba3_new=new_grid(pulse570_aa3_ba3(:,1),nm,pulse570_aa3_ba3(:,2));
pulse570_aa3_ba3_new(2:end)=pulse570_aa3_ba3_new(2:end)-diff(pulse570_aa3_ba3_new(1:2));

yCa=new_grid(Ca.nmMV18V(end-1:-1:1),nm,yCa(end:-1:1));
% yCa(2:end)=yCa(2:end)-diff(yCa(1:2));
yCa=yCa-yCa(end);

max_aa3_alpha=max(aa3(:,2));
k_a_Ca=0.78;
a_Ca_abs_max=k_a_Ca*max_aa3_alpha;

a_Ca=yCa/max(yCa)*a_Ca_abs_max;
a3_Ca=aa3(:,2)-a_Ca;
a3_Ca(nm>635)=0;
k_a3_Ca=max(a3_Ca)/max_aa3_alpha;

aa3_ba3=[zeros(17,1); a3_Ca(1:end-17)];
% kb=max(b)/(max(b)+max(a3_Ca));
kb=0.6;

conv17=-trapz(x,pulse623_aa3_new.*a_Ca);
conv18=-trapz(x,pulse623_aa3_new.*a3_Ca);

sel_623_aa3_Ca2=conv17/conv18

conv19=-trapz(x,pulse608_aa3_new.*a_Ca);
conv20=-trapz(x,pulse608_aa3_new.*a3_Ca);

sel_608_aa3_Ca2=conv19/conv20

conv21=-trapz(x,pulse570_aa3_ba3_new.*a_Ca);
conv22=-trapz(x,pulse570_aa3_ba3_new.*a3_Ca);

sel_570_aa3_Ca2=conv21/conv22

conv9_2=-trapz(x,pulse570_aa3_ba3_new.*aa3_ba3);
conv10_2=-trapz(x,pulse570_aa3_ba3_new.*b);

sel_570_ba3_Ca_2=conv9_2/conv10_2

% figure
% plot(aa3(:,1),aa3(:,2))
% hold on
% plot(nm,a_Ca)
% plot(nm,a3_Ca)
% % plot(nm,na3_Ca,'r')
% xlim([500 700])

% %% plotting according to the alternative scheme
% 
% figure('position', [910 50 900 900], 'paperpositionmode', 'auto');
% axpos{1}=[0.12 0.58 0.85 0.4];
% axpos{2}=[0.12 0.1 0.85 0.4];
% axpos{3}=[0.78 0.72 0.15 0.15];
% axpos{4}=[0.78 0.24 0.15 0.15];
% % axpos{4}=[0.55 0.15 0.4 0.38];
% 
% lw=2; ax=cell(1,length(axpos));
% for k=1:length(axpos)
%     ax{k}=axes('position',axpos{k});
%     hold on
%     set(gca,'fontsize',police,'linewidth',1)
%     box on
% end
% 
% 
% 
% axes(ax{1})
% fill(nm,b/max(b)*kb,'b','EdgeColor','none')
% fill([nm; nm(end:-1:1)],[zeros(length(aa3_ba3),1); aa3_ba3(end:-1:1)]/max(aa3_ba3)*k_a3_Ca,'g','EdgeColor','none')
% % fill(nm2,y1/max(y1)*ky1,'g','EdgeColor','none')
% alpha(0.5)
% % plot(heme_a3_Palmer(:,1),heme_a3_Palmer(:,2)/max(heme_a3_Palmer(:,2))*ky1,'color',[0 0.6 0],'linewidth',lw)
% % nCa_2_ba3/max(nCa_2_ba3)*k_a3_Ca
% plot(ba3(:,1),ba3(:,2)/max(ba3(:,2))*(kb+aa3_ba3(nm>557.9&nm<558.2)/max(aa3_ba3)*k_a3_Ca),'color',0.6*ones(1,3),'linewidth',lw)
% plot(nm,b/max(b)*kb+aa3_ba3/max(aa3_ba3)*k_a3_Ca,'color',[252 119 3]/255,'linewidth',lw)
% plot(pulse623_ba3(:,1),pulse623_ba3(:,2),'r','linewidth',lw)
% plot(pulse570_aa3_ba3(:,1),pulse570_aa3_ba3(:,2),'b','linewidth',lw)
% xlim([500 700])
% ylim([0 1])
% set(gca,'xtick',500:20:700)
% 
% text(560,0.2,'heme','fontsize',police,'HorizontalAlignment','center')
% text(560,0.1,'\itb','fontsize',police,'HorizontalAlignment','center')
% 
% text(622,0.07,'heme {\ita}_3','fontsize',police,'HorizontalAlignment','center')
% 
% text(540,0.92,'570 exc','fontsize',police,'HorizontalAlignment','center','color','b')
% text(540,0.82,['{\ita}_3/{\itb} \approx ', num2str(sel_570_ba3_Ca_2,'%.1f')],...
%     'fontsize',police,'HorizontalAlignment','center','color','b')
% 
% text(650,0.92,'623 exc','fontsize',police,'HorizontalAlignment','center','color','r')
% % text(650,0.82,['{\ita}_3/{\itb} \approx ', num2str(sel_623_ba3_fit_Ca,'%.1f')],...
% %     'fontsize',police,'HorizontalAlignment','center','color','r')
% text(650,0.82,['{\ita}_3/{\itb} \rightarrow \infty'],...
%     'fontsize',police,'HorizontalAlignment','center','color','r')
% 
% % xlabel('Wavelength (nm)','fontsize',police)
% ylabel('Relative units','fontsize',police)
% 
% axes(ax{2})
% % fill(nm2,y4/max(y4)*ky4,'r','EdgeColor','none')
% fill([nm; nm(end:-1:1)],[zeros(length(a_Ca),1); a_Ca(end:-1:1)]/max(a_Ca)*k_a_Ca,'r','EdgeColor','none')
% fill([nm; nm(end:-1:1)],[zeros(length(a3_Ca),1); a3_Ca(end:-1:1)]/max(a3_Ca)*k_a3_Ca,'g','EdgeColor','none')
% % fill(nm2,y3/max(y3)*ky3,'g','EdgeColor','none')
% % 'color',[0 0.6 0]
% 
% alpha(0.5)
% plot(nm,aa3(:,2)/max(aa3(:,2)),'color',0.6*ones(1,3),'linewidth',lw)
% plot(nm,a_Ca/max(a_Ca)*k_a_Ca+a3_Ca/max(a3_Ca)*k_a3_Ca,'color',[252 119 3]/255,'linewidth',lw)
% plot(pulse623_aa3(:,1),pulse623_aa3(:,2),'r','linewidth',lw)
% plot(pulse608_aa3(:,1),pulse608_aa3(:,2),'color',[0 0.7 0],'linewidth',lw)
% plot(pulse570_aa3_ba3(:,1),pulse570_aa3_ba3(:,2),'b','linewidth',lw)
% plot(nmdit,dit/max(dit(nmdit>590)),'color','m','linewidth',lw)
% xlim([500 700])
% ylim([0 1])
% set(gca,'xtick',500:20:700)
% 
% text(540,0.92,'570 exc','fontsize',police,'HorizontalAlignment','center','color','b')
% text(540,0.82,['{\ita}/{\ita}_3 \approx ', num2str(sel_570_aa3_Ca2,'%.2f')],...
%     'fontsize',police,'HorizontalAlignment','center','color','b')
% 
% text(590,0.92,'608 exc','fontsize',police,'HorizontalAlignment','center','color',[0 0.7 0])
% text(590,0.82,'{\ita}/{\ita}_3',...
%     'fontsize',police,'HorizontalAlignment','center','color',[0 0.7 0])
% text(587,0.72,['\approx ', num2str(sel_608_aa3_Ca2,'%.2f')],...
%     'fontsize',police,'HorizontalAlignment','center','color',[0 0.7 0])
% 
% text(650,0.92,'623 exc','fontsize',police,'HorizontalAlignment','center','color','r')
% text(650,0.82,['{\ita}/{\ita}_3 \approx ', num2str(sel_623_aa3_Ca2,'%.2f')],...
%     'fontsize',police,'HorizontalAlignment','center','color','r')
% 
% text(593,0.13,'heme {\ita}_3','fontsize',police,'HorizontalAlignment','center')
% text(530,0.14,'heme {\ita}','fontsize',police,'HorizontalAlignment','center')
% % text(600,0.92,'608 exc','fontsize',police,'HorizontalAlignment','center','color',[0 0.7 0])
% text(520,0.32,'{\itaa}_3','fontsize',police,'HorizontalAlignment','center','color',0.6*ones(1,3))
% 
% 
% 
% xlabel('Wavelength (nm)','fontsize',police)
% ylabel('Relative units','fontsize',police)

% axes(ax{3})
% plot(1e7./dx-0.5,sel_623_ba3_pulseposition_Ca,'linewidth',lw)
% 
% xlabel('Wavelength (nm)','fontsize',police-5)
% ylabel('Selectivity {\ita}_3/{\itb}','fontsize',police-5)
% set(gca,'fontsize',police-5)
% box off
% line([623 623],[0.1 35],'linestyle','--','color','r','linewidth',lw)
% text(620,37,'623 exc','fontsize',police-5,'HorizontalAlignment','right','color','r')
% xlim([550 650])
% ylim([0 40])

% axes(ax{4})
% plot(1e7./dx2-3,sel_608_aa3_pulseposition,'linewidth',lw)
% 
% xlabel('Wavelength (nm)','fontsize',police-5)
% ylabel('Selectivity {\ita}/{\ita}_3','fontsize',police-5)
% set(gca,'fontsize',police-5)
% box off
% line([608 608],[0.1 5],'linestyle','--','color','r','linewidth',lw)
% text(605,4,'608 exc','fontsize',police-5,'HorizontalAlignment','right','color','r')

% save selectivity2
% print('S2_2_new_scheme2','-dpng')

save ..\mat\selectivity

    
