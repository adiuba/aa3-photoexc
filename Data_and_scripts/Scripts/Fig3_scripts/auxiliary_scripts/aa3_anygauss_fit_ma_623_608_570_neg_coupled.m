function err = aa3_anygauss_fit_ma_623_608_570_neg_coupled(l)

global X drw Npos Nneg sp err_all W dit nmdit

a_to_a3_predicted=[0.77 3.15 4.95];

olive=[51 153 51]/255;

axpos{1}=[0.03 0.55 0.2 0.3];
axpos{2}=[0.28 0.55 0.2 0.3];
axpos{3}=[0.53 0.55 0.2 0.3];
axpos{4}=[0.03 0.11 0.2 0.3];
axpos{5}=[0.28 0.11 0.2 0.3];
axpos{6}=[0.53 0.11 0.2 0.3];
axpos{7}=[0.785 0.3 0.2 0.5];

titles{1}='570 exc';
titles{2}='608 exc';
titles{3}='623 exc';

L=l;

ppos=l(1:Npos); l(1:Npos)=[];
wpos=l(1:numel(W)); l(1:numel(W))=[];
apos=l(1:Npos*3); l(1:Npos*3)=[];
apos=reshape(apos,Npos,3);

pneg=l(1:sum(Nneg)); l(1:sum(Nneg))=[];
wneg=l(1:sum(Nneg)); l(1:sum(Nneg))=[];
aneg=l(1:sum(Nneg)*3); %l(1:sum(Nneg))=[];
aneg=reshape(aneg,sum(Nneg),3);

ws=[];
for k=1:numel(W)
    ws=[ws wpos(k)*ones(1,W(k))];
end

% nc=0;
SP=[]; smodel=[];
for m=1:3
    
    s1=zeros(numel(X),Npos+Nneg);
    %----collecting positive gaussians
    for k=1:Npos
        s1(:,k)=apos(k,m)*exp(-(X-1e7/ppos(k)).^2/(2*ws(k)^2));
    end
    
    %----collecting negative gaussians
    for k=Npos+1:Npos+Nneg
%         f=nc+k-Npos;
        f=k-Npos;
        s1(:,k)=aneg(f,m)*exp(-(X-1e7/pneg(f)).^2/(2*wneg(f)^2));
    end
%     nc=nc+Nneg;
    
    S1(:,m)=sum(s1,2);
    err1(m)=sqrt(sum((sp(:,m)-S1(:,m)).^2))/max(sp(:,m));
    %     if deg
    %         s1(:,3)=sum([s1(:,3) s1(:,4)],2);
    %         s1(:,4)=[];
    %     end
    S{m}=s1;
    SP=[SP; sp(:,m)];
    smodel=[smodel; S1(:,m)];
    %     size(s1)
    % %     m
    %     plot(1e7./X,[sp(:,m) S1(:,m) s1])
    %     pause
end

err=sum(err1);
err_all=[err_all err];
chi=sum((SP-smodel).^2./abs(smodel));
degrees_of_freedom=numel(SP)-(numel(L)-1);

% pause

% s=sum(s1,2);
%
% sq1=trapz(X,sum(s1(:,4:5),2))/trapz(X,sum(s1(:,[3 6 7]),2))
% sq2=trapz(X,sum(s1(:,4:5),2))/trapz(X,sum(s1(:,[6 7]),2))
% sq2=trapz(X,sum(s1(:,3:4),2))/trapz(X,sum(s1(:,[1 5]),2))
% sq3=trapz(X,sum(s1(:,[1 3 4]),2))/trapz(X,s1(:,5))


if drw
    drawnow
    figure('position', [10 50 1600 800], 'paperpositionmode', 'auto');
    a3=sum(S{2}(:,[1 3 4]),2); a=sum(S{2}(:,[2 5 6]),2);
    a_a3_608=trapz(X,a)/trapz(X,a3), selrel=a_a3_608/a_to_a3_predicted(2)
    
    for m=1:3
        %         subplot('position',[Mx(m) My(m) 5 5]/19)
        axes('position',axpos{m})
        
        plot(1e7./X, S{m}(:,[1 3 4]),'m','linewidth',2)
        hold on
        plot(1e7./X, S{m}(:,[2 5 6]),'k','linewidth',2)
        plot(1e7./X, S{m}(:,7:end),'color',[0 0.7 0],'linewidth',2)
        plot(1e7./X, sp(:,m),'color',[0.8 0.8 0.8],'linewidth',6)
        plot(1e7./X, S1(:,m),'color','r','linewidth',3)
        %         xlabel('\lambda, nm','Fontsize',26)
        %         ylabel('absorption, r.u.', 'Fontsize',26)
        set(gca,'fontsize',16,'Xlim',[400 480],'Ylim',[-0.2 1.1])%[min(sp(:,m+1))-0.1 1])
        box on
        text(470,0.4,'{\ita}','fontsize',16)
        text(470,0.2,'{\ita}_3','color','m','fontsize',16)
        xlabel('Wavelength (nm)','fontsize',16)
        xlabel('Wavelength (nm)','fontsize',16)
        if m==1
            text(500,1.4,'\bfFitting of the inverted fastest EASs following photoexcitation','fontsize',16)
        end
        title(titles{m},'fontsize',16)
        
%         subplot(2,3,m+3)
        axes('position',axpos{m+3})
        hold on
        a3=sum(S{m}(:,[1 3 4]),2);
        a=sum(S{m}(:,[2 5 6]),2);
%         plot(1e7./X, S{m}(:,[1 3 4]),'m','linewidth',2)
        plot(1e7./X, sum(S{m}(:,[1 3 4]),2),'m','linewidth',2)
        set(gca,'fontsize',16,'Xlim',[400 480],'Ylim',[-0.2 1.1])%[min(sp(:,m+1))-0.1 1])
        box on
%         a_a3=trapz(a)/trapz(a3);
%         text(420,0.6,['{\ita}/{\ita}_3=',num2str(a_a3,'%.1f')],'fontsize',16)
%         text(420,0.5,['real {\ita}/{\ita}_3=',num2str(a_a3/0.88,'%.1f')],'fontsize',16)
        
%         subplot(2,3,m+3)
%         axes('position',axpos{m+3})
%         hold on
%         plot(1e7./X, S{m}(:,[2 5 6]),'m','linewidth',2)
        plot(1e7./X, sum(S{m}(:,[2 5 6]),2),'k','linewidth',2)
        a_a3=trapz(X,a)/trapz(X,a3)
        text(405,0.9,['{\ita}/{\ita}_3=',num2str(a_a3,'%.1f')],'fontsize',16)
%         text(405,0.7,['real {\ita}/{\ita}_3=',num2str(a_a3/a_to_a3,'%.1f')],'fontsize',16)
        text(405,0.7,['real {\ita}/{\ita}_3=',num2str(a_a3/selrel,'%.1f')],'fontsize',16)
        text(470,0.4,'{\ita}','fontsize',16)
        text(470,0.2,'{\ita}_3','color','m','fontsize',16)
        xlabel('Wavelength (nm)','fontsize',16)
    end
    axes('position',axpos{7})
    hold on
    kscat=5e+11; scat=kscat*nmdit.^(-4); dit=dit-scat;
    ditSoretNorm=dit(nmdit>400&nmdit<480); nmditSoret=nmdit(nmdit>400&nmdit<480);
    ditSoretNorm=ditSoretNorm-ditSoretNorm(nmditSoret>470&nmditSoret<470.5);
    ditSoretNorm=ditSoretNorm/max(ditSoretNorm);
%     k=trapz(a/max(a))/trapz(a3/max(a3)); k1=a_to_a3/k;
    k=trapz(a/max(a))/trapz(a3/max(a3)); k1=selrel/k;
    aa3=k1*a/max(a)+1*a3/max(a3); maxaa3=max(aa3); aa3=aa3/maxaa3;
    plot(1e7./X,k1*a/max(a)/maxaa3,'color',olive,'linewidth',2)
    plot(1e7./X,1*a3/max(a3)/maxaa3,'color','m','linewidth',2)
    plot(1e7./X,aa3,'k','linewidth',2)
    plot(nmditSoret,ditSoretNorm,'color',0.6*ones(1,3),'linewidth',2)
    text(405,0.8,['{\ita}:{\ita}_3=',num2str(selrel,'%.2f'),':1'],'fontsize',16)
%     text(405,0.8,['{\ita}:{\ita}_3=',num2str(a_to_a3,'%.2f'),':1'],'fontsize',16)
    set(gca,'fontsize',16,'Xlim',[400 480],'Ylim',[0 1.1])
    box on
    title('Reconstructed {\itaa}_3 Soret','fontsize',16)
    xlabel('Wavelength (nm)','fontsize',16)
    ylabel('Absorption (r.u.)','fontsize',16)
    %     l
    %     lres=reshape(l,3,N);
    %
    al=ws*2*sqrt(2*log(2));
    wls=1e7./(1e7./ppos-al/2)-1e7./(1e7./ppos+al/2)
    ws
%     save macheck L
    sq1=trapz(X,S{1}(:,1))/trapz(X,S{1}(:,3))/(trapz(X,S{2}(:,1))/trapz(X,S{2}(:,3)))
    sq2=trapz(X,S{1}(:,1))/trapz(X,S{1}(:,4))/(trapz(X,S{2}(:,1))/trapz(X,S{2}(:,4)))
    sq3=trapz(X,S{1}(:,1))/trapz(X,S{1}(:,3))/(trapz(X,S{3}(:,1))/trapz(X,S{3}(:,3)))
    sq4=trapz(X,S{1}(:,1))/trapz(X,S{1}(:,4))/(trapz(X,S{3}(:,1))/trapz(X,S{3}(:,4)))
    
    sq5=trapz(X,S{1}(:,2))/trapz(X,S{1}(:,5))/(trapz(X,S{2}(:,2))/trapz(X,S{2}(:,5)))
    sq6=trapz(X,S{1}(:,2))/trapz(X,S{1}(:,6))/(trapz(X,S{2}(:,2))/trapz(X,S{2}(:,6)))
    sq7=trapz(X,S{1}(:,2))/trapz(X,S{1}(:,5))/(trapz(X,S{3}(:,2))/trapz(X,S{3}(:,5)))
    sq8=trapz(X,S{1}(:,2))/trapz(X,S{1}(:,6))/(trapz(X,S{3}(:,2))/trapz(X,S{3}(:,6)))
    %     ctrl1=ppos(3)-ppos(4)
    %     if deg
    %         hemeb=S{4}(:,[1 3]);
    %     else
    %         hemeb=S{4}(:,[1 3 4]);
    %     end
    %     save hemeb X hemeb
    %     pause
%     numel(SP)
%     [h,p,stats] = chi2gof(SP-smodel,'NBins',numel(SP),'NParam',numel(L))
    p = chi2cdf(chi,degrees_of_freedom,'upper');
end
% pause

