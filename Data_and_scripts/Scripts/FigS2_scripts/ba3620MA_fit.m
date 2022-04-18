function err = ba3620MA_fit(l)

global ba3x ba3y

L=l;

ppos=l(1:3); l(1:3)=[];
wpos=l(1:2); l(1:2)=[];
apos=l(1:3); l(1:3)=[];

pneg=l(1:2); l(1:2)=[];
wneg=l(1:2); l(1:2)=[];
aneg=l(1:2);

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

S=sum(s1,2);
err=sqrt(sum((ba3y-S).^2))/max(ba3y);