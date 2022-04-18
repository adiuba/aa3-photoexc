% This script performs the integration of mixed-valence CN +/- Ca2+
% difference spectrum and its decomposition into 7 gaussians
% 
% clear
% close all
load ..\mat\MV18_data % mixed-valence CN +/- Ca2+ difference spectrum
MV0=MV18;

global n y
l=zeros(length(nmMV18),1);

%-----parameters from 31 March 2021-------
n0=650;
i1=1:find(nmMV18<n0,1); i2=find(nmMV18<n0,1):length(nmMV18);
a1=0.00/38; b1=-a1*n0; l(i1)=a1*nmMV18(i1)+b1;
hlim=650;
IV=find(nmMV18<=hlim,1):find(nmMV18>400,1,'last');
%-----------------------------------------
s=0.2*exp(-(1e7/620-1e7./nmMV18).^2/(2*850^2)); %%reference 02042021
a2=0.013/38; b2=-a2*n0; l(i2)=a2*nmMV18(i2)+b2; %%reference 02042021

baseline=-(s+l);
MV18=MV18+s+l;
IS=find(nmMV18<480,1):find(nmMV18>400,1,'last');
nmMV18S=nmMV18(IS);
nmMV18V=nmMV18(IV);
MV18S=MV18(IS)-0.1;
MV18V=MV18(IV);

nfull=1e7./nmMV18V(1:end); %reference
aS=zeros(length(MV18V)-1,1);

for k=2:length(MV18V)
    aS(k-1)=trapz(nfull(1:k),MV18V(1:k));
end
yfull=aS/max(aS);
nm1=nmMV18V(1:end-1); nn1=1e7./nm1;
y=aS(nm1<620)/max(aS(nm1<620));
n=1e7./nm1(nm1<620);

s1LB=440; s1UB=550;
sSLB=250; sSUB=400;
sb1LB=265; sb1UB=600;
sb2LB=265; sb2UB=600;
svLB=140; svUB=400;
slLB=265; slUB=400;
h1LB=0.01; h1UB=0.1;
hLB=0.01; hUB=0.72;
hbUB=0.035;
hv1LB=0.01;

n1LB=1e7/421; n1UB=1e7/420;
n2LB=1e7/442.3; n2UB=1e7/435;
n3LB=1e7/455; n3UB=1e7/450;
n4LB=1e7/519; n4UB=1e7/517;
n5LB=1e7/560; n5UB=1e7/553;
n6LB=1e7/596; n6UB=1e7/593;
n7LB=1e7/608; n7UB=1e7/600;
n8LB=1e7/640; n8UB=1e7/620;

LB1=[n1LB s1LB h1LB]; UB1=[n1UB s1UB h1UB];
LB2=[n2LB sSLB hLB]; UB2=[n2UB sSUB hUB];
LB3=[n3LB hLB]; UB3=[n3UB hUB];
LB4=[n4LB sb1LB hLB]; UB4=[n4UB sb1UB hUB];
LB5=[n5LB sb2LB hLB]; UB5=[n5UB sb2UB hUB];
LB6=[n6LB svLB hv1LB]; UB6=[n6UB svUB hUB];
LB7=[n7LB hLB]; UB7=[n7UB hUB];

LB=[LB1 LB2 LB3 LB4 LB5 LB6 LB7];
UB=[UB1 UB2 UB3 UB4 UB5 UB6 UB7];

x0=LB+(UB-LB).*rand(1,numel(LB));
options.MaxFunEvals = 10000;
[r,fval,exitflag] = fmincon(@MV18fitAll,x0,[],[],[],[],LB,UB,[],options);

res=[1e7/r(1) r(2) r(3) 1e7/r(4) r(5)  r(6) 1e7/r(7) r(8)...
    1e7/r(9) r(10) r(11) 1e7/r(12) r(13)  r(14)...
    1e7/r(15) r(16) r(17) 1e7/r(18) r(19)];

n1 = r(1);
s1 = r(2);
h1 = r(3);

n2 = r(4);
s2 = r(5);
h2 = r(6);

n3 = r(7);
s3 = r(5);
h3 = r(8);

n4 = r(9);
s4 = r(10);
h4 = r(11);

n5 = r(12);
s5 = r(13);
h5 = r(14);

n6 = r(15);
s6 = r(16);
h6 = r(17);

n7 = r(18);
s7 = r(16);
h7 = r(19);

ss1 = h1*exp(-(nn1-n1).^2./(2*s1^2));
ss2 = h2*exp(-(nn1-n2).^2./(2*s2^2));
ss3 = h3*exp(-(nn1-n3).^2./(2*s3^2));
ss4 = h4*exp(-(nn1-n4).^2./(2*s4^2));
ss5 = h5*exp(-(nn1-n5).^2./(2*s5^2));
ss6 = h6*exp(-(nn1-n6).^2./(2*s6^2));
ss7 = h7*exp(-(nn1-n7).^2./(2*s7^2));

SS=ss1+ss2+ss3+ss4+ss5+ss6+ss7;


rms=sqrt(sum(MV0.^2)/numel(MV0));
drms=sqrt(sum((l+s).^2)/numel(MV18));
rms_err=drms/rms*100;

save ..\mat\MV18all_res_v2_02042021 res n1 s1 h1 n2 s2 h2 n3 s3 h3 n4 s4 h4 n5 s5 h5...
    n6 s6 h6 n7 s7 h7...
    ss1 ss2 ss3 ss4 ss5 ss6 ss7 SS nmMV18V MV18V y yfull baseline
