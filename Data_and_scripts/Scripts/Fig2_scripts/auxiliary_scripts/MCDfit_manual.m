% Load and decompose MCD spectrum into 3 gaussians

load ..\mat\mcd %load a2+a32+ MCD spectrum
x=mcd(:,1);
n=1e7./mcd(:,1);

n1 = 1e7/594;
w1 = 240;
h1 = -0.66;

n2 = 1e7/605.5;
w2 = w1;
h2 = -h1;

n3 = 1e7/574;
w3 = 260;
h3 = -0.25;

s1 = h1*exp(-(n-n1).^2./(2*w1^2));
s2 = h2*exp(-(n-n2).^2./(2*w2^2));
s3 = h3*exp(-(n-n3).^2./(2*w3^2));
MCD=mcd(:,2);

MCDdata=[x s1+s2+s3 s1 s2 s3 MCD./(max(MCD)-min(MCD))];