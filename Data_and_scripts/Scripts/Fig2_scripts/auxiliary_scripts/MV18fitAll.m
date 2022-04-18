function err = MV18fitAll(x0)

global n y

n1 = x0(1);
s1 = x0(2);
h1 = x0(3);

n2 = x0(4);
s2 = x0(5);
h2 = x0(6);

n3 = x0(7);
s3 = x0(5);
h3 = x0(8);

n4 = x0(9);
s4 = x0(10);
h4 = x0(11);

n5 = x0(12);
s5 = x0(13);
h5 = x0(14);

n6 = x0(15);
s6 = x0(16);
h6 = x0(17);

n7 = x0(18);
s7 = x0(16);
h7 = x0(19);

sh1 = h1*exp(-(n-n1).^2./(2*s1^2));
sh2 = h2*exp(-(n-n2).^2./(2*s2^2));
sh3 = h3*exp(-(n-n3).^2./(2*s3^2));
sh4 = h4*exp(-(n-n4).^2./(2*s4^2));
sh5 = h5*exp(-(n-n5).^2./(2*s5^2));
sh6 = h6*exp(-(n-n6).^2./(2*s6^2));
sh7 = h7*exp(-(n-n7).^2./(2*s7^2));

S=sh1+sh2+sh3+sh4+sh5+sh6+sh7;

err=sqrt(sum((y-S).^2));