function [c,ceq]=fixed_ratio_aa3_ma_623_608_570(l)


global Npos W

ppos=l(1:Npos); l(1:Npos)=[];
wpos=l(1:numel(W)); l(1:numel(W))=[];
apos=l(1:Npos*3); l(1:Npos*3)=[];
apos=reshape(apos,Npos,3);

%---heme a3----

%----ma---
ceq(1) = apos(1,1)/apos(3,1)-apos(1,2)/apos(3,2);
ceq(2) = apos(1,1)/apos(4,1)-apos(1,2)/apos(4,2);
ceq(3) = apos(1,1)/apos(3,1)-apos(1,3)/apos(3,3);
ceq(4) = apos(1,1)/apos(4,1)-apos(1,3)/apos(4,3);

%--------------------------------------------------------------------------------------------

%---heme a---

%----ma---
ceq(5) = apos(2,1)/apos(5,1)-apos(2,2)/apos(5,2);
ceq(6) = apos(2,1)/apos(6,1)-apos(2,2)/apos(6,2);
ceq(7) = apos(2,1)/apos(5,1)-apos(2,3)/apos(5,3);
ceq(8) = apos(2,1)/apos(6,1)-apos(2,3)/apos(6,3);

c=[];
end