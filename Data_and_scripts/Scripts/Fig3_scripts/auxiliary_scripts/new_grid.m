function y2 = new_grid(a,b,y)

% a - old XData;
% b - new XData;

[a,I]=sort(a);
b=sort(b);
y=y(I);

c=sort(unique([a;b]));
sc=size(c);
dc = diff(c);
K = diff(y)./diff(a);

la=length(a);
f=[];
for i = 1:length(a)
    ai=a(i);
    f=[f; find(c==ai)];
    k=1;
%     pause
end
size(f);
df=diff(f);
sdf=size(df);
sumdf=sum(df)
aend=a(end);
acend=c(f(end));

K2=[];
for i = 1:length(df)
    K2=[K2; K(i)*ones(df(i),1)];
end

sK2=size(K2);

Y2=[y(1); K2(1)*dc(1)];
for i = 2:length(K2)
    Y2=[Y2; Y2(end)+K2(i)*dc(i)];
end
size(Y2);

fb=[];
for i = 1:length(b)
    bi=b(i);
    fb=[fb; find(c==bi)];
end

fb(end);
y2=Y2(fb);
% y2=y2(end:-1:1);
    