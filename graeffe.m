function [k,Q] = graeffe(p,v) %p-polynom (vektor), v-pocet iteracii
if p(1)~= 1
    p= p./p(1);
end
n=length(p);
Q(1,1:n)=p;
for i=2:v
    vect=Q(i-1,1:n);
    Q(i,1:n)=iter(vect);
end
y=zeros(1,n-1);
for i=1:n-1
    y(i)=-(Q(v,i+1))/Q(v,i);
end
for l=1:n-1
    rts=kthroots(y(l),2^(v-1));
    min = abs(polyval(p, rts(1)));
    best = 1;
    for t=2:2^(v-1)
        a=rts(t);
        if abs(polyval(p,a)) < min               
            best = t;
            min = abs(polyval(p,a));
        end
    end
    y(l)=rts(best);
end
k=y;
end

function q = iter(p)
n=length(p);
if mod(n,2)==0
    pe = zeros(1,n/2);
    po = zeros(1,n/2);
    for i=1:n/2
        pe(i)=p(2*i);
        po(i)=p((2*i)-1);
    end
else
    pe = zeros(1,(n+1)/2);
    po = zeros(1,(n-1)/2);
    for i=1:(n-1)/2
        po(i)=p(2*i);
        pe(i)=p((2*i)-1);
    end
    pe((n+1)/2)=p(n);
end
q= ((-1)^(n+1))*sum_poly(conv(pe,pe),-conv([1 0],conv(po,po)));
end

function x = sum_poly(a, b)
m=length(a);
n=length(b);
if m>n 
    s=m-n;
    d=zeros(1,s);
    x=[d,b]+a;
else 
    s=n-m;
    d=zeros(1,s);
    x=[d,a]+b;
end
end

function r=kthroots(x,k)
r=zeros(1,k);
a=abs(x)^(1/k);
ang=angle(x)/k;
for l=0:k-1
    var = (2*pi*l+ang)/k;
    r(l+1)=a*exp(1)^(1i*var);
end
end    
