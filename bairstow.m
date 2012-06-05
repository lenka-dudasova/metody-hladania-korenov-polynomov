function[R] = bairstow(n, p, q, d, a) %n-stupen polynomu, p,q-pociatocne body, d-presnost, a-polynom (vektor)
if n<4
    R=roots(a);
    return
end
[P,Q] = findpq(n, p, q, d, a);
dz=[1,-P,-Q];
R1=roots(dz);
New=deconv(a, dz);
R2=bairstow(n-2, p, q, d, New);
R=[R1;R2];
end

function[P,Q]=findpq(n,p,q,d,a)
b=zeros(n+2,1);
c=zeros(n+2,1);
b(2)=a(1);
c(2)=a(1);
for k=3:(n+1)
    b(k)=a(k-1)+p*b(k-1)+q*b(k-2);
end
b(n+2)=a(n+1)+q*b(n);
for j=3:(n-1)
    c(j)=b(j)+p*c(j-1)+q*c(j-2);
end
c(n)=b(n)+q*c(n-2);
M=p*c(n)-q*c(n-1);
W=q*b(n+1)-p*b(n+2);
N=(c(n))^2 + M*c(n-1);
dp=(c(n-1)*b(n+2)-c(n)*b(n+1))/N;
dq=(c(n-1)*W-c(n)*b(n+2))/N;
P=p+dp;
Q=q+dq;
if max(abs(dp),abs(dq))<d
    return
else
    [P, Q] = findpq(n, P, Q, d, a);
    return
end
end
