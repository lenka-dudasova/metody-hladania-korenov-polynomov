function k = aberth(p,z,iter) %p-polynom (vektor), z-pociatocne aproximacie korenov (vektor), iter-pocet iteracii
k=z;
for i=1:iter
    k=step(p,k);
end
end

function y=step(p,x)
n=length(x);
y=zeros(1,n);
for i=1:n
    s=0;
    for j=1:n
        if j~=i
            s=s+(1/(x(i)-x(j)));
        end
    end
    a=polyval(p,x(i));
    b=polyval(polyder(p),x(i));
    y(i)=x(i)-((a/b)/(1-(a/b)*s));
end
end
