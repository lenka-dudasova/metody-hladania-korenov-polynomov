function k = durand_kerner(p,x0,iter) %p-polynom (vektor), x0- pociatocne aproximacie korenov (vektor), iter-pocet iteracii
k=x0;
for i=1:iter
    k=step(p,k);
end
end

function y=step(p,x)
n=length(x);
y=zeros(1,n);
for i=1:n
    z=1;
    for j=1:n
        if j<i
            z=z*(x(i)-y(j));
        end
        if j>i
            z=z*(x(i)-x(j));
        end
    end
    y(i)=x(i)-(polyval(p,x(i))/z);
end
end
