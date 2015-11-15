function [C_e,nb]=Bezier_extract(xi)

%%Made by Jason Li

clear C Ce alphas
%xi=[0,0,1/3,2/3,1,1];
m=length(xi);
[x,y]=mode(xi);
p=y-1;

a=p+1;
b=a+1;
nb=1;
C(:,:,1)=eye(a,a);
xi(m+1)=0;

while b<m;
    C(:,:,nb+1)=eye(y,y);
    i=b;

    while and(b<m,xi(b+1)==xi(b));
        b=b+1;
    end
    mult=b-i+1;
    if mult<p;
        numer=xi(b)-xi(a);
        for j=p:-1:mult+1
            alphas(j-mult)=numer/(xi(a+j)-xi(a));
        end
        r=p-mult;
        for j=1:r;
            save=r-j+1;
            s=mult+j;
            for k=p+1:-1:s+1;
                alph=alphas(k-s);
                C(:,k,nb)=alph*C(:,k,nb)+(1.0-alph)*C(:,k-1,nb);
            end
            if b<m;
                C(save:j+save,save,nb+1)=C(p-j+1:p+1,p+1,nb);
            end
        end
        nb=nb+1;
        if b<m;
            a=b;
            b=b+1;
        end
    end
end

C_e=C(:,:,1:nb);
