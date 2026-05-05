function [x,w]=lgwt_vpa(N,a,b)

global SIZE

N=N-1;
N1=N+1; N2=N+2;

a = vpa(a);
b = vpa(b);
xu = vpa(linspace(-1,1,N1)');

my_pi = vpa(pi);
y = cos((2*vpa((0:N)')+1)*my_pi/(2*N+2)) + (vpa('0.27')/N1)*sin(my_pi*xu*N/N2);

L = vpa(zeros(N1,N2));
Lp = vpa(zeros(N1,N2));

y0 = vpa(2);
tol = vpa('10')^(-SIZE); 

while max(abs(y-y0)) > tol
    
    L(:,1) = vpa(1);
    Lp(:,1) = vpa(0);
    
    L(:,2) = y;
    Lp(:,2) = vpa(1);
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

x=(a*(1-y)+b*(1+y))/2;      
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;