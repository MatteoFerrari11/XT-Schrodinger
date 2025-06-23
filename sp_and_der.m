function ders = sp_and_der(p,knots,ind_DoF,t_p,k)

%input:
%   p --> splines degree
%
%   knots --> knots vector, in maximal regularity we repeat p+1 the first,
%   then Nt-1 points, and finally p+1 the last. The size is 2p+2+N-1=N+1+2p
%
%   ind_DoF --> which basis function I want to evaluate, from 0 to N+p-1
%
%   t_p --> the point where I evaluate the basis function
%
%   k --> the number of derivative of the basis function I want to compute
%
%
%output:
%   ders --> the evaluation of the k-th derivative of the basis function 
%   with index ind_DoF in t_p 

%the routine is based on the Algorithm A2.4
%there is a small typo in the NURBS BOOK
%Uright = knots(ind_DoF+j+p+jj-kk+2);
%and not  Uright = knots(ind_DoF+j+p+jj+2)

if t_p<knots(ind_DoF+1) || t_p>=knots(ind_DoF+p+2)
    ders = 0;
    return
end

for j = 0:p
    if t_p >= knots(ind_DoF+j+1) && t_p< knots(ind_DoF+j+2)
        N(j+1,1) = 1;
    else
        N(j+1,1) = 0;
    end
end

for kk = 1:p
    if N(1,kk) == 0
        saved = 0;
    else
        saved = (t_p-knots(ind_DoF+1))*N(1,kk)/(knots(ind_DoF+kk+1)-knots(ind_DoF+1));
    end
    for j = 0:p-kk
        Uleft = knots(ind_DoF+j+2);
        Uright = knots(ind_DoF+j+kk+2);
        if N(j+2,kk) == 0
            N(j+1,kk+1) = saved;
            saved = 0;
        else
            temp = N(j+2,kk)/(Uright-Uleft);
            N(j+1,kk+1) = saved + (Uright-t_p)*temp;
            saved = (t_p-Uleft)*temp;
        end
    end
end

ders(1)= N(1,p+1);

for kk = 1:k
    for j = 0:kk
        ND(j+1) = N(j+1,p-kk+1);
    end
    for jj = 1:kk
        if ND(1) == 0
            saved = 0;
        else
            saved = ND(1)/(knots(ind_DoF+p-kk+jj+1)-knots(ind_DoF+1));
        end
        for j = 0:kk-jj
            Uleft = knots(ind_DoF+j+2);
            Uright = knots(ind_DoF+j+p+jj-kk+2);
            if ND(j+2) == 0
                ND(j+1) = (p-kk+jj)*saved;
                saved = 0;
            else
                temp = ND(j+2)/(Uright-Uleft);
                ND(j+1) = (p-kk+jj)*(saved-temp);
                saved = temp;
            end
        end
    end
    ders(kk+1) = ND(1);
end

ders = ders(end);
