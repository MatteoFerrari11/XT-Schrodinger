function ders = sp_and_der_vpa(p,knots,ind_DoF,t_p,k, span_idx)

N_mat = vpa(zeros(p+2, p+1)); 

j_active = span_idx - ind_DoF - 1;
if j_active >= 0 && j_active <= p
    N_mat(j_active+1, 1) = vpa(1);
end

for kk = 1:p
    if N_mat(1,kk) == 0
        saved = vpa(0);
    else
        saved = (t_p-knots(ind_DoF+1))*N_mat(1,kk)/(knots(ind_DoF+kk+1)-knots(ind_DoF+1));
    end
    for j = 0:p-kk
        Uleft = knots(ind_DoF+j+2);
        Uright = knots(ind_DoF+j+kk+2);
        if N_mat(j+2,kk) == 0
            N_mat(j+1,kk+1) = saved;
            saved = vpa(0);
        else
            temp = N_mat(j+2,kk)/(Uright-Uleft);
            N_mat(j+1,kk+1) = saved + (Uright-t_p)*temp;
            saved = (t_p-Uleft)*temp;
        end
    end
end

ders_array = vpa(zeros(k+1, 1));
ders_array(1)= N_mat(1,p+1);
ND = vpa(zeros(p+2,1));

for kk = 1:k
    for j = 0:kk
        ND(j+1) = N_mat(j+1,p-kk+1);
    end
    for jj = 1:kk
        if ND(1) == 0
            saved = vpa(0);
        else
            saved = ND(1)/(knots(ind_DoF+p-kk+jj+1)-knots(ind_DoF+1));
        end
        for j = 0:kk-jj
            Uleft = knots(ind_DoF+j+2);
            Uright = knots(ind_DoF+j+p+jj-kk+2);
            if ND(j+2) == 0
                ND(j+1) = vpa(p-kk+jj)*saved;
                saved = vpa(0);
            else
                temp = ND(j+2)/(Uright-Uleft);
                ND(j+1) = vpa(p-kk+jj)*(saved-temp);
                saved = temp;
            end
        end
    end
    ders_array(kk+1) = ND(1);
end

ders = ders_array(end);