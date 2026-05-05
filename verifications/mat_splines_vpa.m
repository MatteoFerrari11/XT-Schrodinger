function K = mat_splines_vpa(Nt,p,k1,k2,knots,a,b)

global xs_g ws_g

x = xs_g;
w = ws_g;

a = vpa(a);
b = vpa(b);

t = a + (vpa(0:Nt)/Nt)*(b-a);
K = vpa(zeros(Nt+p,Nt+p));

for k = 1:Nt
    span_idx = k + p; 
    
    xs = (t(k+1)-t(k))/vpa(2)*x + (t(k+1)+t(k))/vpa(2);
    ws = (t(k+1)-t(k))/vpa(2)*w;
    
    active_dofs = k:(k+p); 
    
    phi_1 = vpa(zeros(p+1,p+1));
    phi_2 = vpa(zeros(p+1,p+1));
    
    for f_idx = 1:(p+1)
        ind_DoF = active_dofs(f_idx) - 1; 
        for q = 1:p+1
            phi_1(q, f_idx) = sp_and_der_vpa(p,knots,ind_DoF,xs(q),k1,span_idx);
            phi_2(q, f_idx) = sp_and_der_vpa(p,knots,ind_DoF,xs(q),k2,span_idx);
        end
    end
    
    for f_idy = 1:(p+1)
        for f_idx = 1:(p+1)
            jj = active_dofs(f_idy);
            ii = active_dofs(f_idx);
            local_integral = sum(phi_2(:,f_idy) .* phi_1(:,f_idx).*ws);
            K(jj,ii) = K(jj,ii) + local_integral;
        end
    end
end
end