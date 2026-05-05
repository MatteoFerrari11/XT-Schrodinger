% This code is associated with the paper 
% "A matrix-based approach to the stability of a space-time isogeometric method for the linear Schrödinger equation"
%
% by M. Ferrari and S. Gómez

% This code reproduces Figure 6

clc
clear
close all
format long

p1 = 2; % Degree of spline functions
p2 = 5; % Degree of spline functions
cont_p = 0;

T = 1; % Final time
a = -3; % Space interval [a,b]
b = 3;

N = 3; % Resolution level (time and space intervals 2^(N+1))

Nq = 16; % Number of quadrature points

syms t x

%quantum harmonic oscillator
w = 10;
F = @(t,x) 0*t*x;
n = 2;
Psi = @(t,x) 1/sqrt(2^n*factorial(n))*(w/pi)^(1/4)*exp(-1/2*(2*n+1)*1i*w*t)*(hermiteH(n,sqrt(w)*x).*exp(-1/2*(w*x.^2)));
dtPsi = @(t,x) -1i*(2*n+1)*w/2* Psi(t,x);
dxPsi = @(t,x) 1/sqrt(2^n*factorial(n))*(w/pi)^(1/4)*sqrt(w)*exp(-1i*(2*n+1)*w*t/2).* ...
    ((2*n*hermiteH(n-1,sqrt(w)*x)-sqrt(w)*x.*hermiteH(n,sqrt(w)*x)).*exp(-w*x.^2/2));

Psi_0 = @(x) Psi(0,x);
Psi_ta = @(t) Psi(t,a);
Psi_tb = @(t) Psi(t,b);
V = @(x) -w^2*x.^2/2;

for p = p1:p2

    cont_p = cont_p + 1;

    N_x = 2^(N+1)*6;
    N_t = 2^(N+1);

    knots_t = [zeros(p,1); linspace(0,T,N_t+1)'; ones(p,1)*T]';
    knots_x = [a*ones(p,1); linspace(a,b,N_x+1)'; ones(p,1)*b]';

    t = ((0:N_t)/N_t)*T;
    x = a+((0:N_x)/N_x)*(b-a);
    h_t = t(2);
    h_x = x(2)-x(1);

    B_t = mat_splines(N_t,p,1,1,knots_t,0,T);
    C_t = mat_splines(N_t,p,0,1,knots_t,0,T);
    M_t = mat_splines(N_t,p,0,0,knots_t,0,T);
    B_x = mat_splines(N_x,p,1,1,knots_x,a,b);
    M_x = mat_splines(N_x,p,0,0,knots_x,a,b);
    M_x2 = mat_splines(N_x,p,0,0,knots_x,a,b);
    M_v_x = mat_splines_V(N_x,p,0,0,knots_x,a,b,V);

    F_vec = zeros((N_x+p)*(N_t+p),1);
    for j_t = 1 : N_t+p
        for j_x = 1 : N_x+p
            for k_t = max(1,j_t-p):min(j_t,N_t)
                phi_t = zeros(Nq,1);
                for k_x = max(1,j_x-p):min(j_x,N_x)
                    phi_x = zeros(Nq,1);
                    [xs_t,ps_t] = lgwt(Nq,t(k_t),t(k_t+1));
                    [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
                    for iii = 1:Nq
                        phi_t(iii) = sp_and_der(p,knots_t,j_t-1,xs_t(iii),0);
                        phi_x(iii) = sp_and_der(p,knots_x,j_x-1,xs_x(iii),0);
                    end
                    F_vec((j_t-1)*(N_x+p)+j_x) = F_vec((j_t-1)*(N_x+p)+j_x) + (phi_t.*ps_t)'*F(xs_t,xs_x')*(phi_x.*ps_x);
                end
            end
        end
    end

    F_proj_x = zeros(N_x+p,1);
    for j_x = 1 : N_x+p
        for k_x = max(1,j_x-p):min(j_x,N_x)
            phi_x = zeros(Nq,1);
            [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
            for iii = 1:Nq
                phi_x(iii) = sp_and_der(p,knots_x,j_x-1,xs_x(iii),0);
            end
            F_proj_x(j_x) = F_proj_x(j_x) + sum(Psi_0(xs_x).*phi_x.*ps_x);
        end
    end

    F_proj_ta = zeros(N_t+p,1);
    F_proj_tb = zeros(N_t+p,1);
    for j_t = 1 : N_t+p
        for k_t = max(1,j_t-p):min(j_t,N_t)
            phi_t = zeros(Nq,1);
            [xs_t,ps_t] = lgwt(Nq,t(k_t),t(k_t+1));
            for iii = 1:Nq
                phi_t(iii) = sp_and_der(p,knots_t,j_t-1,xs_t(iii),0);
            end
            F_proj_ta(j_t) = F_proj_ta(j_t) + sum(Psi_ta(xs_t).*phi_t.*ps_t);
            F_proj_tb(j_t) = F_proj_tb(j_t) + sum(Psi_tb(xs_t).*phi_t.*ps_t);
        end
    end

    Psi_proj_ta = M_t\F_proj_ta;
    Psi_proj_tb = M_t\F_proj_tb;
    Psi_proj_x = M_x\F_proj_x;

    S = 1i*kron(sparse(B_t),sparse(M_x)) + kron(sparse(C_t),sparse(-1/2*B_x+M_v_x));

    Psi_0_ab_proj = zeros(size(F_vec));
    for j_x = 1
        for j_t = 1 : N_t+p
            Psi_0_ab_proj((j_t-1)*(N_x+p)+j_x) = Psi_proj_ta(j_t);
        end
    end
    for j_x = N_x+p
        for j_t = 1 : N_t+p
            Psi_0_ab_proj((j_t-1)*(N_x+p)+j_x) = Psi_proj_tb(j_t);
        end
    end
    for j_t = 1
        for j_x = 1 : N_x+p
            Psi_0_ab_proj((j_t-1)*(N_x+p)+j_x) = Psi_proj_x(j_x);
        end
    end

    F_vec = F_vec-S*Psi_0_ab_proj;

    %removing initial condition in time and boundary in space
    cont_v = [];
    for j_t = 1 : N_t+p-1
        for j_x = 2 : N_x+p-1
            cont_v = [cont_v (j_t-1)*(N_x+p)+j_x];
        end
    end
    F_vec = F_vec(cont_v);

    B_t = B_t(1:end-1,2:end); %zero initial condition
    C_t = C_t(1:end-1,2:end); %zero initial condition
    B_x = B_x(2:end-1,2:end-1); %zero boundary conditions
    M_x = M_x(2:end-1,2:end-1); %zero boundary conditions
    M_v_x = M_v_x(2:end-1,2:end-1); %zero boundary conditions

    S = 1i*kron(sparse(B_t),sparse(M_x)) + kron(sparse(C_t),sparse(-1/2*B_x+M_v_x));
    Psi_app_coeff = S\F_vec;

    %we add zero values at zero in time and at the boundary in space
    Psi_app_coeff_w_b = zeros((N_x+p)*(N_t+p),1);
    cont_v = 0;
    for j_t = 1 : N_t+p
        for j_x = 1 : N_x+p
            if j_t == 1
                Psi_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = 0;
            elseif j_x == 1
                Psi_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = 0;
            elseif j_x == N_x+p
                Psi_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = 0;
            else
                cont_v = cont_v+1;
                Psi_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = Psi_app_coeff(cont_v);
            end
        end
    end
    Psi_app_coeff = Psi_app_coeff_w_b;

   
    tic
    Nx_p = N_x + p;
    Nt_p = N_t + p;
    n_basis = Nx_p * Nt_p;

    t_plot = linspace(2*eps,1-2*eps,100);
    
    mass_2 = zeros(length(t_plot),1);
    mass_2_V = zeros(length(t_plot),1);
    energy_2 = zeros(length(t_plot),1);
    for ind_t = 1:length(t_plot)
        t_curr = t_plot(ind_t);

        for int_x = 1:N_x
            [xs_x, ps_x] = lgwt(p+2, x(int_x), x(int_x+1));
            loc_sum = zeros(size(ps_x));
            loc_sum_d = zeros(size(ps_x));

            for iii = 1:length(xs_x)
                phi_x_vals = zeros(1, Nx_p);
                dphi_x_vals = zeros(1, Nx_p);

                for j_x = 1:Nx_p
                    phi_x_vals(j_x) = sp_and_der(p, knots_x, j_x-1, xs_x(iii), 0);
                    dphi_x_vals(j_x) = sp_and_der(p, knots_x, j_x-1, xs_x(iii), 1);
                end

                for j_t = 1:Nt_p
                    phi_t = sp_and_der(p, knots_t, j_t-1, t_curr, 0);
                    jt_offset = (j_t - 1) * Nx_p;

                    for j_x = 1:Nx_p
                        idx = jt_offset + j_x;
                        psi_val = Psi_0_ab_proj(idx) + Psi_app_coeff(idx);
                        phi_val = phi_t * phi_x_vals(j_x);
                        dphi_val = phi_t * dphi_x_vals(j_x);

                        loc_sum(iii) = loc_sum(iii) + psi_val * phi_val;
                        loc_sum_d(iii) = loc_sum_d(iii) + psi_val * dphi_val;
                    end
                end
            end

            integrand = abs(loc_sum).^2 .* ps_x;
            mass_2(ind_t) = mass_2(ind_t) + sum(integrand);
            mass_2_V(ind_t) = mass_2_V(ind_t) + sum(integrand .* V(xs_x));
            energy_2(ind_t) = energy_2(ind_t) + sum(abs(loc_sum_d).^2 .* ps_x);
        end
    end
    toc

    mass(cont_p,:) = mass_2/2;
    energy(cont_p,:) = -energy_2/2 + mass_2_V;

    figure(1)
    semilogy(t_plot,abs(mass(cont_p,:)-mass(cont_p,1)),'o-','LineWidth',2)
    hold on

    figure(2)
    semilogy(t_plot,abs(energy(cont_p,:)-energy(cont_p,1)),'o-','LineWidth',2)
    hold on
end

if cont_p>1
    figure(1)
    legend('p=2','p=3','p=4','p=5')
    title('mass')

    figure(2)
    legend('p=2','p=3','p=4','p=5')
    title('energy')
end
