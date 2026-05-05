% This code is associated with the paper 
% "A matrix-based approach to the stability of a space-time isogeometric method for the linear Schrödinger equation"
%
% by M. Ferrari and S. Gómez

% This code reproduces Figure 4

clc
clear
close all
format long

p1 = 3; %Degree of spline functions
p2 = 5;
cont_p = 0;

T = 1; %Final time
a = -3; %Space interval [a,b]
b = 3;

N1 = 3; % Starting resolution level (time and space intervals 2^(N1+1))
N2 = 10; % Ending resolution level (time and space intervals 2^(N2+1))

N_plot = 100; % Number of points for reconstructing the solution
Nq = 16; % Number of quadrature points

err_L2_app = zeros(p2-p1+1,N2-N1+1);
err_H1_app = zeros(p2-p1+1,N2-N1+1);

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

N_t = 16;

for p = p1:p2 

    cont_p = cont_p + 1;
    cont_h = 0;

    for ii = N1:N2
        cont_h = cont_h + 1;
        N_x = 2^(ii+1)*6;

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

        t_plot = linspace(0.0001,T-0.0001,N_plot);
        x_plot = linspace(a+0.0001,b-0.0001,N_plot);

        Psi_ex = Psi(t_plot',x_plot);
        dtPsi_ex = dtPsi(t_plot',x_plot);
        dxPsi_ex = dxPsi(t_plot',x_plot);
        Psi_app = zeros(size(Psi_ex));
        Psi_app_0 = zeros(size(Psi_ex));
        dtPsi_app_0 = zeros(size(Psi_ex));
        dxPsi_app_0 = zeros(size(Psi_ex));
        dtPsi_app = zeros(size(Psi_ex));
        dxPsi_app = zeros(size(Psi_ex));

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

        for i_plot_t = 1 : N_plot
            for i_plot_x = 1 : N_plot
                for ind_t = 1 : N_t+p
                    if t_plot(i_plot_t) >= knots_t(ind_t) && t_plot(i_plot_t) < knots_t(ind_t+p+1)
                        for ind_x = 1 : N_x+p
                            if x_plot(i_plot_x) >= knots_x(ind_x) && x_plot(i_plot_x) < knots_x(ind_x+p+1)
                                Psi_app(i_plot_t,i_plot_x) = Psi_app(i_plot_t,i_plot_x)...
                                    + Psi_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                                Psi_app_0(i_plot_t,i_plot_x) = Psi_app_0(i_plot_t,i_plot_x)...
                                    + Psi_0_ab_proj((ind_t-1)*(N_x+p)+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                                dtPsi_app_0(i_plot_t,i_plot_x) = dtPsi_app_0(i_plot_t,i_plot_x)...
                                    + Psi_0_ab_proj((ind_t-1)*(N_x+p)+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                                dxPsi_app_0(i_plot_t,i_plot_x) = dxPsi_app_0(i_plot_t,i_plot_x)...
                                    + Psi_0_ab_proj((ind_t-1)*(N_x+p)+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),1);

                                dxPsi_app(i_plot_t,i_plot_x) = dxPsi_app(i_plot_t,i_plot_x)...
                                    + Psi_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),1);

                                dtPsi_app(i_plot_t,i_plot_x) = dtPsi_app(i_plot_t,i_plot_x)...
                                    + Psi_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);
                            end
                        end
                    end
                end
            end
        end
        Psi_app = Psi_app + Psi_app_0;
        dtPsi_app = dtPsi_app + dtPsi_app_0;
        dxPsi_app = dxPsi_app + dxPsi_app_0;

        %errors
        err_L2_app(cont_p,cont_h) = sqrt(sum(sum(abs(Psi_app-Psi_ex).^2))/(N_plot^2))/sqrt(sum(sum(abs(Psi_ex).^2))/(N_plot^2))
        err_H1_app(cont_p,cont_h) = (sqrt(sum(sum(abs(dtPsi_app-dtPsi_ex).^2))/(N_plot^2))...
            + sqrt(sum(sum(abs(dxPsi_app-dxPsi_ex).^2))/(N_plot^2)))./...
            (sqrt(sum(sum(abs(dtPsi_ex).^2))/(N_plot^2)) + sqrt(sum(sum(abs(dxPsi_ex).^2))/(N_plot^2)))
    end

    %error estimates
    if N1 < N2
        plot_error({}, err_L2_app(cont_p,:),1./2.^(N1+1:N2+1),1);
        plot_error({}, err_H1_app(cont_p,:),1./2.^(N1+1:N2+1),2);
        xlabel('h')
        hold on
        pause(0.01)
    end
end

