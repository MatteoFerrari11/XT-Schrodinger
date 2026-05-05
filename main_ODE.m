% This code is associated with the paper 
% "A matrix-based approach to the stability of a space-time isogeometric method for the linear Schrödinger equation"
%
% Authors: M. Ferrari and S. Gómez
%
% This code computes an approximate solution of the ODE associated to the Schrödinger problem:
%
%     i d_t u(t) + u(t) = f(t),  t in (0,T)
%     u(0) = u0, 
%
% using maximal regularity splines and the variational formulation:
%
% find u_h in u_0 + S_{h,0,\bullet}^p(0,T) 
%
%     i (d_t u_h, d_t w_h)_{L^2(0,T)} + mu * (u_h, d_t w_h)_{L^2(0,T)} = (f, d_t w_h)_{L^2(0,T)}
%
% for all w_h in S_{h,\bullet,0}^p(0,T).
%
% S_{h,0,\bullet}^p(0,T) is the spline space of maximal regularity p-1 and
% polynomial degree p over a uniform mesh with mesh size h and zero
% initial condition
%
% S_{h,0,\bullet}^p(0,T) is the spline space of maximal regularity p-1 and
% polynomial degree p over a uniform mesh with mesh size h and zero
% final condition
%
% We calculate the errors in the norms L^2 and H1, L^inf e W^{1,inf}

clc
clear
close all
format long

p = 2; % Degree of spline functions
T = 1; % Final time
mu = 100000; % Parameter in the differential equation

N_plot = 10000; % Number of points for reconstructing the solution
Nq = 64; % Number of quadrature points
N1 = 3; % Starting resolution level (time intervals 2^(N1+1))
N2 = 3; % Ending resolution level (time intervals 2^(N2+1))

err_L2_app = zeros(1,N2-N1+1);
err_Linf_app = zeros(1,N2-N1+1);
err_W1inf_app = zeros(1,N2-N1+1);
err_H1_app = zeros(1,N2-N1+1);

syms t
u = exp(-t); % Analytical solution
du = matlabFunction(diff(u));
f = matlabFunction(1i*diff(u,1)+mu/2*u); 
u = matlabFunction(u); 
u0 = u(0);

cont = 0;

for ii = N1:N2
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);
    
    B = mat_splines(Nt,p,1,1,knots,0,T);
    C = mat_splines(Nt,p,1,0,knots,0,T);
    
    F = zeros(Nt+p,1);
    for jj = 1:Nt+p
        phi = zeros(Nq,1);
        for k = max(1,jj-p):min(jj,Nt)
            [xs,ps] = lgwt(Nq,t(k),t(k+1));
            for iii = 1:Nq
                phi(iii) = sp_and_der(p,knots,jj-1,xs(iii),1);
            end
            F(jj) = F(jj) + sum(f(xs).*phi.*ps) - mu/2*u0*sum(phi.*ps);
        end
    end

    S = 1i*B - mu/2*C;
    S = S(1:end-1,2:end); %zero initial conditions
    F = F(1:end-1); %zero initial conditions

    u_app_coeff = [0; S\F];
  
    x_plot = linspace(0.0001, T-0.0001, N_plot);
    u_ex = u(x_plot);
    du_ex = du(x_plot);
    u_app = zeros(size(u_ex));
    du_app = zeros(size(u_ex));

    for i_plot = 1:N_plot
        for ind = 1:Nt + p
            u_app(i_plot) = u_app(i_plot) + u_app_coeff(ind) * sp_and_der(p, knots,ind-1,x_plot(i_plot),0);
            du_app(i_plot) = du_app(i_plot) + u_app_coeff(ind) * sp_and_der(p,knots,ind-1,x_plot(i_plot),1);
        end
    end
    u_app = u_app+u0;

    if N1 == N2
        % Plot u_app and u_ex (solution)
        figure(1)
        plot(x_plot,real(u_app),'r','LineWidth',1.5)
        hold on
        plot(x_plot,real(u_ex),'b--','LineWidth',1.5)
        hold off
        title(['$u_h(t)$ and $u(t)$ for maximal regularity splines with p=' num2str(p)],'FontSize',14,'Interpreter','latex')
        xlabel('$t$','FontSize',12,'Interpreter','latex')
        legend({'$u_h(t)$', '$u(t)$'},'FontSize',16,'Location','best',...
            'Interpreter','latex')
        grid on

        % Plot du_app and du_ex (first derivative)
        figure(2)
        plot(x_plot,real(du_app),'r','LineWidth',1.5)
        hold on
        plot(x_plot,real(du_ex),'b--','LineWidth',1.5)
        hold off
        title(['$\partial_t u_h(t)$ and $\partial_t u(t)$ for maximal regularity splines with p=' num2str(p)],'FontSize',14,'Interpreter','latex')
        xlabel('$t$','FontSize',12,'Interpreter','latex')
        legend({'$\partial_t u_h(t)$','$\partial_t u(t)$'},'FontSize', 16,...
            'Location','best','Interpreter','latex')
        grid on

        % Plot the absolute error in u
        figure(3)
        loglog(x_plot,abs(real(u_app)-real(u_ex)),'r','LineWidth',1.5)
        title(['Absolute error $|u_h(t)-u(t)|$ for maximal regularity splines with p=' num2str(p)],'FontSize',14,'Interpreter','latex')
        xlabel('$t$','FontSize',12,'Interpreter','latex')
        ylabel('$|u_h(t)-u(t)|$','FontSize',12,'Interpreter','latex')
        grid on
    end

    err_L2_app(cont) = sqrt(T/length(x_plot)*sum(abs(u_app-u_ex).^2))/sqrt(T/length(x_plot)*sum(abs(u_ex).^2));
    err_Linf_app(cont) = max(abs(u_app-u_ex))/max(abs(u_ex));
    err_W1inf_app(cont) = max(abs(du_app-du_ex))/max(abs(du_ex));
    err_H1_app(cont) = sqrt(T/length(x_plot)*sum(abs(du_app-du_ex).^2))/sqrt(T/length(x_plot)*sum(abs(du_ex).^2));
end

if N1 < N2
    plot_error({}, err_L2_app,T./2.^(N1:N2),1);
    plot_error({}, err_H1_app,T./2.^(N1:N2),1);
    plot_error({}, err_Linf_app, T ./ 2.^(N1:N2),1);
    plot_error({'L2 error', 'H1 error','Linf error', 'W1inf error',}, err_W1inf_app, T ./ 2.^(N1:N2),1);
    title(['Various errors for maximal regularity splines with p=' num2str(p)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
end
