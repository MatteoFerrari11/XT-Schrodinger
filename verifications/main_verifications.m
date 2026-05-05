% This code is associated with the paper 
% "A matrix-based approach to the stability of a space-time isogeometric method for the linear Schrödinger equation"
%
% by M. Ferrari and S. Gómez

% This code reproduces Figure 1

clc
clear
close all
format long

global SIZE
SIZE = 1000;

digits(SIZE)

p = 2;
N1 = 2;
N2 = 3;

global xs_g ws_g
Nq = p+1; 
[xs_g,ws_g] = lgwt_vpa(Nq,vpa(-1),vpa(1));

T = vpa('1');

cont = 0;

for ii = N1:N2
    cont = cont+1;

    Nt = 2^(ii+1);

    zero_vec = vpa(zeros(p,1));
    one_vec  = vpa(ones(p,1));

    mid_knots = (vpa(0:Nt)'/Nt)*T;

    knots = [zero_vec;mid_knots;one_vec*T]';
    h = T/vpa(Nt);

    B = mat_splines_vpa(Nt,p,1,1,knots,vpa('0'),T);
    C = mat_splines_vpa(Nt,p,1,0,knots,vpa('0'),T);

    B = B(1:end-1,2:end)*h;
    C = C(1:end-1,2:end);

    M_vpa = B\C;
    eig_n = eig(M_vpa);

    cmap = lines(N2-N1+1);
    color = cmap(cont,:);
    plot(double(real(eig_n)), double(imag(eig_n)), 'o', ...
         'MarkerFaceColor', color, 'MarkerEdgeColor', color)
    xlabel('$\Re(\lambda)$', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$\Im(\lambda)$', 'Interpreter', 'latex', 'FontSize', 14)
    title(sprintf(['Nt = %d, $p=%d$, ' ...
        '$N_1=%d$, $N_2=%d$, $T=%s$'], Nt, p, N1, N2, char(T)), ...
        'Interpreter', 'latex', 'FontSize', 14)

    grid on
    axis equal
    box on

    ax = gca;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.TickLabelInterpreter = 'latex';
    hold on

end
