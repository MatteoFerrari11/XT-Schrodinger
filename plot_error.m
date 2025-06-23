function [imm,mylegendL2] = plot_error(mylegendL2,Err,h_max,N)

figure(N)
xv = h_max;
yv = Err;

slope = [
    'for i=1:length(xv)-1;'...
    'x = sqrt(xv(i)*xv(i+1));'...
    'y = sqrt(yv(i)*yv(i+1));'...
    'm = (log(yv(i+1))-log(yv(i)))/(log(xv(i+1))-log(xv(i)));'...
    'text(x,y,sprintf(''%.2f'',(m)),''HorizontalAlignment'',''center'',''BackgroundColor'',''y'',''FontSize'',15);'...
    'end'
    ];

imm = loglog(xv,yv,'-d','LineWidth',2);
hold on
eval(slope);

grid on
legend(mylegendL2,'Location','SouthEast','FontSize',14)
drawnow
