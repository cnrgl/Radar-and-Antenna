unction hplot = Plot3d(Txpos,Rxpos,Tgtposs,R,Rang)
%PLOT3D Summary of this function goes here
figure('Color','w');
Txpos = Txpos/1000;
Rxpos = Rxpos/1000;
Tgtpos = Tgtposs/1000;

plot3(Txpos(1),Txpos(2),Txpos(3),...
    'Color','k','Marker','o','MarkerSize',10,'LineStyle','none');
text(Txpos(1)+0.5,Txpos(2)+0.3,Txpos(3)+0.15,'Transmitter');
hold on;

plot3(Rxpos(1),Rxpos(2),Rxpos(3),...
    'Color','m','Marker','>','MarkerSize',10,'LineStyle','none');
text(Rxpos(1)+0.5,Rxpos(2)+0.3,Rxpos(3)+0.15,'Reciever');
hold on;

Ntgt=size(Tgtpos,2);
for n=1:Ntgt
    plot3(Tgtpos(1,n),Tgtpos(2,n),Tgtpos(3,n),...
        'Color','r','Marker','d','MarkerSize',10,'LineStyle','none');
    label=sprintf('Target %d',n);
    text(Tgtpos(1,n)+0.5,Tgtpos(2,n)+0.3,Tgtpos(3,n)+0.15,label);
    
    hold on;
end
    L=Rxpos-Txpos;
    b=sqrt((R^2)-((L(1)^2)/4));
    [ex, ey, ez]=ellipsoid(L(1)/2,L(2)/2,L(3)/2,R,b,b);
    s=surf(ex,ey,ez);
    s.EdgeColor = 'none';
    %rotate(s,[1 0 0],Rang);
    alpha 0.1
    axis([-2 20 -5 10 -5 10]);
    set(gca,'Color','none');
    xlabel('x (km)');
    ylabel('y (km)');
    zlabel('z (km)');
    title('Sistemin 3B Grafi i');
    view(50,10);
    grid on;
end
