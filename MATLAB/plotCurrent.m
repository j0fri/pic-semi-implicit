clear;

config = readmatrix("../outputs/config.txt");
T = config(1);
dtSim = config(2);
dt = config(3);
xmin = config(4);
xmax = config(5);
Nx = config(6);
ymin = config(7);
ymax = config(8);
Ny = config(9);
Ns = config(10);

if Ns > 1
    disp("Pick which component to plot of Current, 1-indexed");
end
si = 1;

dx = (xmax-xmin)/Nx;
dy = (ymax-ymin)/Ny;
x = linspace(xmin+dx/2,xmax-dx/2,Nx);
y = linspace(ymin+dy/2,ymax-dy/2,Ny);
[X,Y] = meshgrid(x,y);

xSplit = diff(x)/2;
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
%YGrid = flipud(YGrid);            
%reference: https://uk.mathworks.com/matlabcentral/answers/238757-create-checkerboard-plot-from-a-matrix-where-each-point-is-represented-by-a-colored-square

data = readmatrix("../outputs/current.txt",'NumHeaderLines',0);
lda = (3*Nx);
Nt = size(data,1)/lda;
distribution = cell(Nt,1);

for i = 1:Nt
    distribution{i}=data(lda*(i-1)+(Nx)*(si-1)+1:lda*(i-1)+(Nx)*(si-1)+Nx, 1:Ny);
end

figure(1);
clf;

for i = 1:Nt
    C = distribution{i}';
    C = [[C zeros(size(C,1),1)] ; zeros(1,size(C,2)+1)];% Last row/col ignored
    pcolor(XGrid,YGrid,C)
    hold on;
    %pcolor(X,Y,distribution{i}');
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    xlabel("x");
    ylabel("y");
    colorbar;
    plot3(X,Y,distribution{i}','x','color','k');

    pause(dt);
end

