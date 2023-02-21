Nt = 6;

data = cell(Nt,1);
for i = 1:Nt
    data{i} = readmatrix(['../outputs/presetDistributionTest',num2str(i)]);
end

figure(1);
clf

[p,x] = hist(data{1}(:,1),50);
int = trapz(x,p);
plot(x,p/int,'.-');
hold on;
[p,x] = hist(data{2}(:,1),50);
int = trapz(x,p);
plot(x,p/int,'*-');

a = 5; %m/2KbT0
plot(-1:0.01:1, sqrt(a/pi)*exp(-a*(-1:0.01:1).^2));
title('1D boltzmann comparison')
legend('Np=1000','Np=10000','Theoretical');
xlabel('velocity')
ylabel('density function');

figure(2);
clf
vel3rms1 = sqrt(data{3}(:,1).*data{3}(:,1)+data{3}(:,2).*data{3}(:,2)+data{3}(:,3).*data{3}(:,3));
[p,x] = hist(vel3rms1,50);
int = trapz(x,p);
plot(x,p/int,'.-');
hold on;
vel3rms2 = sqrt(data{4}(:,1).*data{4}(:,1)+data{4}(:,2).*data{4}(:,2)+data{4}(:,3).*data{4}(:,3));
[p,x] = hist(vel3rms2,50);
int = trapz(x,p);
plot(x,p/int,'*-');

a = 5; %m/2KbT0
plot(0:0.01:1, (a/pi)^1.5*4*pi*(0:0.01:1).^2.*exp(-a*(0:0.01:1).^2));
title('3D boltzmann comparison')
legend('Np=1000','Np=10000','Theoretical');
xlabel('RMS Speed')
ylabel('density function');

figure(3);
clf
hist3(data{6}(:,1:2),'CdataMode','auto');
colormap('hot');
colorbar;
xlabel('x');
ylabel('y');
title('Contour plot of density function')