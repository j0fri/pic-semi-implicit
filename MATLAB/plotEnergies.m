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

fieldData = readmatrix("../outputs/fieldEnergy.txt",'NumHeaderLines',0);
speciesData = readmatrix("../outputs/speciesEnergy.txt",'NumHeaderLines',0);

PE = fieldData(:);
KE = zeros(size(PE));
for si = 1:Ns
    %disp(speciesData(2+2*(si-1):(2*Ns):end,1));
    KE = KE + speciesData(2+2*(si-1):(2*Ns):end,1);
end

Nt = size(PE,1);
t = linspace(0,T,Nt);

figure(1);
clf;
plot(t,KE);
hold on;
plot(t,PE);
plot(t,PE+KE);
legend('Kinetic Energy','Field Energy','Total Energy');
xlabel('Time');
ylabel('Energy');

figure(2);
clf;
semilogy(t,KE);
hold on;
semilogy(t,PE);
semilogy(t,KE+PE);
legend('Kinetic Energy','Field Energy','Total Energy');
xlabel('Time');
ylabel('Energy');


figure(3);
clf;
semilogy(t,PE);
hold on
legend('Field Energy','Linear damping');
xlabel('Time');
ylabel('Energy');

figure(4);
clf
semilogy(t,KE+PE);
hold on
title('Total energy evolution in time');
xlabel('Time')
ylabel('Energy')

disp('Max percentage difference in energy')
disp((max(KE+PE)-min(KE+PE))/(mean(KE+PE)))