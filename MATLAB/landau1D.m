fieldEnergy = readmatrix('../outputs/fieldEnergy.txt');
Nt = length(fieldEnergy);
T = 10;
dt = T/Nt;

speciesEnergy = readmatrix('../outputs/speciesEnergy.txt');
speciesEnergy = speciesEnergy(1:2:end) + speciesEnergy(2:2:end);

times = linspace(0,10,Nt);

figure(1);
plot(times, log(fieldEnergy));

figure(2);
plot(times,fieldEnergy);
hold on 
plot(times,speciesEnergy);
plot(times,fieldEnergy+speciesEnergy);
