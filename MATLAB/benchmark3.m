clear;
semiimplicitDirs = dir('../outputs/benchmark3/semiimplicit/');
semiimplicitNts = [];
semiimplicitEnergies = [];
semiimplicitRuntimes = [];
for i = 3:length(semiimplicitDirs)
    Nti = uint64(str2double(semiimplicitDirs(i).name));
    semiimplicitNts = [semiimplicitNts;Nti];
    data = readmatrix(['../outputs/benchmark3/semiimplicit/',num2str(Nti),'/runtime.txt']);
    semiimplicitRuntimes = [semiimplicitRuntimes;data(2)/10];
    fieldData = readmatrix(['../outputs/benchmark3/semiimplicit/',num2str(Nti),'/fieldEnergy.txt'],'NumHeaderLines',0);
    speciesData = readmatrix(['../outputs/benchmark3/semiimplicit/',num2str(Nti),'/speciesEnergy.txt'],'NumHeaderLines',0);
    
    PE = fieldData(:);
    KE = zeros(size(PE));
    for si = 1:2
        %disp(speciesData(2+2*(si-1):(2*Ns):end,1));
        KE = KE + speciesData(2+2*(si-1):(2*2):end,1);
    end
    deltaE = calcDeltaE(PE+KE);
    semiimplicitEnergies = [semiimplicitEnergies; deltaE];
end
[semiimplicitNts,sortIdx] = sort(semiimplicitNts,'ascend');
semiimplicitRuntimes = semiimplicitRuntimes(sortIdx);
semiimplicitEnergies = semiimplicitEnergies(sortIdx);


explicitDirs = dir('../outputs/benchmark3/explicit/');
explicitNts = [];
explicitEnergies = [];
explicitRuntimes = [];
for i = 3:length(explicitDirs)
    Nti = uint64(str2double(explicitDirs(i).name));
    explicitNts = [explicitNts;Nti];
    data = readmatrix(['../outputs/benchmark3/explicit/',num2str(Nti),'/runtime.txt']);
    explicitRuntimes = [explicitRuntimes;data(2)/10];
    fieldData = readmatrix(['../outputs/benchmark3/explicit/',num2str(Nti),'/fieldEnergy.txt'],'NumHeaderLines',0);
    speciesData = readmatrix(['../outputs/benchmark3/explicit/',num2str(Nti),'/speciesEnergy.txt'],'NumHeaderLines',0);
    
    PE = fieldData(:);
    KE = zeros(size(PE));
    for si = 1:2
        %disp(speciesData(2+2*(si-1):(2*Ns):end,1));
        KE = KE + speciesData(2+2*(si-1):(2*2):end,1);
    end
    deltaE = calcDeltaE(PE+KE);
    explicitEnergies = [explicitEnergies; deltaE];
end
[explicitNts,sortIdx] = sort(explicitNts,'ascend');
explicitRuntimes = explicitRuntimes(sortIdx);
explicitEnergies = explicitEnergies(sortIdx);


figure('DefaultAxesFontSize',22);
clf;
loglog(double(semiimplicitNts),semiimplicitEnergies,'-','color','k');
hold on;
%loglog(semiimplicitNts*10,double(semiimplicitNts).^2*semiimplicitRuntimes(end)/double(semiimplicitNts(end))^2,'--','color','k')
loglog(double(explicitNts),explicitEnergies,'--','color','k');
%loglog(explicitNts*10,double(explicitNts)*explicitRuntimes(end)/double(explicitNts(end)),'-.','color','k')
grid on;
axis square;
xlabel('Time steps');
ylabel('Energy Change')
legend('Semi-implicit scheme', 'Explicit scheme');

figure('DefaultAxesFontSize',22);
clf;
loglog(semiimplicitRuntimes,semiimplicitEnergies,'-','color','k');
hold on;
%loglog(semiimplicitNts*10,double(semiimplicitNts).^2*semiimplicitRuntimes(end)/double(semiimplicitNts(end))^2,'--','color','k')
loglog(explicitRuntimes,explicitEnergies,'--','color','k');
%loglog(explicitNts*10,double(explicitNts)*explicitRuntimes(end)/double(explicitNts(end)),'-.','color','k')
grid on;
axis square;
xlabel('Runtime');
ylabel('Energy Change')
legend('Semi-implicit scheme', 'Explicit scheme');

function deltaE = calcDeltaE(E)
    deltaE = abs(E(1)-E(end));
end