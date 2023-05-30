semiimplicitDirs = dir('../outputs/benchmark2/semiimplicit/');
semiimplicitNps = [];
semiimplicitRuntimes = [];
semiimplicitInitialisationTimes = [];
for i = 3:length(semiimplicitDirs)
    Npi = uint64(str2double(semiimplicitDirs(i).name));
    semiimplicitNps = [semiimplicitNps;Npi];
    data = readmatrix(['../outputs/benchmark2/semiimplicit/',num2str(Npi),'/runtime.txt']);
    semiimplicitRuntimes = [semiimplicitRuntimes;data(2)/10];
    semiimplicitInitialisationTimes = [semiimplicitInitialisationTimes;data(1)];
end
[semiimplicitNps,sortIdx] = sort(semiimplicitNps,'ascend');
semiimplicitRuntimes = semiimplicitRuntimes(sortIdx);
semiimplicitInitialisationTimes = semiimplicitInitialisationTimes(sortIdx);


explicitDirs = dir('../outputs/benchmark2/explicit/');
explicitNps = [];
explicitRuntimes = [];
explicitInitialisationTimes = [];
for i = 3:length(explicitDirs)
    Npi = uint64(str2double(explicitDirs(i).name));
    explicitNps = [explicitNps;Npi];
    data = readmatrix(['../outputs/benchmark2/explicit/',num2str(Npi),'/runtime.txt']);
    explicitRuntimes = [explicitRuntimes;data(2)/10];
    explicitInitialisationTimes = [explicitInitialisationTimes;data(1)];
end
[explicitNps,sortIdx] = sort(explicitNps,'ascend');
explicitRuntimes = explicitRuntimes(sortIdx);
explicitInitialisationTimes = explicitInitialisationTimes(sortIdx);

figure('DefaultAxesFontSize',22);
clf;
loglog(semiimplicitNps*10,semiimplicitRuntimes,'-','color','k');
hold on;
%loglog(semiimplicitNps*10,double(semiimplicitNps)*semiimplicitRuntimes(end)/double(semiimplicitNps(end)),'--','color','k')
loglog(explicitNps*10,explicitRuntimes,'--','color','k');
loglog(explicitNps*10,double(explicitNps)*explicitRuntimes(end)/double(explicitNps(end)),'-.','color','k')
grid on;
axis square;
xlabel('N{p}');
ylabel('Runtime per step (s)')
legend('Semi-implicit scheme', 'Explicit scheme', 'O(Np)');