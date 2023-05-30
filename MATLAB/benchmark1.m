semiimplicitDirs = dir('../outputs/benchmark1/semiimplicit/');
semiimplicitNxs = [];
semiimplicitRuntimes = [];
semiimplicitInitialisationTimes = [];
for i = 3:length(semiimplicitDirs)
    Nxi = uint64(str2double(semiimplicitDirs(i).name));
    semiimplicitNxs = [semiimplicitNxs;Nxi];
    data = readmatrix(['../outputs/benchmark1/semiimplicit/',num2str(Nxi),'/runtime.txt']);
    semiimplicitRuntimes = [semiimplicitRuntimes;data(2)/10];
    semiimplicitInitialisationTimes = [semiimplicitInitialisationTimes;data(1)];
end
[semiimplicitNxs,sortIdx] = sort(semiimplicitNxs,'ascend');
semiimplicitRuntimes = semiimplicitRuntimes(sortIdx);
semiimplicitInitialisationTimes = semiimplicitInitialisationTimes(sortIdx);


explicitDirs = dir('../outputs/benchmark1/explicit/');
explicitNxs = [];
explicitRuntimes = [];
explicitInitialisationTimes = [];
for i = 3:length(explicitDirs)
    Nxi = uint64(str2double(explicitDirs(i).name));
    explicitNxs = [explicitNxs;Nxi];
    data = readmatrix(['../outputs/benchmark1/explicit/',num2str(Nxi),'/runtime.txt']);
    explicitRuntimes = [explicitRuntimes;data(2)/10];
    explicitInitialisationTimes = [explicitInitialisationTimes;data(1)];
end
[explicitNxs,sortIdx] = sort(explicitNxs,'ascend');
explicitRuntimes = explicitRuntimes(sortIdx);
explicitInitialisationTimes = explicitInitialisationTimes(sortIdx);

figure('DefaultAxesFontSize',22);
clf;
loglog(semiimplicitNxs*10,semiimplicitRuntimes,'-','color','k');
hold on;
loglog(semiimplicitNxs*10,double(semiimplicitNxs).^2*semiimplicitRuntimes(end)/double(semiimplicitNxs(end))^2,'--','color','k')
loglog(explicitNxs*10,explicitRuntimes,'.','color','k');
loglog(explicitNxs*10,double(explicitNxs)*explicitRuntimes(end)/double(explicitNxs(end)),'-.','color','k')
grid on;
axis square;
xlabel('N{g}');
ylabel('Runtime per step (s)')
legend('Semi-implicit scheme', 'Explicit scheme');