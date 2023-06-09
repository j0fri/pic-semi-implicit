load('rooflineData.mat');
figure('DefaultAxesFontSize',22);
clf;
loglog(rooflineData(:,1),1./rooflineData(:,2),'-x','color','k');
hold on
loglog(rooflineData(:,1),rooflineData(:,1)./rooflineData(1,2),'--','color','k');
legend('Scaling','Theoretical scaling')
grid on;
xlabel('Processes');
ylabel('1/Runtime');


figure('DefaultAxesFontSize',22);
clf;
semilogx(rooflineData(:,1),rooflineData(:,3),'-x','color','k');
hold on
semilogx(rooflineData(:,1),log(rooflineData(:,1))*log(rooflineData(3,2))/log(2),'--','color','k');
legend('Scaling','O(log(np))')
grid on;
xlabel('Processes');
ylabel('1/Runtime');