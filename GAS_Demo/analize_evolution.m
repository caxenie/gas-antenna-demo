% simple script to plot the candidate function for optimization using GAS
figure(1);
title('The best maximum across training - the fittest chromosome value');
load('gas_demo_log.txt');
plot(gas_demo_log(:,1), gas_demo_log(:,2));