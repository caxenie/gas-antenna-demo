% simple script to plot the candidate function for optimization using GAS
LOW = -1;
HI = 2;
% pointer to candidate function
fct=@cand_func;
fplot(fct, [LOW HI]);
figure(2);
title('The best maximum across training - the fittest chromosome value');
load('gas_demo_log.txt');
plot(gas_demo_log(:,1), gas_demo_log(:,2));