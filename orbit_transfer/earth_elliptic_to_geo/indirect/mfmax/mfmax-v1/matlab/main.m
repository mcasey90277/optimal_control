
clear all;

data;

nb_caches = 5;

nntwarn off;

P = data(:,1:2)';
T = data(:,10)';

[w1,b1,w2,b2] = initff([0.1 10.; 1.1 3.5],nb_caches,'tansig',1,'logsig');

[w1,b1,w2,b2,ep,tr] = trainbp(w1,b1,'tansig',w2,b2,'logsig',P,T);

simuff(P,w1,b1,'tansig',w2,b2,'logsig')
