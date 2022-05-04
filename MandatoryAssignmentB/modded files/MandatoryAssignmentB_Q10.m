%% 31320 RAFTC
% Question 10

MandatoryAssignmentB_Q9

%% Unfinished..

%% Result validationvalidate: 
figure ,bode(Pcl('y','r'), (G*K)/(1+G*K),'--r',{0.1,1000})
disp ('validate that G*K/(1+G*K) == transfer function from ref to y')
figure ,bode(Pcl('z1','r'), W1/(1+G*K),'--r',{0.1,1000})
disp ('validate that W1/(1+G*K) == transfer function from ref to z1')

% This closed loop validation is not right:
figure ,bode(Pcl('z2','r'), W2/(1+G*K),'--r',{0.1,1000})
disp ('validate that W2/(1+G*K) == transfer function from ref to z2')
figure ,bode(Pcl('e','r'), 1/(1+G*K),'--r',{0.1,1000})
disp ('validate that 1/(1+G*K) == transfer function from ref to e')

% Controller assumes positive feedback

hinfnorm(Pcl('y','r'))
hinfnorm(Pcl('z1','r'))
hinfnorm(Pcl('z2','r'))
hinfnorm(Pcl('e','r'))

figure;
step(Pcl)


%% Poles and zeros of controller:
pzG = pzmap(G);

pzK = pzmap(K);

pzW1 = pzmap(W1);

pzW2 = pzmap(W2);
