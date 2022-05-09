%close all
mdl = 'coloumb_bandwitdh';
load_system(mdl);
io(1) = linio('coloumb_bandwitdh/Gain1',1,'input');
io(2) = linio('coloumb_bandwitdh/Gain',1,'output');
linsys1 = linearize(mdl,io);
sys = tf(linsys1)
bode(sys)


%The eigenvalue of the coloumb friction is 0.3993, so it should be fast
%enough to handle the coloumb friction.
