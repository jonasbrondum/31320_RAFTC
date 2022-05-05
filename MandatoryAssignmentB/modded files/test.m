close all
% W0 = (1-(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2))*10
% [A,B,C,D] = tf2ss(W0.Numerator{1},W0.Denominator{1});
% 
% wb=31.67*2;
% Q=10;
% 
% W2 = (W0*(1-(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)))*50
% [A,B,C,D] = tf2ss(W2.Numerator{1},W2.Denominator{1});
% 
% W2 = ss(A,B,C,D);
% 
% bodemag(W2)


wb=31.67;
Q=5;

W0 = (s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)

wb=65.8;

W0 = W0*(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)


[A,B,C,D] = tf2ss(W0.Numerator{1},W0.Denominator{1});

notches=ss(A,B,C,D)

bodemag((1-notches)*5)
W1 = notches*makeweight(2,10,0.01);


%Setup 3

%wb=31.67;
wb=35.6747;
Q=7;

W0 = (s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)

%wb=65.8;

wb=65.7936;

W0 = W0*(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)

%wb=25;

%W0 = W0*(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)

wb=25;

W0 = W0*(s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)


[A,B,C,D] = tf2ss(W0.Numerator{1},W0.Denominator{1});

W2=(1-ss(A,B,C,D))*5
%Setup 3
wb=25;

W0 = (s^2 + wb^2)/(s^2+(wb/Q)*s+wb^2)
M=100; wb=3; A=1.e-4; % Hvad betyder disse og hvor i bogen kommer de fra?

W1 = W0*tf([1/M wb], [1 wb*A]);


[A,B,C,D] = tf2ss(W1.Numerator{1},W1.Denominator{1});

W1 = ss(A,B,C,D);


