% Clear all variables and close all figures
clear all;
close all;

% Generate a random bit stream
L=10; %length of stream 10samples
BITS=randi([0,1],1,L);
ts=1; %sympoltime =1
t=linspace(0,L,1000); % generates  points between X1 and X2. no of samples =1000/10


%% ASK Modulation
fc=4; %carrier frequancy 
carrier_mes=sin(2*pi*fc*t);
mBits=reshape(repmat(BITS,100,1),1,1000)+1;
modSignal=mBits.*carrier_mes;
subplot(4,1,1);
stem(BITS);
subplot(4,1,2);
plot(t,carrier_mes);
subplot(4,1,3);
plot(t,modSignal);

%% MSK Modulation
fs=1;
fm=6; % any two frequencies with spacing (5)using buad rate=10 so fs=0.5fm 10/2=5
t=linspace(0,10,1000);
mBits=reshape(repmat(BITS,100,1),1,1000)*5+1;
modSignal=sin(2*pi.*mBits.*t);
%figure;
subplot(4,1,4)
plot(t,modSignal);
