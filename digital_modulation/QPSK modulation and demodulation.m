clc;
clear all;
close all;
%% Bit stream
L=20;
stream_real=randi([0,1],1,L/2);
stream_imag=randi([0,1],1,L/2);
output= myQPSK(stream_real,stream_imag);
%% carrier
Ts=100;
fc=1;
t=linspace(0,L/2,L/2*Ts);
real_carrier=cos(2*pi*fc.*t);
imag_carrier=sin(2*pi*fc.*t);
%% QPSK modulation
signal= reshape(repmat(output,100,1),1,[]);
mod_sig= abs(signal).*(real(signal).*real_carrier+ imag(signal).*imag_carrier);
figure
subplot(5,1,1)
stem(stream_real)
subplot(5,1,2)
stem(stream_imag)
subplot(5,1,3)
plot(t,mod_sig)

%% QPSK demodulation
% demod_real= real(mod_sig).*real_carrier;
% demod_imag= imag(mod_sig).*imag_carrier;
demod=mod_sig.*real_carrier+1i.*mod_sig.*imag_carrier;
detected=[ceil(intdump(real(demod),100));ceil(intdump(imag(demod),100))];
real_ber= sum(detected(1,:)== stream_real)
imag_ber= sum(detected(2,:)== stream_imag)
subplot(5,1,4)
stem(detected(1,:))
subplot(5,1,5)
stem(detected(2,:))


