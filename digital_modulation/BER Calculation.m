%%bitstream generation and BPSK symbols
Bits=randi([0,1],1,200000);
Symb=Bits*2-1;
%%Noise generation
ni=randn(1,200000);
nj=randn(1,200000);
N=ni+j*nj;
V= [0.1, 0.1259, 0.1585, 0.1995, 0.2512, 0.3162, 0.3981, 0.5012, 0.6310, 0.7943, 1];%% variance values
%%simulate for BER CALC
ber=[];
for v=V
    n=N.*sqrt(v);
    signal=Symb+n;
    %%symbols detection and demodulation;
    reSymb=ones(1,200000);
    reSymb(find(real(signal)<=0))=-1;
    reBits=ones(1,200000);
    reBits(find(reSymb==-1))=0; 
    %%ber calc
    b=sum(Bits~=reBits)/200000;
    ber=[ber b];  
end
snr=10*log10(1./V); %% snr_dB
semilogy(snr,ber);

%%%%%%%%%%% QPSK %%%%%%%%%%%%%
%% generate bit stream
bitsNum=20000;
I = 2*randi([0,1],1,bitsNum)-1;
Q = 2*randi([0,1],1,bitsNum)-1;

%% QPSK Modulation
Syms=sqrt(0.5)*(I+j*Q); %% 0*2-1 = -1 AND 1*2-1 = 1 (THUS from 0/1 to -1/1)

%% AWGN
noise=randn(1,bitsNum)+i*randn(1,bitsNum);

%% BER Calculation
ber=[]; %empty vector for concatination ber reults
var=[0.1, 0.1259, 0.1585, 0.1995, 0.2512, 0.3162, 0.3981, 0.5012, 0.6310, 0.7943, 1]; 
for v=var;
    n=sqrt(v/2).*noise;
    
    % recived signal:
    y=Syms+n;
    
    % ML rule:
    deI=ones(1,bitsNum); % initiat detection vector with all ones.
    deI(find(real(y)<=0))=-1; % replace the value by -1 in the corresponding locations to real(y)<0
                                                      % imagiry part is ignored as it's pure noise                           
    
    deQ=ones(1,bitsNum); % initiat detection vector with all ones.
    deQ(find(imag(y)<=0))=-1;
    
    % Demodulation
      
    % BER
    b=(sum(I~=deI)+sum(Q~=deQ))/(2*bitsNum);
    ber=[ber b];                                                   
    
end
snr=10*log10(1./var);
hold on;
semilogy(snr,ber,'r');

