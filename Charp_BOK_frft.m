clear all;
close all;
clc;
N=511;         %Sampling points
fi=0;              %Initial phase
f0=4e3;        %Initial frequency
fs=24e3;       %Sampling frequency
k=2e5;
%k=-2e5;
%k=0;
%%
for n=-255:255
    x(n+256)=exp(j*2*pi*fi+j*2*pi*f0*n/fs+j*pi*k*((n/fs)*(n/fs)));    %generate signal t=n/fs
    t(n+256)=n/fs;        %time 
    f(n+256)=n*fs/N;      %freq
end
%% 
figure(1);
plot(t,x);
title('k= - 2e5');

%%
tt=t(N)-t(1);            %Uniformly set signal interception time
kk=tt/fs;                 %Transformation of time scale
pg=acot(-k*kk)*2/pi;     %k=-cot(p*pi/2)/kk£¬   kk,It's a time scale£¬pg Is the best order of calculation
%%   
X=awgn(x,20);                      %Add white noise to signal x (n)£¬SNR=20
r=0.01;
a=[-4:r:4];                              %fractional power£¬reversal 0~2
G=zeros(length(a),length(X));      %for different power do frft
H=12;                              %Modular cache variable
ph=zeros(1,4);
nm=1;
for l=1:length(a)                %Cycle 801 times to search for peak value
    F=frft(X,a(l));
    G(l,:)=abs(F(:));
  if(H<=max(abs(F(:))))
    ph(nm)=a(l);                  %The order corresponding to the maximum modulus  a(l)£¬pl
    nm=nm+1;
  end
end
for m=1:length(a);                %Search space of order p
    y(m)=max(abs(G(m,:)));
end
figure(2);
stem(a,y);                         %stem You can look at the curve projected on the order axis
title('fo=4e3, fs=24e3, k=2e5');
grid on

%%
figure(3);
mesh(G);                           %mesh
title('fo=4e3, fs=24e3, k= 2e5');
grid on
%%
F=frft(X,ph(3));               %Call frft()
[aa,bb]=max(abs(F));
 figure(4);
 plot(t,F);
 grid on
 title('k= 2e5 p= 0.89 FrFT Energy Focus ');
%% 
Fg=frft(X,pg);    %Fourier transform the signal according to the order PG measured by frequency modulation frequency, and observe the focusing condition
[aag,bbg]=max(abs(Fg));
 figure(5);
 plot(t,Fg);
 title('k= 2e5 p= 0.8884 FrFT Energy Focus ');
grid on
  
%% 
if 1
%%The above code explores the feasibility of this detection method
%%The following is the simulation code
f=zeros(1,N);               %Generate a frequency vector space with n zeros
t=zeros(1,N);               
tongbu_data=[1 1];          %Synchronization signal
fasong_data=randsrc(1,1024,[1,0;0.5,0.5]);       %Generate 500 1 or 0 random transmission data sequences
%fasong_data=[1 0 1];
SS=length(fasong_data);      %SS  Length of data sent
zero=zeros(N,1);                    % 0 vector of N rows and 1 columns

%% Interleaving, scrambling
M=size(fasong_data,2);       %N-th power of 2      
s=2^(ceil(log2(M)));              %Take the nearest integer ceil
q=s/4-1;
M=length(fasong_data);
p=zeros(1,M);
p(1)=0;
for jj=2:M
    p(jj)=mod(21*p(jj-1)+q,s);
    while (p(jj)>=M)
        p(jj)=mod(21*p(jj)+q,s);
    end
end
data_interweave=zeros(1,M);   %interwaeav
for ii=1:M
    data_interweave(p(ii)+1)=fasong_data(ii);
end
figure(6);
plot(fasong_data,'*r-');                         %Interleaved data
hold; plot(data_interweave,'kO-');        %send data
title('Comparison of original data and interleaved data')

%% modulation
hn=1;
for n=-255:255
   %tongbu(n+256)=cos(2*pi*f0*n/fs+2*pi*k*((n/fs)*(n/fs))); 
   %The generated signal T = n / Fs, and the synchronous signal uses twice the frequency modulation rate for fractional Fourier transform
   tongbu(n+256)=exp(j*2*pi*fi+j*2*pi*f0*n/fs+j*2*pi*k*((n/fs)*(n/fs)));       
    t(n+256)=n/fs;          %time 
    f(n+256)=n*fs/N;      %freq
end
tt=t(N)-t(1);   %Uniformly set signal interception time
kk=tt/fs;        %Transformation of time scale
p_tb=acot(-2*k*kk)*2/pi;      %Fractional order of predictive synchronization signal

%%
r=0.01;
a=[-4:r:4];                               %fractional power£¬reversal 0~2
G=zeros(length(a),length(tongbu));      %for different power do frft
H=12;                                   %Modular cache variable
ph=zeros(1,4);
nm=1;
for l=1:length(a)                       %Cycle 100 times
     T=frft(tongbu,a(l));
     G(l,:)=abs(T(:));
  if(H<=max(abs(T(:))))
     ph(nm)=a(l);                       %The order corresponding to the maximum modulus a(l)£¬pl
     nm=nm+1;
  end
end
for m=1:length(a);                 %Search space of order p
    y(m)=max(abs(G(m,:)));
end
figure(7);
stem(a,y);                         %stem£¬You can look at the curve projected on the order axis 
title('k=2e5 tongbu shousuo');

%%
F_tb=frft(tongbu,p_tb);
[aa_tb,bb_tb]=max(abs(F_tb));
 figure(8);
 plot(t,F_tb);
 title('The synchronization signal can be most focused......');

%%     
%k=data;
for n=-255:255
    up_LFM(n+256)=exp(j*2*pi*fi+j*2*pi*f0*n/fs+j*pi*k*((n/fs)*(n/fs)));      %Generate signal t=n/fs
    t(n+256)=n/fs;           %time 
    f(n+256)=n*fs/N;      %freq  
end
%figure
%plot(t,up_LFM);
%title('up_LFM');
p_up=acot(-k*kk)*2/pi;              %%Predicted chirp_ Fractional order p of up signal_ up
F_up=frft(up_LFM,p_up);  
   figure(9);
   plot(t,F_up);
   title('up_LFMEnergy focusing effect of signal ......');
for n=-255:255
    down_LFM(n+256)=exp(j*2*pi*fi+j*2*pi*f0*n/fs-j*pi*k*((n/fs)*(n/fs)));  %Generate signal t=n/fs
    %down_LFM(n+256)=cos(2*pi*f0*n/fs-pi*k*((n/fs)*(n/fs)));  %Generate signal t=n/fs  
     t(n+256)=n/fs;        %time 
    f(n+256)=n*fs/N;      %freq
end
%figure
%plot(down_LFM);
%title('down_LFM');
%%Prediction fractional order p_down
p_down=acot(k*kk)*2/pi;
F_down=frft(down_LFM,p_down);  
  figure(10);
  plot(t,F_down);
  title('down_LFM Energy focusing effect of signal...');
  grid on
%tongbu_modulate=tongbu_data'*tongbu;                     
%% Operators' and. ' It's the same. The matrix transposes a '* B. A transposes and then multiplies B
tongbu_matrix=tongbu_data'*tongbu;
up_modulate=fasong_data'*up_LFM;
down_modulate=(~fasong_data)'*down_LFM;
modulate_matrix=up_modulate+down_modulate;    %After the alignment positions are added, they are correct

for n=-255:255
   t(n+256)=n*(1/fs);
end
%%     
%modulate_tongbu=reshape(tongbu_matrix',[],1);
%modulate_signal=reshape(modulate_matrix',[],1);
%modulate_signal1=[modulate_signal(1:961)',zero',modulate_signal(962:9*961)'];
%signal_tongbu=modulate_tongbu';  
%signal_data=modulate_signal';
signal_tongbu=tongbu_matrix;  %Signal = synchronous signal, plus 0, plus signal
signal_data=modulate_matrix;
%jieshou_data1=conv(signal_data,hn);   
%The received data is equal to the convolution of the signal and the channel HN
SNR=-24:0;   %Transmit in different noise environments in order to count the bit error rate
for iii=1:length(SNR)           %lth= 16
jieshou_tongbu=awgn(signal_tongbu,SNR(iii));
jieshou_data=awgn(signal_data,SNR(iii));            %Add white noise to the signal
 
 %% demodulation
 N_jietiao=N;      %Points, sampling points per symbol
 h1=0;
  for n=-255:255
      t(n+256)=n*(1/fs);
  end
 code_jietiao=[];             %Put the demodulated code here
 for re=1:SS                  %SS=Length of transmitted signal
   %rd_start=(rd-1)*N_jietiao+1; 
   %rd_end=rd*N_jietiao;             %Rd symbol£¬rd_end-rd_start=N_jietiao=511
   X=jieshou_data(re,:);
   up_FF=frft(X,p_up);
   [aU,bU]=max(abs(up_FF));
 
   down_FF=frft(X,p_down);
   [aD,bD]=max(abs(down_FF));
 
   MaxAll=[aU aD];
   MAX=max(MaxAll);
   if(aU==MAX)
       code_jietiao(re)=1;
   elseif(aD==MAX) 
       code_jietiao(re)=0;
   end
   MAX
 end
 [number(iii),ratio(iii)]=biterr(fasong_data,code_jietiao);
end 
ratio
save WML_BOK_frft.txt ratio -ascii;
end

















