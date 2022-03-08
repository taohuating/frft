clear all;
close all;
clc;
N=1023;             %Sampling point
fi=0;                    %Initial frequency
f0=4e3;             %Initial frequency
fs=48e3;            %Sampling frequency
k=2e5;               %Frequency modulation rate  
f=zeros(1,N);
t=zeros(1,N);
tongbu_data=[1 1];      %Synchronous data
zero=zeros(N,1);
shuju_data=randsrc(1,1024,[1,0;0.5,0.5]);
%shuju_data=[1 0 0 1 1 1 0 0];
%%Interleaving, scrambling
M=size(shuju_data,2);
s=2^(ceil(log2(M)));
q=s/4-1;
M=length(shuju_data);
p=zeros(1,M);
p(1)=0;
for jj=2:M
    p(jj)=mod(21*p(jj-1)+q,s);
    while (p(jj)>=M)
        p(jj)=mod(21*p(jj)+q,s);
    end
end
data_interweave=zeros(1,M);    %Set space for interleaving vectors
for ii=1:M
    data_interweave(p(ii)+1)=shuju_data(ii);
end
figure(1);
plot(data_interweave,'*r-');
hold; plot(shuju_data,'k>-');
title('interweave');
%% modulation
%tongbu_data=[1 1];
%zero=zeros(T*fs+1,1);
%xunlian_data=[1,1,1];     %Training data
%shuju_data=randsrc(1,1024,[1,0;0.5,0.5]);    
%shuju_data=[1 0 0 1 1 1 0 0];
SS=length(shuju_data);         %Length of data sent
fasong_data=[shuju_data];      %send data
for n=-511:511
    tongbu(n+512)=cos(2*pi*f0*n/fs+2*pi*k*((n/fs)*(n/fs)));  %Generate signal t=n/fs
    t(n+512)=n/fs;      %time               
    f(n+512)=n*fs/N;  %frequency         
end
figure(2)
plot(t,tongbu);
%Add coordinates here
title('upLFM');
xlabel('t');
ylabel('Amplitude')     %amplitude
for n=-511:511
    upLFM(n+512)=cos(2*pi*f0*n/fs+pi*k*((n/fs)*(n/fs)));  %Generate signal t=n/fs
    %Compared with the synchronous signal, the frequency modulation rate is k, while the frequency modulation rate of the synchronous signal is 2K
    t(n+512)=n/fs;      %time               
    f(n+512)=n*fs/N;  %frequency          
end
figure(3)
plot(t,upLFM);
title('upLFM£¬Compared with the synchronous signal, the frequency modulation rate is half of it');
for n=-511:511
    downLFM(n+512)=cos(2*pi*f0*n/fs-pi*k*((n/fs)*(n/fs)));  %Generate signal t=n/fs
    t(n+512)=n/fs;      %time               
    f(n+512)=n*fs/N;    %frequency           
end
figure(4)
plot(t,downLFM);
title('downLFM');

tongbu_modulate=tongbu_data'*tongbu;        
%%Synchronous data modulation, transpose and multiply the two signals. After the two signals are transposed, they can be multiplied by 1023. Prime should mean transpose

up_modulate=fasong_data'*upLFM;                    %The transmitted data uses upchirp modulation 1023
down_modulate=(~fasong_data)'*downLFM;      %The transmitted data is modulated by downchirp, ~ which means the inverse 1023
modulate_matrix=up_modulate+down_modulate;  %When the corresponding bits are added, it is still 1023 elements

modulate_signal=reshape(modulate_matrix',[],1);   
tongbu_modulate1=reshape(tongbu_modulate',[],1);

signal=[tongbu_modulate1',zero',modulate_signal'];  %Synchronous signal - 0-modulated signal

signal_data=signal'
%figure;plot(signal_data);
hn=[1];     %channel
jieshou_data1=conv(signal_data,hn);      %Receive data, signal through channel, convolution
SNR=[-24:0];                             %Signal to noise ratio
for iii=1:length(SNR)
jieshou_data=awgn(jieshou_data1,SNR(iii));
N_dianshu=N;
%% relevant
%xg=xcorr(BC,tongbu_modulate1');
%figure;plot(xg);
%[y_max2,x_max2]=max(xg);
%qishi=x_max2-length(BC);
qishi=0;                      %start
jieshou_shuju=jieshou_data(qishi+3*N_dianshu+1:qishi+(SS+3)*N_dianshu);
%% Matched filtering
code_jietiao=[];
for kk=1:SS
   n_start=(kk-1)*N+1;
   n_end=kk*N;

   X=jieshou_shuju(n_start:n_end);
   
   XX1=xcorr(X,upLFM);
   MAX1=max(XX1);
   
   XX2=xcorr(X,downLFM);
   MAX2=max(XX2);
   
   MAX3=[MAX1 MAX2];
   
   MAX=max(MAX3);
   if(MAX1==MAX)
       code_jietiao=[code_jietiao,1];
   elseif(MAX2==MAX) code_jietiao=[code_jietiao,0];
   end
   MAX    
end
[number(iii),ratio(iii)]=biterr(shuju_data,code_jietiao);
end
ratio
save WML_BOK_pipei.txt ratio -ascii;
   

















