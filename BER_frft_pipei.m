clear all;
close all;
clc;
%SNR=[-14:2:1];
SNR=-24:0;
PIPEI=[4.7460938e-01   4.2089844e-01   4.1210938e-01   3.9355469e-01   3.5449219e-01 ...,
    3.1152344e-01   2.7246094e-01   2.1093750e-01   1.2207031e-01   7.8125000e-02 ...,
    3.0273438e-02   1.1718750e-02   2.9296875e-03   9.7656250e-04   0.0000000e+00 ...,
    0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00 ...,
    0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00 
 ];

FRFT=[4.7949219e-01   4.5898438e-01   4.4238281e-01   4.0917969e-01   3.5644531e-01 ...,   
2.9589844e-01   2.1582031e-01   1.2988281e-01   6.2500000e-02   2.0507813e-02 ...,   
9.7656250e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00 ...,   
0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00 ...,   
0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00    
];
figure
semilogy(SNR,FRFT,'k-*');
axis([-24,-8,0,0.6])
hold on;
semilogy(SNR,PIPEI,'r-o');
grid on;
title('Bit Error Rate Curve','fontname','Song style','fontsize',11)
xlabel('SNR /(dB)','fontname','Song style','fontsize',11)
ylabel('SBR','fontname','Song style','fontsize',11)
legend('FrFT Demodulation','Matching Filter Demodulation','fontname','Song style','fontsize',11)