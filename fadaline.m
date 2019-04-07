%
clear all; close all; clc
fs=1e4; %Sampling frequency
t=0:1/fs:0.1-1/fs; %time vector
f_fundamental=50; %fundamental frequency in Hertz
y=zeros(size(t));
for k=1:2:12
    y=y+ (1/k)*sin(2*pi*f_fundamental*k*t);  %Input Signal that you can play with
end

%%
w_temel=2*pi*f_fundamental; %Fundamental Frequency 
M=20; % Number of harmonics that you want to detect
eta=0.01; %(learning rate)
W=zeros(2*M,length(t)); % Weights 
y_est=zeros(length(t),1); %neural network's output
er=y_est; %error vector
b=zeros(length(t),1);
for n=1:length(t)-1
    x=[];
    for k=1:M
    x=[x;sin(k*w_temel*t(n));cos(k*w_temel*t(n))]; %input
    end
y_est(n)=W(:,n)'*x +b(n); % output
er(n)=(y(n)-y_est(n));   % Error
W(:,n+1)= W(:,n)+eta*er(n)*x; %Update weights 
b(n+1)=b(n)+eta*er(n); %Update Bias
end

A=W(1:2:end,end); %A coefficients of Fourier Transform (sine)
B=W(2:2:end,end); %B coefficients of Fourier Transform (cosine)
G=sqrt(A.^2 +B.^2);
% G=G/max(G);
%% Plotting
figure
plot(t,W)% Change of Weights in time
title('Change of Weights')
saveas(gcf,'Weights.png')


figure
plot(t,b)
title('Change of Bias in time')
saveas(gcf,'Bias.png')


figure
plot(t,er)% 
title('Change of Error')
saveas(gcf,'Error.png')


figure
subplot(3,1,1)
plot(t,y)
title('Input Signal')
xlabel('Time(s)')

subplot(3,1,2)
bar(f_fundamental:f_fundamental:M*f_fundamental,G)
xlim([0 (M)*f_fundamental])
set(gca, 'XTick',(0:f_fundamental:M*f_fundamental))
title('Forier-ADALINE')

subplot(3,1,3)
bar(0:fs/length(t):(fs-fs/length(t)),2*abs(fft(y))/length(t),length(t)/f_fundamental/M*4)
xlim([0 (M)*f_fundamental])
set(gca, 'XTick',(0:f_fundamental:M*f_fundamental))
title('FFT')
xlabel('Frequency [Hz]')
saveas(gcf,'FFT_FADALINE.png')
