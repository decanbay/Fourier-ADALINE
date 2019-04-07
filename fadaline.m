%%
% clear all; close all; clc
% fs=1e5; %örnekleme frekansý (Sampling frequency)
% t=0:1/fs:0.3-1/fs; %zaman vektörü
% ftemel=80; %Giriþ sinaylinin temel frekensý 1/s veya Hz
% y=zeros(size(t));
% for k=1:2:10
%     y=y+ (1/k)*sin(2*pi*ftemel*k*t);  %giriþ sinyali
% end


clear all; clc
fs=1e4;
t0=0:1/fs:(0.2-1/fs);
f0 = 49.5
input0 = 0.25+sin(2*pi*f0*t0+deg2rad(10))+ 0.5*sin(3*2*pi*f0*t0+deg2rad(-23)) + ...
    0.25*sin(5*2*pi*f0*t0+deg2rad(30)) +0.13*sin(7*2*pi*f0*t0+deg2rad(-15))...
    +0.06*sin(9*2*pi*f0*t0+deg2rad(50))+0.03*sin(11*2*pi*f0*t0+deg2rad(60)) ;
% input0 = awgn(input0,20,'measured');
t1 = 0.2:1/fs:(0.4-1/fs);
f1 = 50
input1 = 0.15+0.75*sin(2*pi*f1*t1+deg2rad(20))+ 0.1*sin(3*2*pi*f1*t1+deg2rad(-20)) + ...
    0.45*sin(5*2*pi*f1*t1+deg2rad(20))+0.64*sin(7*2*pi*f1*t1+deg2rad(-20))...
    +0.3*sin(9*2*pi*f1*t1+deg2rad(45))+0.25*sin(11*2*pi*f1*t1+deg2rad(65)) ;
% input1 = awgn(input1,10,'measured');

t2 = 0.4:1/fs:0.6;
f2 = 50.25
input2 = 0.35+sin(2*pi*f2*t2+deg2rad(30))+ 0.48*sin(3*2*pi*f2*t2+deg2rad(-23)) + ...
    0.05*sin(5*2*pi*f2*t2+deg2rad(30))+0.23*sin(7*2*pi*f2*t2+deg2rad(-15))...
    +0.13*sin(9*2*pi*f2*t2+deg2rad(35))+0.02*sin(11*2*pi*f2*t2+deg2rad(50)) ;
% input2 = awgn(input2,5,'measured');

y = [input0, input1, input2];
t=[t0,t1,t2];

ftemel=50

%  y=10*sin(2*pi*ftemel*t)+5*cos(2*pi*5*ftemel*t);
%%
w_temel=2*pi*ftemel; %Temel frekans rad/s
M=20; % belirlenmek istenen harmonik sayýsý
eta=0.01; %öðrenme oraný (learnin rate)
W=zeros(2*M,length(t)); % Aðýrlýk matrisi
y_est=zeros(length(t),1); %neural network un çýktýsý
er=y_est; %hata vektörü
b=zeros(length(t),1);
for n=1:length(t)-1
    x=[];
    for k=1:M
    x=[x;sin(k*w_temel*t(n));cos(k*w_temel*t(n))]; %input
    end
y_est(n)=W(:,n)'*x +b(n); % output
er(n)=(y(n)-y_est(n));   %hata
W(:,n+1)= W(:,n)+eta*er(n)*x; %aðýrlýklarýn güncellenmesi
b(n+1)=b(n)+eta*er(n); %bias ýn güncellenmesi
end

A=W(1:2:end,end); %A katsayýlarý
B=W(2:2:end,end); %B katsayýlarý
G=sqrt(A.^2 +B.^2);
% G=G/max(G);
%% Plotting
figure
plot(t,W)% aðýrlýklarýn zamanla deðiþimi
title('Aðýrlýklarýn Zamanla Deðiþimi')

figure
plot(t,b)% aðýrlýklarýn zamanla deðiþimi
title('Bias''ýn Zamanla Deðiþimi')

figure
plot(t,er)% aðýrlýklarýn zamanla deðiþimi
title('Hatanýn Zamanla Deðiþimi')

figure
subplot(3,1,1)
plot(t,y)
title('Input Signal')
xlabel('Zaman(s)')

subplot(3,1,2)
bar(ftemel:ftemel:M*ftemel,G)
xlim([0 (M)*ftemel])
set(gca, 'XTick',(0:ftemel:M*ftemel))
title('Fadaline Algoritmasý')

subplot(3,1,3)
bar(0:fs/length(t):(fs-fs/length(t)),2*abs(fft(y))/length(t),length(t)/ftemel/M)
xlim([0 (M)*ftemel])
set(gca, 'XTick',(0:ftemel:M*ftemel))
title('FFT')
xlabel('Frekans [Hz]')