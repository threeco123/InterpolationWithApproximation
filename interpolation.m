clc;
factor=2;
%Read sample wav file
[original,Fs,nbits]=wavread('testing6.wav'); 
original=original(:,1);
%store into original array
%count number of elements in the original sample.wav
n=length(original);
Wave_sample_rate=Fs;
Wave_bits=nbits;
%duration of the original sample in seconds
duration=n/Fs;
%x is an evenly spaced from 0 to duration
%x has n elments
x=linspace(0,duration,n)'; 
%
%
%
%Interpolation with natural cubic spline 
%Use the pseudo code as shown in the txt book
a=zeros(n,1);
l=zeros(n,1);
u=zeros(n,1);
z=zeros(n,1);
c=zeros(n,1);
b=zeros(n,1);
d=zeros(n,1);
for i=1:n;
    a(i)=original(i);
end
%h is the constant step 
h=x(2)-x(1);
alpha=zeros(n-1,1);
for i=2:1:n-1
     alpha(i)=(3/h)*(a(i+1)-a(i))-(3/h)*(a(i)-a(i-1));
end
l(1)=1;
u(1)=0;
z(1)=0;
for i=2:1:n-1;
    l(i)=2*(x(i+1)-x(i-1))-h*u(i-1);
    u(i)=h/l(i);
    z(i)=(alpha(i)-h*z(i-1))/l(i);
end
l(n)=1;
z(n)=0;
c(n)=0;
for j=n-1:-1:1;
    c(j)=z(j)-u(j)*c(j+1);
    b(j)=(a(j+1)-a(j))/h-h*(c(j+1)+2*c(j))/3;
    d(j)=(c(j+1)-c(j))/(3*h);
end
s=zeros(n-1,1);
for j=1:n-1;
    k=x(j)+h/factor;
    while k<x(j+1);
        s(j)=a(j)+b(j)*(k-x(j))+c(j)*(k-x(j))^2+d(j)*(k-x(j))^3;
        k=k+h/factor;
    end
end
%
%Time-stretching
interpo_n=factor*n;
interpo_wave=zeros(interpo_n,1);
%x_interpo is an evenly spaced from 0 to new duration
%x has interpo_n elments
x_interpo=linspace(0,factor*duration,interpo_n)'; 
i=1;
j=1;
%interpolation wave
while i<=interpo_n;
    interpo_wave(i)=original(j);
    interpo_wave(i+1)=s(j);
    i=i+factor;
    j=j+1;
    if j>=n-1;
        break
    end
end
%
%Plot
subplot(2,1,1);
plot(x,original);
xlabel('Time(s)');
ylabel('Amplitude');
title('Original Signal before Time-stretching');
axis([0 nbits -1 1]);
subplot(2,1,2);
plot(x_interpo,interpo_wave);
xlabel('Time(s)');
ylabel('Amplitude');
title('Singal after Natural Cubic Spline Interpolation');
axis([0 nbits -1 1]);
%Write the interpo_wave to Interpolation.wav at original rate and numbers of bits
wavwrite(interpo_wave,Fs,nbits,'Interpolation.wav')
%play the interpolated file after time-stretching
soundsc(interpo_wave,Fs)