clear
clc

y = dlmread('f14.txt');
dt = 0.01;
T = 5;
t = 0:dt:T;

figure(1,"position",get(0,"screensize"))
plot(t,y),grid;
%FFT : Standard
f = fft(y);

n = length(y);
c1(1:n) = 0;
c2(1:n) = 0;
c(1:n) = 0;

%Written FFT
for k = 1:n
    cos_sum = 0;
    sin_sum = 0;
    for  m = 1:n
        cos_sum = cos_sum + y(m)*cos(2*pi*(k-1)*(m-1)/n);
        sin_sum = sin_sum + y(m)*sin(2*pi*(k-1)*(m-1)/n); 
    end
    c1(k) = cos_sum/n;
    c2(k) = sin_sum/n;
    c(k) = sqrt(c1(k)*c1(k)+c2(k)*c2(k));
end

% plt1 : Written
figure(2,"position",get(0,"screensize"))
plot(t,c),grid;

% plt2 : Standart
figure(3,"position",get(0,"screensize"))
plot(abs(f)),grid;

k=0;

for i=3:(n/2)
    if((c(i-1)<c(i))&(c(i+1)<c(i))&&(c(i-2)<c(i))&(c(i+2)<c(i))) 
            k = k+1;
            max(k) = i;
    end
end
w = (max-1)/T;

% Display w
w
% Display max
max
% Display freq
c(w)


b(1) = 0;
for i= 1: n
    b(1) = b(1) + y(i)*t(i)*t(i)*t(i);
end

b(2) = 0;
for i= 1: n
    b(2) = b(2) + y(i)*t(i)*t(i);
end

b(3) = 0;
for i= 1: n
    b(3) = b(3) + y(i)*t(i);
end

b(4) = 0;
for i= 1: n
    b(4) = b(4) + y(i)*sin(2*pi*w(1)*t(i));
end

b(5) = 0;
for i= 1: n
    b(5) = b(5) + y(i);
end

% Display b
b

powT(1:6)=0;
for k=1:6
    powT(k)=sum(t.^k);
end

powS(1:3)=0;
for k=1:3
    powS(k)=sum((t.^k).*sin(2*pi*w(1).*t));
end

A(1:5,1:5)=[powT(6), powT(5), powT(4), powS(3), powT(3);
            powT(5), powT(4), powT(3), powS(2), powT(2);
            powT(4), powT(3), powT(2), powS(1), powT(1);
            powS(3), powS(2), powS(1), sum(sin(2*pi*w(1).*t).^2), sum(sin(2*pi*w(1).*t));
            powT(3), powT(2), powT(1), sum(sin(2*pi*w(1).*t)), 501];

% Display A
A

x=inv(A)*b';

approximation = x(1).*t.^3+x(2).*t.^2+x(3).*t+x(4)*sin(2*pi*20.*t)+x(5);
figure(4,"position",get(0,"screensize"))
plot(t,approximation),grid;

figure(4,"position",get(0,"screensize"))
plot(t,approximation),grid;
hold on
plot(t,y);

err = 0;
for i = 1:n 
    err = err + (approximation(i)-y(i))*(approximation(i)-y(i));
end;

% Display x
x

%Display error
err