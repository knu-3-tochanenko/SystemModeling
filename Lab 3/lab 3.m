clear
clc

yExperimental = dlmread('y4.txt',' ');
plot(0:0.2:50, yExperimental(1,:));
[M,N] = size(yExperimental);

c(1) = 0.14;
c(2) = 0.3;
c(3) = 0.2;
c(4) = 0.1;

m(1) = 13;
m(2) = 28;
m(3) = 23;

h = 0.2;
t = 0.0;
tk = 50.0;
eps = 1e-6;
I = inf;
numIter = 0;
yGained = zeros(M,N);

function res = getA(m, c)
    A = zeros(6,6);
    A(1, 2) = 1;
    A(2, 1) = -(c(1) + c(2)) / m(1);
    A(2, 3) = c(2) / m(1);
    A(3, 4) = 1;
    A(4, 1) = c(2) / m(2);
    A(4, 3) = -(c(2) + c(3)) / m(2);
    A(4, 5) = c(3) / m(2);
    A(5, 6) = 1;
    A(6, 3) = c(3) / m(3);
    A(6, 5) = -(c(3) + c(4)) / m(3);
    res = A;
end

function res = fy(A, y)
    res = A*y;
end

function res = fU(A, U, y, m, c)
    dAy = zeros(6,3);
    dAy(6, 1) = -(y(5) / m(3));
    dAy(2, 2) = -(c(2) * y(3) - (c(1) + c(2)) * y(1)) / (m(1) * m(1));
    dAy(6, 3) = -(c(3) * y(3) - (c(3) + c(4)) * y(5)) / (m(3) * m(3));
    res = A * U + dAy;
end

function res = U_RungeKutt(A, U, h, y, m, c)
    k1 = h * fU(A, U, y, m, c);
    k2 = h * fU(A, U + k1 / 2.0, y, m, c);
    k3 = h * fU(A, U + k2 / 2.0, y, m, c);
    k4 = h * fU(A, U + k3, y, m, c);
    res = U + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
endfunction

function res = Y_RungeKutt(A, y, h, m, c)
    k1 = h * fy(A, y, m, c);
    k2 = h * fy(A, y + k1 / 2.0, m, c);
    k3 = h * fy(A, y + k2 / 2.0, m, c);
    k4 = h * fy(A, y + k3, m, c);
    res = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
endfunction


while (I > eps)
    numIter = numIter+1;
    y = yExperimental(:,1);
    dy = zeros(M,1);
    U = zeros(6,3);
    B = zeros(3,3);
    b = zeros(3,1);
    I = 0.0;
    i = 2;

    for iter = 1:M
        yGained(iter,1) = y(iter,1);
    end

    while(i <= N)
        A = getA(m,c);
        Unew = U_RungeKutt(A, U, h, y, m, c);
        yNew = Y_RungeKutt(A, y, h, m, c);
        dyNew = yExperimental(:,i) - yNew;
        for iter = 1:M
            yGained(iter,i)=y(iter,1);
        end

        B = B + h*(U'*U + Unew'*Unew) / 2.0;
        b = b + h*(U'*dy + Unew'*dyNew) / 2.0;
        I = I + h*(dy'*dy + dyNew'*dyNew) / 2.0;
        U = Unew;
        y = yNew;
        dy = dyNew;
        i = i + 1;
    end
    delta = pinv(B)*b;
    c(4) = c(4) + delta(1);
    m(1) = m(1) + delta(2);
    m(3) = m(3) + delta(3);
    disp(I);
    figure
    plot(0:0.2:50, yExperimental(1,:),0:0.2:50, yGained(1,:));
end

disp(size(yGained));
figure
plot(0:0.2:50, yExperimental(1,:),0:0.2:50, yGained(1,:));
disp('c4: ');
disp(c(4));
disp('m1: ');
disp(m(1));
disp('m3: ');
disp(m(3));
disp('Iterations: ');
disp(numIter);
