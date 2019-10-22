clear
clc

eps = 1e-12;

X = imread("x3.bmp");
figure(1,"position",get(0,"screensize"))
imshow(X);
title("X input");
X = double(X);
X=[X;ones(1,columns(X))];
 
Y = imread("y8.bmp");
figure(2,"position",get(0,"screensize"))
imshow(Y);
title("Y input");
Y = double(Y);
 
function invertedA = pseudoInverse(A, delta)
  columns = columns(A);
  rows = rows(A);
  if(rows > columns)
    invertedA = inv(A'*A + delta^2*eye(columns))*A';
  else
    invertedA = A'*inv(A*A' + delta^2*eye(rows));
  endif
endfunction

function retval = Greville(A, eps)
  vector = A(1,:)';
  rows = rows(A);
  columns = columns(A);
 
  if(vector.*vector' < eps)
    invertedA = vector;
  else
    invertedA = vector/(vector'*vector);
  endif

  for i = 2:rows
    vector = A(i,:)';
    Z = eye(columns) - invertedA * A(1:i-1,:);
    norm = vector'*Z*vector;
    if(norm < eps)
      Z = invertedA*invertedA';
      norm = 1+vector'*Z*vector;
    endif
    invertedA=invertedA-Z*vector*vector'*invertedA/norm;
    invertedA=[invertedA,Z*vector/norm];
  endfor
  retval = invertedA;
endfunction

V = rand(rows(Y), rows(X));

Greville_X = Greville(X, eps);
Greville_Z = eye(rows(X)) - X*Greville_X;
Greville_A = Y*Greville_X + V*Greville_Z';
Greville_Y = Greville_A*X;
  
figure(4,"position",get(0,"screensize"))
imshow(uint8(Greville_Y));
title("Greville");
 
function retval = MoorePenrose(A, eps)
  delta = 100;
  eps = 1e-12;
  diff = 1;

  invertedA1 = pseudoInverse(A,delta);
 
  while(diff > eps)
    delta = delta / 2;
    invertedA2 = pseudoInverse(A,delta);
    diff = norm(invertedA1-invertedA2);
    invertedA1 = invertedA2;
  endwhile

  retval = invertedA2;
endfunction
  
MoorePenrose_X = MoorePenrose(X,eps);
MoorePenrose_Z = eye(rows(X)) - X*MoorePenrose_X;
MoorePenrose_A = Y*MoorePenrose_X + V*MoorePenrose_Z';
MoorePenrose_Y = MoorePenrose_A*X;
  
figure(3,"position",get(0,"screensize"))
imshow(uint8(MoorePenrose_Y));
title("Moore-Penrose");