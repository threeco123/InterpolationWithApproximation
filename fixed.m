function [] = fixed (x0,y0)
% Tolerance and format are set.
% Some preliminary variables are defined as well.
format long;
tol=10e-8;
LeNorme=1;
x0=2;
y0=0;
x_input=x0;
y_input=y0;
X0=[x_input;y_input];
k=0;
% Implemented the Fixed point algorithm discussed in class.
% To find the fixed point function, the equations were
% rearranged by hand to isolate x and y.
%
% 1st Iteration
% +sqrt and -sqrt
while (LeNorme > tol)
    k=k+1;
    X1=[-sqrt(4-y0*y0);1-sqrt(4-(4/9 *(x0-1)^2))];
    deltaX=X1-X0;
    LeNorme=norm(deltaX);
    X0=X1;
    x0=X1(1);
    y0=X1(2);
    if (k==1000)
        fprintf('Script unresponsive; too many iterations.');
        break;
    end 
end
S1=single(X1)
k1=k
%
% 2nd Iteration
% -sqrt and -sqrt
LeNorme=1;
X0=[x_input;y_input];
k=0;
%
while (LeNorme > tol)
    k=k+1;
    X1=[sqrt(4-y0*y0);1-sqrt(4-(4/9 *(x0-1)^2))];
    deltaX=X1-X0;
    LeNorme=norm(deltaX);
    X0=X1;
    x0=X1(1);
    y0=X1(2);
    if (k==1000)
        fprintf('Script unresponsive; too many iterations.');
        break;
    end 
end
S2=single(X1)
k2=k
%
return;
