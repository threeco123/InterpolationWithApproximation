function [] = broyden2()    
    %
    % Define a function that carries out
    % Gaussian Elimination with Partial Pivoting
    % and Back Substitution
    %
    function [F] = LeFonction(u0,v0)
    fu0=u0^2+v0^2-4;
    fv0=((u0-1)^2)/9 + ((v0-1)^2)/4 -1;
    F=[fu0;fv0];
    end
    %
    % Define a function that carries out
    % Gaussian Elimination with Partial Pivoting
    % and Back Substitution
    %
    function [X] = Gauss (A,b)
    %
    % Partial Pivot: If the second row is bigger than 
    % the first row, swap them. The code is not scalable to
    % dimensions greater than 2.
    %
    if (A(2,1) > A(1,1))
        A([1,2],:) = A([2,1],:);
        b([1,2],:) = b([2,1],:);
    end
    %
    % Gaussian Elimination using the algorithm found
    % in the class notes
    %
    mult=A(2,1)/A(1,1);
    A(2,1)=A(2,1)-mult*A(1,1);
    A(2,2)=A(2,2)-mult*A(1,2);
    b(2)=b(2)-mult*b(1);
    %
    % Back Substitution. Again, the code is not very scalable.
    %
    x2=b(2)/A(2,2);
    x1=(b(1)-A(1,2)*x2)/A(1,1);
    X=[x1 x2];
    end
%
% Solutions: These were found using Wolfram Alpha.
% It is important to know the solutions in order to 
% determine if we are getting close to them.
%
x_ans1=-1.9229733785204901811;
y_ans1=0.54970299753729885611;
x_ans2=1.7689438563170703178;
y_ans2=-0.93318681580811678487;
%
% The format is set to long to improve accuracy. The norm
% and tolerance are initialized. 
%
format long
LeNorme=1;
tol=10^(-8);
%
% 'm' is the number of steps, while b is the length
% of the square we will consider. Thus the step size will
% be 2b/m. L will contain the color information for each
% point in the square. x and y are all the points that
% will be iterated on.
%
m=256;
b=3;
L=zeros(m,m,3);
x=linspace(-b,b,m);
y=linspace(-b,b,m);
% 
% The square region from -3 to 3 is scanned with a step size
% of 2b/m.
% 
for p=1:m
    for q=1:m
        l=0;
        x0=x(p);
        y0=y(q);
        LeNorme=1;
        % Define the Jacobian and the iteration matrix A0
        A0=[1 0;0 1];
        %J0=[2*x0, 2*y0; 2*(x0-1)/9, 2*(y0-1)/4];
        %A0=J0;
        % The algorithm for Broyden's method found in the class notes.
        while (1)
            l=l+1;
            F0=LeFonction(x0,y0);
            DeltaX=Gauss(A0,-F0);
            x0=x0+DeltaX(1);
            y0=y0+DeltaX(2);
            F1=LeFonction(x0,y0);
            DeltaF=F1-F0;
            A0=A0+(((DeltaF-A0*DeltaX')*(DeltaX))/(norm(DeltaX)));
            LeNorme=norm(DeltaX);
            % Max. numebr of iterations is 200.
            if (l==200)
                L(q,p,3)=1;
                break
            end
            % If the tolerance is reached...
            if (LeNorme< tol)  
                % ...and if we have converged to the first solution...
                if (abs(x0-x_ans1)<tol && abs(y0-y_ans1)<tol)
                % ...set the color to red, depending on the speed
                % of convergence. Notice that 'l' is the number of
                % iterations.
                    L(q,p,1)=(l/200)^(1/3);
                elseif (abs(x0-x_ans2)<tol && abs(y0-y_ans2)<tol)
                    L(q,p,2)=(l/200)^(1/3);
                else
                % Otherwise, mark it as blue.
                L(q,p,3)=1;
                end
                break
            end
        end
    end
end
% 
% Initialize colormap matrix. Define r, which is the 
% darkest color in the colormap. Set k to m/2 for convenience.
% 
colors=zeros(m-1,3);
r=(1/200)^(1/3);
k=m/2;
% 
% Between 1 and m/2, set the colors to red. The formula used
% below is found by knowing that at 1, color must be 1, (brightest),
% while at g=m/2 we must have the darkest shade, which has a value of r
% 
for g=1:k;
    colors(g,1)=(g*(r-1))/(k-1)+1-((r-1)/(k-1));
end
% 
% Similar to above, only now we need dark green at the beginning,
% and bright green at the end.
% 
for h=k+1:m;
    colors(h,2)=(h*(r-1))/(k+1-m)+r-(k+1)*(r-1)/(k+1-m);
end
% 
% Display the computed answers using imagesc
% 
imagesc(x,y,L)
axis xy
colormap(colors)
colorbar
caxis([-200, 200])
end