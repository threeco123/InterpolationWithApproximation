function [] = granular (filename,N2,N3,M)
% Syntax: granular('testing.wav',10000,20000,15000)
%
%
% This is the window function. It takes an empty array as
% input and fills it with elements from the original sample.
% The information for this window function is found in the
% project description.
% 
function [C] = window(C)
% The first interval
lower=N1;
higher=N1+N2;
for i=lower:higher;
    % This 'if' statement is needed in order to avoid
    % going beyond the limits of the array 'Y'
    if i>=numel(Y)
        break
    end
    C(i-M*(k-1)) = Y(i)*((i-N1)/N2);
end
% The second interval
lower=N1+N2;
higher=N1+N2+N3;
for i=lower:higher;
    if i>=numel(Y)
        break
    end
    C(i-M*(k-1)) = Y(i);
end
% The third interval
lower=N1+N2+N3;
higher=N1+2*N2+N3;
for i=lower:higher;
    if i>=numel(Y)
        break
    end
    C(i-M*(k-1)) = Y(i)*((N1+2*N2+N3-i)/N2);
end  
end
%
% Reads the file and makes sure it is monophonic
% by only retaining the first column of samples
%
Y=wavread(filename);
Y=Y(:,1);
% The windows have the same length
wlength=2*N2+N3;
% 
wavsize=size(Y);
wavsize=wavsize(1,1);
% The number of windows is given by the following formula.
% This formula was determined by inspection.
N_windows=(wavsize-(wlength))/M +1;
N_windows=N_windows(:,1);
N_windows=round(N_windows);
% 
% Initializes the matrix that contains all the windows. 
% Each column represents the data points of 1 window.
%
W=zeros(wlength,N_windows);
% 
% In order to start reading windows from the beginning of the audio signal,
% we need to set N1=1. This value will then be updated by M
% during the following loop in order to get subsequent windows.
%
N1=1;
%
for k=1:N_windows;
    % Initializes a temporary array that contains the data points
    % for the current window.
    A=zeros(wlength,1);
    % This empty array is filled up by calling on the window function.
    A=window(A);
    % The elements contained in A are now copied into the appropriate 
    % column of W.
    for l=1:wlength;
        W(l,k)=A(l);
    end
    % The next window is found by shifting the index of N1 by M
    N1=N1+M;
end
% 
% The number 'z' is equal to the total length of the stretched audio
% sample. This formula was found by inspection.
% 
z=N2+(N2+N3)*N_windows;
% 
% Initializes the stretched vector. It is 1-dimensional just like
% the original sample 'Y'.
% 
B=zeros(z,1);
% These values represent the starting and ending points in 'B'.
% Initially, we will put the first window at the beginning of 'B'
smallest=1;
biggest=wlength;
% The following loops will recombine all the windows into 
% the stretched signal.
for m=1:N_windows;
    % The following formula was again foudn by inspection.
    for j=smallest:biggest;
        B(j)=B(j)+W(j-(m-1)*(N2+N3),m);
    end
    % The next window is not put at the beginning of 'B', but 
    % rather where the last window left off.
    smallest=smallest+N2+N3;
    biggest=biggest+N2+N3;
end
%
% Define time domains
time1=(1/44100)*(wavsize);
time2=(1/44100)*z;
t1=linspace(0,time1,wavsize);
t2=linspace(0,time2,z);	
% Plots original sample
subplot(2,1,1)
plot(t1,Y)
xlabel('Time(s)');
ylabel('Amplitude');
title('Original Signal before Time-stretching');
axis([0 time1 -1 1])
% Plots stretched sample
subplot(2,1,2)
plot(t2,B)
xlabel('Time(s)');
ylabel('Amplitude');
title('Signal after Granular Synthesis');
axis([0 time1 -1 1])
% Outputs the stretched audio file
wavwrite(B,44100,'granular.wav')
% Outputs the stretching factor alpha
alpha=size(B)/wavsize;
% Outputs the time-stretching factor
alpha=alpha(1,1)
end
