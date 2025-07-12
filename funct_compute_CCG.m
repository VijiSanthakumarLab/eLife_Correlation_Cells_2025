function [CCG,bins] = compute_CCG(event_times_A,event_times_B)

% Computes the cross-correlogram between two event time vectors and plots
% the +/- 100 ms normalized cross-correlation. Discrete events are
% convolved with a gaussian of width 20 ms (see Parameters).
% Also computes the shuffled control, but does not plot or output (can be
% modified as needed)
%
% Inputs:
% event_times_A - vector of timestamps for event A (in seconds)
% event_times_B - vector of timestamps for event B (in seconds)
%
% Ouputs:
% CCG - output of the MATLAB 'xcorr' function, normalized cross-correlation
% bins - bin lags in milliseconds


% Parameters
maxlag=1000;
gauss_width=20;
norm=1;
iter=100;

x=event_times_A;
y=event_times_B;

N_x=length(x);
N_y=length(y);

time=1000*(round(max([x;y]))+1); %in milliseconds

x_dig=round(1000*x);
y_dig=round(1000*y);

x_samples=zeros(1,time);
x_samples(x_dig)=1;
y_samples=zeros(1,time);
y_samples(y_dig)=1;

gausskernel=gausswin(gauss_width);

x_sm=conv(x_samples,gausskernel);
y_sm=conv(y_samples,gausskernel);

x_sm=x_sm(1:time);
y_sm=y_sm(1:time);

if norm==1
    C = xcorr(y_sm,x_sm,maxlag,'normalized');
else
    C = xcorr(y_sm,x_sm,maxlag);
end

C_all=zeros(iter,maxlag*2+1);

for (i=1:iter)

    x_unique=randperm(time);
    y_unique=randperm(time);

    x=zeros(1,time);
    x(x_unique(1:N_x))=1;

    y=zeros(1,time);
    y(y_unique(1:N_y))=1;

    gausskernel=gausswin(gauss_width);

    x_sm=conv(x,gausskernel);
    y_sm=conv(y,gausskernel);

    x_sm=x_sm(1:time);
    x_sm=x_sm-mean(x_sm);
    y_sm=y_sm(1:time);
    y_sm=y_sm-mean(y_sm);

    if norm==1
        C_err = xcorr(x_sm,y_sm,maxlag,'normalized');
    else
        C_err = xcorr(x_sm,y_sm,maxlag);
    end

    C_all(i,:)=C_err;

end

figure
plot(-1000:1:1000,C)
xlim([-100 100])

CCG=C;
bins=-1000:1:1000;

end