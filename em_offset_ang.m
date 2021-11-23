function [e1off,e2off,e1fit,e2fit,anghxhy] = ...
    em_offset_ang(N,t,nstep,navg,E1,E2,Hx,Hy)
%EM_OFFSET_EST is used to do the least-squarcs fit within small windows,
%in order to estimate the time series of voltage offset
%The concept is the same as the default data processing method
%but focus on estimating the offset
%This function will generate the required input for the function
%em_vel_rotate.m

%e1off is the result of voltage offset on E1
%e2off is the result of voltage offset on E2
%e1fit is the fitted result of voltage on E1
%e2fit is the fitted result of voltage on E1
%anghxhy is the results of orientation of magnetometer


%N is the order of the polynomial fit, which should not be more than 4 in
%50-s windows

%t is time series of em measurements
%navg is the length of the individual window, = 50 s in the paper
%nstep is the step of moving windows, = 5 s for generating 10 realizations
%In other words, the number of overlapped measurements equals navg - nstep
%e1 is the raw data on E1
%e2 is the raw data on E2
%hx is the raw hx
%hy is the raw hy

%The current codes are written by Dr. J-Y Hsu in National Taiwan University
%on 11/23/2021
%%


mout = ceil((length(t) - navg) / nstep);
e1off   = NaN(length(t),mout);
e2off   = NaN(length(t),mout);
e1fit   = NaN(length(t),mout);
e2fit   = NaN(length(t),mout);

anghxhy   = NaN(length(t),mout);
%%
nout=1;
for i1 =1:nstep:length(t)-navg
    j = i1:i1+navg-1;
    tj = t(j);
    e1 = E1(j);
    e2 = E2(j);
    %%
    hx = Hx(j);
    hy = Hy(j);
    hx = hx - mean(hx); hx = hx / std(hx) / sqrt(2);
    hy = hy - mean(hy); hy = hy / std(hy) / sqrt(2);
    %%
    e1=e1(:);e2=e2(:);
    hx=hx(:);hy=hy(:);tj=tj(:);e1=e1(:);e2=e2(:);
    %%
    angs=atan2(hx,hy);
    anghxhy(j,nout)=angs;      
    %%
    con = ones(size(tj));
    trnd = tj - tj(1); % trend
    trnd = trnd / (max(trnd)-min(trnd)); 
    BASIS=cat(2,hx,hy,con);    
    for k=1:N
        BASIS=cat(2,BASIS,trnd.^k);        
    end
    COEF1 = BASIS \ e1; % least squares fit
    COEF2 = BASIS \ e2;
    e1b=con*COEF1(3); %Fitted result of background constant
    e2b=con*COEF2(3);    
    e1f=hx.*COEF1(1)+hy.*COEF1(2)+con*COEF1(3); %Fitted result of the curve
    e2f=hx.*COEF2(1)+hy.*COEF2(2)+con*COEF2(3);
    
    for k=1:N
        e1b=e1b+COEF1(3+k)*trnd.^k; %the revised results of voltage offset by adding trend
        e2b=e2b+COEF2(3+k)*trnd.^k;
        e1f=e1f+COEF1(3+k)*trnd.^k;
        e2f=e2f+COEF2(3+k)*trnd.^k;        
    end
    e1off(j,nout)=e1b;
    e2off(j,nout)=e2b;
    e1fit(j,nout)=e1f;
    e2fit(j,nout)=e2f;  

    %%
    nout=nout+1;
end

return

end

