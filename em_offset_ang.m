function [e1off,e2off,e1fits,e2fits,me1r,me2r,anghxhy] = ...
    em_offset_ang(N,t,nstep,navg,E1,E2,Hx,Hy)
%EM_OFFSET_EST is used to do the harmonic within small windows,
%in order to estimate the time series of offset and instrumental noise
%The concept is the same as the default data processing method
%but focus on estimating the offset

%N is the order of the polynomial fit

%t is time series of em measurements
%navg is the length of the individual window
%nstep is the step of moving windows
%e1 is the raw data on E1
%e2 is the raw data on E2
%hx is the raw hx
%hy is the raw hy
%sm_scale is the scale for the Gaussian smoothing (sec)
%drange set up the depth criteria for estimating noise level

no_a1a2=0; %whether to fit the sinusoidal function of constant current



%%
mout = ceil((length(t) - navg) / nstep);
e1off   = NaN(length(t),mout);
e2off   = NaN(length(t),mout);
e1fits   = NaN(length(t),mout);
e2fits  = NaN(length(t),mout);
timeref  = NaN(1,mout);
me1ref  = NaN(1,mout);
me2ref   = NaN(1,mout);
anghxhy   = NaN(length(t),mout);
%%
nout=1;
for i1 =1:nstep:length(t)-navg
    j = i1:i1+navg-1;
    tj = t(j);
    timeref(nout)=mean(tj);
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
    if no_a1a2
        BASIS=con;
    else
        BASIS=cat(2,hx,hy,con);
    end
    for k=1:N
        BASIS=cat(2,BASIS,trnd.^k);
    end
    COEF1 = BASIS \ e1; % least squares fit
    COEF2 = BASIS \ e2;
    if no_a1a2 == 1
        e1b=con*COEF1(1); %Fitted result of background
        e2b=con*COEF2(1);
        e1f=con*COEF1(1); %Fitted result of the curve
        e2f=con*COEF2(1);
        for k=1:N
            e1b=e1b+COEF1(1+k)*trnd.^k;
            e2b=e2b+COEF2(1+k)*trnd.^k;
            e1f=e1f+COEF1(1+k)*trnd.^k;
            e2f=e2f+COEF2(1+k)*trnd.^k;
        end
    elseif no_a1a2==0
        e1b=con*COEF1(3); %Fitted result of background
        e2b=con*COEF2(3);
        e1f=hx.*COEF1(1)+hy.*COEF1(2)+con*COEF1(3); %Fitted result of the curve
        e2f=hx.*COEF2(1)+hy.*COEF2(2)+con*COEF2(3);
        for k=1:N
            e1b=e1b+COEF1(3+k)*trnd.^k;
            e2b=e2b+COEF2(3+k)*trnd.^k;
            e1f=e1f+COEF1(3+k)*trnd.^k;
            e2f=e2f+COEF2(3+k)*trnd.^k;
        end
    end



    e1off(j,nout)=e1b;
    e2off(j,nout)=e2b;
    e1fits(j,nout)=e1f;
    e2fits(j,nout)=e2f;    
    me1ref(nout)=var(e1-e1f);
    me2ref(nout)=var(e2-e2f);
    %%
    nout=nout+1;
end
if isempty(timeref) ==0
    me1r=interp1(timeref,me1ref,t);
    me2r=interp1(timeref,me2ref,t);
else
    me1r=NaN;
    me2r=NaN;
end

return

end

