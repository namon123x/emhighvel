function [tverr,pverr,verr1,verr2] = verr_construct(fz,e1,e2,e1fit,e2fit,t,Pef,navg,nstep)
%VERR_CONSTRUCT is used for constructing the verr from the rotation axes
%method
%verr is the variance of the residuals in the least-squares fit
%See functions em_vel_rotate.m and em_offset_ang.m for the descriptions of
%the variable
%

%The current codes are written by Dr. J-Y Hsu in National Taiwan University
%on 11/23/2021



c1=0.5;
esep1 = (8+5/8)*0.0254; % m, the length of electrode pairs
esep2 = (8+5/8)*0.0254; % m
sfv1 = 1e3/(fz*esep1*(1.0+c1));
sfv2 = 1e3/(fz*esep2*(1.0+c1));


mout = ceil((length(t) - navg) / nstep);
e1r   = NaN(1,mout);
e2r   = NaN(1,mout);
pverr= NaN(1,mout);
tverr= NaN(1,mout);
%%
nout=1;
for i1 =1:nstep:length(t)-navg
    j = i1:i1+navg-1;
    e1 = e1(j);
    e2 = e2(j);
    
    e1f=e1fit(j);
    e2f=e2fit(j);
    e1r(nout)=var(e1(:)-e1f(:));
    e2r(nout)=var(e2(:)-e2f(:));    
    pverr(nout)=mymean(Pef(j));
    tverr(nout)=mymean(t(j));
    nout=nout+1;
end

verr1=e1r*(sfv1.^2);
verr2=e2r*(sfv2.^2);


return

end

