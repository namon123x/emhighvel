function [u_r,v_r,verr1,verr2,E_S,E_E] = ...
    em_vel_spec_volt(mltref,fh,fz,e1,e2,wef,off_time,e1off,e2off,me1r,me2r,ang_time,anghxhy)
%em_vel_rotate is used to rotate the E1 and E2 to Cartesian coordinates to
%derive the true velocity field

% electrode separation
esep1 = (8+5/8)*0.0254; % m
c1 = 0.5;
c2 = -0.2;
% angles between compass and electrode axes
alpha2 = 1.95;
alpha1 = alpha2 - pi/2;
sfv1 = 1e3/(fz*esep1*(1.0+c1));
sfw = fh/fz*(1.0+c2)/(1.0+c1);
%%
%Change all variable into the same dimension
var_name=char('mltref','e1','e2','wef','off_time','e1off','e2off','ang_time','anghxhy');
for N=1:length(var_name(:,1))
    vars=strtrim(var_name(N,:));
    eval([vars '=' vars '(:);'])
end
%%
ks=~isnan(off_time)&~isnan(e1off)&~isnan(e2off);
% off_time=off_time(ks);e1off=e1off(ks);e2off=e2off(ks);

% e1off(ks)
% 
% stop
e1_rv=interp1(off_time(ks),e1off(ks),mltref);
e2_rv=interp1(off_time(ks),e2off(ks),mltref);
% stop
e1f=e1-e1_rv;

e2f=e2-e2_rv;
ks=~isnan(ang_time)&~isnan(anghxhy);
angi=unwrap(anghxhy);
ang_rs=interp1(ang_time(ks),angi(ks),mltref);
ang_r=atan2(sin(ang_rs),cos(ang_rs));


%%
%Rotate electrode coordinates to the Cartesian coordinates
Ex_r=e1f.*cos(ang_r)-e2f.*sin(ang_r);
Ey_r=e1f.*sin(ang_r)+e2f.*cos(ang_r);


%%
%Adjust the angle between Hx and E1 in order put it back to
%Cartesian coordinates
%H_S means this Jx is toward south
s = sin(alpha1);
c = cos(alpha1);
E_S=Ex_r*c-Ey_r*s;
E_E=Ex_r*s+Ey_r*c;


%%

u_r=-E_E*sfv1;
v_r=-E_S*sfv1 + wef * sfw;
if isempty(find(isnan(me1r), 1)) || isempty(find(isnan(me2r), 1)) 
    verr1=NaN;verr2=NaN;
    return

else
verr1=me1r*(sfv1.^2);
verr2=me2r*(sfv1.^2);
end

return
end

