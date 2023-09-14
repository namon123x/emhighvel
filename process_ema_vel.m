%This script is an example for processing the EM voltage measurements via
%rotating axes method (RAM)
%If you need codes for decoding the .csv files of APF-11 floats at matlab files, 
%please contact Dr. Hsu at jyahsu@ntu.edu.tw
%
%In the output file, 
%prof includes the horizontal current velocity results after bin average,
%i.e., the official u and v from RAM
%
%rotpr contains information before bin average
%
%This version is written by Dr. J.-Y. Hsu at 09.13.2023


clear

%WMM coefficient can be found at https://www.ncei.noaa.gov/products/world-magnetic-model
wmmfile='E:\Main_Work\Data_Reposit\WMM\WMM2020\WMM.COF';


%%
navg=50; %50-s windows
nstep=5; %5-s increment
tstep=12; %temporal-averaging interval



tsbin=0; %if choose not using the bin result of science log for processing
%%



%%
b=1;

emaname='f9207_009_ema_log.mat'; % ema log file as example
sciname='f9207_009_science_log.mat'; %science log file as example
E1_as=NaN;
load(emaname)

time_ema=mtime;
id1=strfind(emaname,'_');
cycles=emaname(id1+1:id1+3);
cycles
load(sciname)
%%

if tsbin %use bin profile of T and S
    timeds=mtime_bin_ds;
    Pds=P_bin_ds;
    Tds=T_bin_ds;
    Sds=S_bin_ds;
    timeas=mtime_bin_as;
    Pas=P_bin_as;
    Tas=T_bin_as;
    Sas=S_bin_as;
else
    timeds=mtime_ts_ds;
    Pds=P_ts_ds;
    Tds=T_ds;
    Sds=S_ds;
    timeas=mtime_ts_as;
    Pas=P_ts_as;
    Tas=T_as;
    Sas=S_as;

end
times=mean(time_gps,'omitnan');
B= wrldmagm(0, mean(lat,'omitnan'),mean(lon,'omitnan'),...
    decyear(datetime(datestr(times))),'CUSTOM',wmmfile);
Bx=B(1);By=B(2);Bz=B(3);
fz=-mean(Bz,'omitnan');
fh=mean(sqrt(Bx.^2+By.^2),'omitnan');

Wctd=gradient(-Praw,mtime_Praw)/86400;
%%
for p=1:2 %1 for descending and 2 for ascending
    if p==1
        Pef=interp1(mtime_Praw_ds,Praw_ds,mtime); %Get ema's depth based on pressure-time measurements
        km=~isnan(Pef);
        mlt_efr=mtime(km);
        e1=E1_mean(km);e2=E2_mean(km);
        hx=Hx_mean(km);hy=Hy_mean(km);
        Pef=Pef(km);
        timectd=timeds;
        Pctd=Pds; %the difference to Praw_ds is that Pctd is the depth with temperature and salinity
        T=Tds;
        S=Sds;
    elseif p==2
        Pef=interp1(mtime_Praw_as,Praw_as,mtime);
        km=~isnan(Pef);
        mlt_efr=mtime(km);
        e1=E1_mean(km);e2=E2_mean(km);
        hx=Hx_mean(km);hy=Hy_mean(km);
        Pef=Pef(km);
        timectd=timeas;
        Pctd=Pas;
        T=Tas;
        S=Sas;
    end


    if isempty(Pef)
        continue
    end



    %%
    W=gradient(-Pef,mlt_efr)/86400; %compute vertical speed of floats


    [sv_e1fp,sv_e2fp,~,~,me1r,me2r,angr] = ...
        em_offset_ang(4,mlt_efr,nstep,navg,e1,e2,hx,hy); %estimate voltage offset and orientation

    mang=nan(length(mlt_efr),1);
    for i=1:length(mlt_efr)
        angi=mean(unwrap(angr(i,:)),'omitnan');
        mang(i)=atan2(sin(angi),cos(angi));
    end
    emnum=nan(length(mlt_efr),1);
    for i=1:length(mlt_efr)
        e_ens=sv_e1fp(i,:);
        emnum(i)=length(e_ens(~isnan(e_ens)));
    end
    me1p=meanNaN(sv_e1fp,2);
    se1p=stdNaN(sv_e1fp,2);
    me2p=meanNaN(sv_e2fp,2);
    se2p=stdNaN(sv_e2fp,2);

    %                         stop
    %%
    %Aritrarily remove results in the top and bottom due to
    %insufficient realizations for estimated offset
    koff=emnum>=5;
    paraname=char('mlt_efr','Pef','e1','e2','W',...
        'me1p','me2p','se1p','se2p','me1r','me2r','mang');
    for i=1:length(paraname(:,1))
        para=strtrim(paraname(i,:));
        eval([para '=' para '(koff);'])
    end


    if isempty(Pef(Pef<100))
        continue
    end
    if isempty(e1)
        continue
    end
    %%
    [u_r,v_r,verr1i,verr2i,e_S,e_E] = ...
        em_vel_spec_volt(mlt_efr,fh,fz,e1,e2,W,mlt_efr,...
        me1p,me2p,me1r,me2r,mlt_efr,mang); %RAM for deriving 1-s seawater velocity



    %%
    %Bin-averaging results in independent bins of 12 s
    %variable with sm_xxx means the temporal-averaging results
    EN=-e_S;
    EE=e_E;
    time_sm=mlt_efr(1):tstep/86400:mlt_efr(end);
    [sm_me1p,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,me1p,0,'std');
    [sm_me2p,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,me2p,0,'std');
    [sm_se1p,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,se1p,0,'std');
    [sm_se2p,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,se2p,0,'std');
    [sm_ur,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,u_r,0,'std');
    [sm_Pr,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,Pef,0,'std');
    [sm_vr,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,v_r,0,'std');
    [meant,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,mlt_efr,0,'std');



    Verr1=interp1(mlt_efr,verr1i,meant);
    Verr2=interp1(mlt_efr,verr2i,meant);

    ge1p=(gradient(me1p,mlt_efr))/86400;
    ge2p=(gradient(me2p,mlt_efr))/86400;
    [sm_ge1p,sm_sge1p,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,ge1p,0,'std');
    [sm_ge2p,sm_sge2p,~]=bin_aver(time_sm,0.5*tstep/86400,mlt_efr,ge2p,0,'std');


    div_ids(b)=str2double(cycles);
    nids(b)=b; %define as the ID of profiles
    prof(b).fh=fh; %horizontal magnetic field
    prof(b).fz=fz; %vertical magnetic field

    prof(b).time=time_sm;
    prof(b).u=sm_ur;
    prof(b).v=sm_vr;
    prof(b).P=sm_Pr;



    %The four following variables are used for quality control
    prof(b).sm_sge1p=sm_sge1p;% std of temporal change of E1 offset
    prof(b).sm_sge2p=sm_sge2p;% std of temporal change of E2 offset
    prof(b).sm_se1p=sm_se1p;% scattering of E1 offset
    prof(b).sm_se2p=sm_se2p;% scattering of E2 offset


    %CTD measurements
    prof(b).timectd=timectd;
    prof(b).Pctd=Pctd;
    prof(b).T=T;
    prof(b).S=S;




    rotpr(b).EE=EE; %motion-induced voltage in zonal
    rotpr(b).EN=EN;%motion-induced voltage in meridional
    rotpr(b).mlt=mlt_efr; %time
    rotpr(b).ur=u_r; %u
    rotpr(b).vr=v_r; %v
    rotpr(b).pr=Pef; %pressure
    rotpr(b).me1offset=me1p; %mean offset on e1
    rotpr(b).me2offset=me2p; %mean offset on e2
    rotpr(b).se1offset=se1p; %std offset on e1
    rotpr(b).se2offset=se2p; %std offset on e2
    rotpr(b).verr1=Verr1; %residuals of the fit (not used)
    rotpr(b).verr2=Verr2; %residuals of the fit (not used)
    rotpr(b).mang=mang; %mean of orientation




    b=b+1;

end




