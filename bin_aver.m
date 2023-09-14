%type='se' Standard error
%type='std' Standard deviation
%y is the data will be bin averaged
%x is the domain falls on the bin
%Nz is the number in each bin

function [mey,erry,Nz] = bin_aver( bin,db,x,y,lcutx,type)

[xsor,ind]=sort(x);

ysor=y(ind);
nb=length(bin);

km=xsor<lcutx;
xsor(km)=NaN;
ysor(km)=NaN;

mey=nan(nb-1,1);
erry=nan(nb-1,1);
Nz=nan(nb-1,1);
for i=1:nb

    kr=xsor>=bin(i)-db&xsor<bin(i)+db;
    miux=xsor(kr);
    miy=ysor(kr);
    if isempty(miux)==1
        continue
    end

    mey(i)=mymean(miy);

    nz=length(miux);
    Nz(i)=nz;
    if strcmp(type,'se')
        erry(i)=stdNaN(miy(:),1)/(sqrt(nz));
        %             erry(i)=stdNaN(col(inty),1)/(sqrt(length(inty)));
    elseif strcmp(type,'std')
        erry(i)=stdNaN(miy(:),1);
    end
end

% mey=interp1(bin(~isnan(mey)),mey(~isnan(mey)),bin);
% erry=interp1(bin(~isnan(erry)),erry(~isnan(erry)),bin);
end

