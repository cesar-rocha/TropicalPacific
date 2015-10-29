% define desired region
region_name='TropicalPacific';
minlat=-25;
maxlat=-5;
minlon=200;
maxlon=280;

% extract indices for desired region
nx=4320;
prec='real*4';
gdir='/nobackupp8/dmenemen/llc/llc_4320/grid/';
fnam=[gdir 'Depth.data'];
[tmp fc ix jx] = ...
    quikread_llc(fnam,nx,1,prec,gdir,minlat,maxlat,minlon,maxlon);
m(1)=0;
for f=1:length(fc)
    m(f+1)=length(ix{fc(f)});
end
n=length(jx{fc(1)});
fld=zeros(sum(m),n);
for f=1:length(fc)
    fld((sum(m(1:f))+1):sum(m(1:(f+1))),:)=tmp{fc(f)};
end
quikpcolor(fld')

%% Get and save grid information
%pin='/nobackupp8/dmenemen/llc/llc_4320/grid/';
%%pout=['./grid/'];
%pout=['/nobackupp8/dmenemen/llc/llc_4320/regions/' region_name '/grid/'];
%
%eval(['mkdir ' pout])
%eval(['cd ' pout])
%for fnm={'AngleCS','AngleSN','DXC','DXG','DYC','DYG','Depth', ...
%         'RAC','RAS','RAW','RAZ','U2zonDir','V2zonDir', ...
%         'XC','XG','YC','YG','hFacC','hFacS','hFacW'}
%    fin=[pin fnm{1} '.data'];
%    fout=[fnm{1} '_' int2str(sum(m)) 'x' int2str(n)];
%    for f=1:length(fc)
%        fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
%            read_llc_fkij(fin,nx,fc(f),1,ix{fc(f)},jx{fc(f)});
%    end
%    writebin(fout,fld);
%end

% get and save regional fields
pin='/nobackupp8/dmenemen/llc/llc_4320/MITgcm/run/';
pout=['/nobackupp8/dmenemen/llc/llc_4320/regions/' region_name '/'];
%pout=['../'];

kx=1:1;
%for fnm={'Eta'}
for fnm={'Eta','U','V'}
    eval(['mkdir ' pout fnm{1}])
    eval(['cd ' pout fnm{1}])
    %for ts=10368:144:279360, disp(ts)
    for ts=279360:144:485568, disp(ts)

        fin=[pin fnm{1} '.' myint2str(ts,10) '.data'];
        dy=ts2dte(ts,25,2011,9,10,30);
        fout=[fnm{1} '_' int2str(sum(m)) 'x' int2str(n) 'x' int2str(length(kx)) '.' dy];
        for k=1:length(kx); mydisp(k)
            for f=1:length(fc)
                fld((sum(m(1:f))+1):sum(m(1:(f+1))),:) = ...
                    read_llc_fkij(fin,nx,fc(f),kx(k),ix{fc(f)},jx{fc(f)});
            end
            writebin(fout,fld,1,'real*4',k-1);
        end
    end
end
