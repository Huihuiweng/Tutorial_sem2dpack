% SEM2D_READ_FAULT reads fault outputs from SEM2DPACK
%
% SYNTAX	data = sem2d_read_fault(name)
%
% INPUT		name	[Flt*] 	prefix of header and data files (name_sem2d.*) 
%				The default is the first FltXX_sem2d.* found 
%				in the current directory.
%
% OUTPUT	data.nx		number of fault nodes
%		data.nt		number of time samples
%		data.dt		time step
%		data.x,data.z	coordinates of fault nodes
%		data.d		slip [nx,nt]
%		data.v		slip rate
%		data.st		shear stress
%		data.sn		normal stress
%		data.mu		friction coefficient
%		data.st0	initial value of shear stress
%		data.sn0	initial value of normal stress
%		data.mu0	initial value of friction coefficient
%		data.P   	pore pressure
%		data.T   	temperature
% 
%		If output on each side of the fault (osides=T):
%  		data.d1t	displacement on side 1, fault parallel component
%  		data.d2t	displacement on side 2, fault parallel component
%  		data.d1n	displacement on side 1, fault normal component
%  		data.d2n	displacement on side 2, fault normal component
%  		data.v1t	velocity on side 1, fault parallel component
%  		data.v2t	velocity on side 2, fault parallel component
%  		data.v1n	velocity on side 1, fault normal component
%  		data.v2n	velocity on side 2, fault normal component
%
% NOTE		Fault normal components on each side of the fault are exported only in P-SV
%
function data = sem2d_read_fault(name)

% length of the tag at the begining and end of a binary record
% in number of single precision words (4*bytes)
LENTAG = 2; % gfortran older versions
LENTAG = 1;

% assumes header file name is FltXX_sem2d.hdr
if ~exist('name','var')
  list = dir('Flt*.hdr');
  list = {list.name};
  if isempty(list)
    name = '';
  else
    name=list{1}(1:5);
  end
end

% Read parameters from header file
hdr = strcat(name,'_sem2d.hdr');
if ~exist(hdr,'file')
  data=[]; 
  warning(['File ' hdr ' not found'])
  return
end
[data.nx,ndat,data.nt,data.dt] = textread(hdr,'%n%n%n%n',1,'headerlines',1);
variables = textread(hdr,'%s',1,'headerlines',2,'delimiter','\n');
variables = split(variables,':');
[data.x,data.z] = textread(hdr,'%f%f','headerlines',4);

% Read initial fault data
if exist([name '_init_sem2d.tab'],'file')
  raw = load([name '_init_sem2d.tab']);
  data.st0 = raw(:,1);
  data.sn0 = raw(:,2);
  data.mu0 = raw(:,3);
end

% Read fault data in a big matrix
dat  = strcat(name,'_sem2d.dat');
fid=fopen(dat); 
raw = fread(fid,[data.nx+2*LENTAG,inf],'single') ; 
fclose(fid);
%raw = reshape(raw(2:data.nx+1,:),[data.nx ndat data.nt]);
raw = reshape(raw(LENTAG+1:LENTAG+data.nx,:),data.nx,ndat,[]);


% Reformat each field [nx,nt]
for i=1:numel(variables)
    variable = char(variables(i));
    switch variable
        case 'Slip'
            data.d  = squeeze(raw(:,i,:)); 
        case 'Slip_Rate'
            data.v  = squeeze(raw(:,i,:)); 
        case 'Shear_Stress'
            data.st = squeeze(raw(:,i,:)); 
        case 'Normal_Stress'
            data.sn = squeeze(raw(:,i,:)); 
        case 'Friction'
            data.mu  = squeeze(raw(:,i,:)); 
        case 'D1t'
            data.d1t  = squeeze(raw(:,i,:)); 
        case 'D2t'
            data.d2t  = squeeze(raw(:,i,:)); 
        case 'V1t'
            data.v1t  = squeeze(raw(:,i,:)); 
        case 'V2t'
            data.v2t  = squeeze(raw(:,i,:));
        case 'D1n'
            data.d1n  = squeeze(raw(:,i,:)); 
        case 'D2n'
            data.d2n  = squeeze(raw(:,i,:)); 
        case 'V1n'
            data.v1n  = squeeze(raw(:,i,:)); 
        case 'V2n'
           data.v2n  = squeeze(raw(:,i,:)); 
        case 'Pore_Pressure'
           data.P = squeeze(raw(:,i,:)); 
        case 'Temperature'
           data.T = squeeze(raw(:,i,:));
    end
end