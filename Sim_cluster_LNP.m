
clear all
%clear mex
mex LNP_ClusterNet3.c
fp = mfilename('fullpath');

if ~isempty(fp)
    % If this file is located in the repo root:
    repo_root = fileparts(fp);

    % If this file is inside a subfolder (e.g., repo_root/matlab_code/this_file.m),
    % then use one more fileparts:
    % repo_root = fileparts(fileparts(fp));

else
    % Run Selection / interactive execution case: mfilename is empty.
    % Use the call stack to locate the file being executed.
    st = dbstack('-completenames');

    if ~isempty(st) && isfield(st(1), 'file') && ~isempty(st(1).file)
        % st(1).file points to the current file on disk
        repo_root = fileparts(st(1).file);

        % If this file is inside a subfolder (e.g., repo_root/matlab_code/),
        % then use one more fileparts:
        % repo_root = fileparts(fileparts(st(1).file));
    else
        % Last-resort fallback: current working directory
        repo_root = pwd;
    end
end

% Optional one-time debug (uncomment to verify, then re-comment)
% disp("repo_root = " + string(repo_root))

data_folder   = fullfile(repo_root, 'unstruc_no_overlap/');
stimID     = fullfile(data_folder, 'stimID.txt');
contextID  = fullfile(data_folder, 'contextID.txt');
stimstr = fullfile(data_folder, 'stim_str.txt');  % random photostim strength
I_cont = fullfile(data_folder, 'I_cont.txt');  % context current

sound_str = 5;

% I_cont = gamrnd(5,1,2000,1); 
% stimstr = gamrnd(7,1,1600,1); 
% save("I_cont.txt","I_cont","-ascii","-double")
% save("stim_str.txt","stimstr","-ascii","-double")

for condition =1:4 

switch condition
    case 1
        task='LeftSoundI'; modulation=''; % sound only, passive 
    case 2
        task='LeftSSI'; modulation=''; % sound + stim, passive
    case 3
        task='LeftSoundI'; modulation='Ibias'; % sound only, active 
    case 4
        task='LeftSSI'; modulation='Ibias'; % sound + stim, active
end


if strcmp(modulation,'')
filename=@(trial) strrep(sprintf('%s2clusterLNP_%s_ID%d',data_folder,task,trial),'.','d'),
else
filename=@(trial) strrep(sprintf('%s2clusterLNP_%s_%s_ID%d',data_folder,modulation,task,trial),'.','d'),
end

% W_fname=sprintf('%sweight_2cluster_orig',data_folder), 
%W_fname=sprintf('%sweight_2clusterLNP_2',data_folder),
W_fname='weight_2clusterLNP_unstruc';
%W_fname = sprintf('%sweight_2clusterLNP_uns',data_folder);
load(W_fname,'PostID','out_degree','Wrr_c');

Ne=6800;   %6800;   %4000; 
Ni=1200;   %1200;   %1000
N=Ne+Ni; 

switch modulation
    case ''
        param.bias= [5+0*randn(Ne,1); 5+0*randn(Ni,1)];  %passive, originally 30,30
    case 'Ibias'
        param.bias= [5+0*randn(Ne,1); 5+0*randn(Ni,1)] + 1*randinput(contextID, N, I_cont);  %active, originally 20, 10 
end

for trial = 1:1000 %1:100 % (1:10) + 10*(job_dex-1)
fprintf('condition %d, trial %d\n',condition,trial)

taue=20; %20; % membrane time constant
taui=10; %10; 

param.Ne=Ne; 
param.Ni=Ni; 
param.N=N; 
param.tau_m= [taue*ones(Ne,1); taui*ones(Ni,1)];  % [E; I] (ms)

param.tau_syn= 2;  % (ms)  
% param.phi=@(v) 0.3*v.^2.*(v>Vth); % Hz  (need to change in C code) 
% param.phi=@(v) 30./(1+exp(-(v-5))); 
param.phi=@(v) 30./(1+exp(-((v+60)-15)/3));

Prr=[0. 0.6; 0.1 0.4]/3; 
% W = [1.25 -0.65; 1.2 -0.5];  
% W = [1.8 -1; 1.8 -0.5];  
% W=[1.25 -0.7; 1.2 -0.5]*0.8; 
% W=[0. -0.5; 1 -1]*5; 
%W=[0. -2.5; 5 -5];   %original
W = [0. -2.5; 5 -5]/2;   %/5
%W = [0. -5; 5 -5]/2;
J=W./Prr./[Ne Ni;Ne Ni]/param.tau_syn*1e3;    
J(1,1)=0; 

Ncluster=2; 
Nepop = 2040;   %2040;        %1200; 
Nipop = 360;    %360;         %300;
param.Nepop = Nepop; 
param.Nipop = Nipop; 
Nstim = length(stimID); %int32(0.1*N);    % 10 percent of N
N_cont = int32(0.25*N);     %no. of neurons getting I_cont
%stimulation
%stimstr = 2;  original

if strcmp(task,'randInput')
   % Nstim = Nepop;    
    T=3000;   % (ms)
    t_on=2000;
    t_off=t_on+100;
    %stimID = load('only_photostim/stimID.txt');
    %if trial ==1
        %if strcmp(modulation,'randBias')||strcmp(modulation,'Ibias')
            %load('data/2clusterLNP_randInput_ID1','stimID')
        %else
            %stimID = sort(randsample(N,Nstim),'asc');
            %disp('generate stimID')
        %end
    %else
        %load(filename(1),'stimID')
    %end
else
    %Nstim = Nepop;   
    %stimID = 1:Nstim;
    %stimID = load('dataE2/stimID.txt');
%     T=4000;   % (ms)
T = 1500; 
    t_on=1000;
    t_off=t_on+1000;  %sound applied from 2000 to 3000
end

param.T=T; 


% earlier was 0.6
switch task
    case 'constInput'
        Input=@(t) [stimstr*(t>=t_on&t<=t_off).*ones(Nstim,1); zeros(N-Nstim,1)] + param.bias; 
    
    case 'LeftSound'
        Input=@(t) [1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(N-2*Nepop,1)] + param.bias;
    
    case 'RightSound'
        Input=@(t) [0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(N-2*Nepop,1)] + param.bias;
    
    case 'randInput'
        Input=@(t) randinput(stimID,N,stimstr*(t>=t_on&t<=t_off).*ones(Nstim,1)) + param.bias;   
    
    case 'LeftSS'
        Input=@(t) [1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(N-2*Nepop,1)] + randinput(stimID,N,stimstr*(t>=t_on&t<=t_on+100).*ones(Nstim,1)) + param.bias; 
    
    case 'RightSS'
        Input=@(t) [0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(N-2*Nepop,1)] + randinput(stimID,N,stimstr*(t>=t_on&t<=t_on+100).*ones(Nstim,1)) + param.bias; 
    
    case 'LeftSoundI'  %left sound input to E1, E2, I1, I2
        Input=@(t) [1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(Ne - 2*Nepop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); zeros(Ni-2*Nipop,1)] + param.bias;

    case 'RightSoundI'  %right sound input to E1, E2, I1, I2
        Input=@(t) [0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(Ne - 2*Nepop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); zeros(Ni-2*Nipop,1)] + param.bias;

    case 'LeftSSI'  %photostim + left sound to E1, E2, I1, I2
        Input=@(t) [1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(Ne - 2*Nepop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); zeros(Ni-2*Nipop,1)] + randinput(stimID,N,stimstr*(t>=t_on&t<=t_on+100).*ones(Nstim,1)) + param.bias; 
    
    case 'RightSSI'  %photostim + right sound to E1, E2, I1, I2
        Input=@(t) [0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nepop,1); zeros(Ne - 2*Nepop,1); 0.85*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); 1*sound_str*(t>=t_on&t<=t_off).*(exp(-(t-t_on)/500)).*ones(Nipop,1); zeros(Ni-2*Nipop,1)] + randinput(stimID,N,stimstr*(t>=t_on&t<=t_on+100).*ones(Nstim,1)) + param.bias; 
   
end

dt=0.1; % (ms) 
param.Nskip=10;  %  save variables every Nskip time steps (Nskip=1 means save every time step) 
param.dt=dt; 

% param.idx_save=[randsample(Ne,80); Ne+randsample(Ni,20)]; 
% param.idx_save=[(1:5:Nstim)'; (Nstim+Nepop/2:Nepop:Ne)'];
%param.idx_save=[];

% param.idx_save=[randsample(Nstim,80); Nstim+randsample(Ne-Nstim,20)]; 
% param.Nsave=N; 
param.idx_save=[];   %1:N; 
param.Nsave=length(param.idx_save); 
maxns=N*T*.05; 
param.maxns=maxns; 

rng('shuffle')
    % initial condition
    Vm0=[-57 + randn(Ne,1); -57 + randn(Ni,1)];

%     tic

    [s1,Vm_save,Vm_ave,IsynE_save,IsynI_save,time_save]= ...
        LNP_ClusterNet3(Vm0,Input,PostID,out_degree,Wrr_c,param);
% 
%     elapsetime=toc;
%     fprintf('elapsetime=%.2f sec\n',elapsetime)
    
    s1(:,s1(2,:)==0)=[];
    save(filename(trial),'Input','s1','T','param','Ne','Ni','N','W_fname',...
        't_on','t_off','stimID','stimstr','sound_str')

end
end 

%{
figure
plot(s1(1,:),s1(2,:),'.')
%%
figure
stimID = sort(stimID,'asc'); 
nonstimID = setdiff(1:N,stimID); 
for ii =1:Nstim 
    ts = s1(1,(s1(2,:)==stimID(ii)));
    plot(ts, ii*ones(size(ts)),'r.')
    hold on 
end
for ii = 1:(N-Nstim)
    ts = s1(1,(s1(2,:)==nonstimID(ii)));
    plot(ts, (ii+Nstim)*ones(size(ts)),'b.')
end
%%
figure
RateE1=hist(s1(1,s1(2,:)<=Nepop&s1(2,:)>0),1:T)*1e3/Nepop;
RateE2=hist(s1(1,s1(2,:)<=2*Nepop&s1(2,:)>=Nepop+1),1:T)*1e3/Nepop;
% RateI=hist(s1(1,s1(2,:)>Ne&s1(2,:)<N),1:T)*1e3/Ni;
t_h = -50:50;
h = normpdf(t_h,0,10);
plot(imfilter(RateE1,h,'conv'))
hold on
plot(imfilter(RateE2,h,'conv'))

%%
rate_pre = hist(s1(2,s1(1,:)<t_on&s1(1,:)>200),1:N)*1e3/(t_on-200);
rate_post = hist(s1(2,s1(1,:)>t_on&s1(1,:)<t_off),1:N)*1e3/(t_off-t_on);
figure
subplot(1,2,1)
plot(rate_pre,rate_post,'.')
hold on 
plot(rate_pre(stimID),rate_post(stimID),'.')
plot(xlim,xlim,'k-') 

MI = (rate_post-rate_pre)./(rate_post+rate_pre); 
subplot(1,2,2)
hist(MI(1:Ne),30)
nanmean(abs(MI(1:Ne)))

%%
% id=1;
% figure
% 
% subplot(2,1,1)
% plot(time_save,IsynI_save(id,:))
% hold on 
% plot(time_save,imfilter(IsynI_save(id,:),ones(1,50)/50))
% title(['neuron ID: ' num2str(param.idx_save(id))])
% subplot(2,1,2)
% plot(time_save,Vm_save(id,:))
%}


