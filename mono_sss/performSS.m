% Script to perform the Source Separation - Score Informed Approach
%
% Further info can be found in:
% "Online/Offline Score Informed Music Signal Decomposition: Application to
% Minus One"
% Antonio Jesus Munoz-Montoro, Julio Jose Carabias-Orti, Pedro
% Vera-Candeas, Francisco Jesus Canadas-Quesada and Nicolas Ruiz-Reyes
% EURASIP Journal on Audio, Speech, and Music Processing.
%
% Input arguments:
%	see userparameter.m
%
% Output:
%   separated signals stored in the path given by userparameters
%
% UJA 2019, Antonio Muñoz

fprintf('\t\t--- SOURCE SEPARATION IN ORCHESTRA ---');
%% Add routines to the path
addpath(genpath([pwd filesep '../Toolboxes']));
addpath(genpath([pwd filesep 'NMF_INC']));

%% Load parameters
load wind_excit; % Window Transform
NMFparams.wind_excit = wind_excit;

%% ANALIZE THE INPUT SIGNALS
dwav = dir(param.ruta_bbdd);
dwav(1:3) = [];
nwav = length(dwav);

for file = 1:nwav
    % Audio file detail
    filename = dwav(file).name;
    sep = strfind(filename,'_');
    inst = [];
    for ii = 1:length(sep)-2
        inst{ii} = filename(sep(ii+1)+1:sep(ii+2)-1);
    end
    inst{length(sep)-1} = filename(sep(end)+1:end);
    param.filename = filename;
    param.inst = inst;
    
    % Load Real Audio
    fprintf('\n\nStarting Real Signal Analysis for %s ... ',filename );
    fprintf('\n\tReading WAV file ... ');
    [x,fs] = audioread([param.ruta_bbdd filename filesep 'AuMix_' filename '.wav']);
    nsamples = length(x);
    fprintf('Done');
    
    % Compute the spectrogram
    fprintf('\n\tComputing the magnitude spectrogram ... ');
    [param] = getParametrosMusica(param,nsamples,fs);
    Xr_ft = computeCfreqFast(x,param);
    fprintf('Done');
    
    %% Training Stage
    [NMF_inic,NMFparams,param] = trainingStage(NMFparams,param);
    
    %% Alignment Stage
    [NMF_inic,NMFparams,param] = alignmentStage(Xr_ft,NMF_inic,NMFparams,param);
    
    %% Separation Stage
    % Estimating parameters
    fprintf('\n\tComputing SS ... ');
    [S_fpj,A_ptj]=NMF_factorization(Xr_ft,NMFparams,NMF_inic);
    fprintf('Done');
    
    % Perform Reconstruction - Wiener Filtering
    fprintf('\n\tPerforming Separation ... ');
    [senales]=senal_extraction_INC(x,NMFparams.fftparams,A_ptj,S_fpj);
    
    % Writing estimated signals
    for jj=1:NMFparams.j_max
        audiowrite([param.outpath filesep filename(1:sep(2)-1) '_' num2str(jj) '.wav'],senales(:,jj),fs);
    end
    
    clear A_ptj S_fpj;
    fprintf('Done');
end

fprintf('\nThanks for using me');
