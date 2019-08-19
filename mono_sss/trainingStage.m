function [NMF_inic,NMFparams,param] = trainingStage(NMFparams,param)
% [NMF_inic,NMFparams,param] = trainingStage(NMFparams,param)
% Compute the training stage
%
% Input arguments:
%   NMFparams --> store the NMF parameters
%	param --> see userparameter.m to set initial values
%
% Output:
%   NMF_inic --> Initialized NMF parameters
%   NMFparams
%	param

%% Load Instrument Models from RWC Database
fprintf('\n\tLoading instrument models ... ');
C_pmj_aux = zeros(114,20,length(param.inst));
SV = [];
WV = [];

for ii = 1:length(param.inst)
    switch param.inst{ii}
        case 'vn'
            load([param.ruta_modelos '151VNNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
            SV(ii) = ii;
        case 'va'
            load([param.ruta_modelos '161VLNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
            SV(ii) = ii;
        case 'vc'
            load([param.ruta_modelos '171VCNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
            SV(ii) = ii;
        case 'db'
            load([param.ruta_modelos '181CBNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
            SV(ii) = ii;
        case 'fl'
            load([param.ruta_modelos '331FLNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
            WV(ii) = ii;
        case 'ob'
            load([param.ruta_modelos '291OBNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
        case 'cl'
            load([param.ruta_modelos '311CLNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
            WV(ii) = ii;
        case 'sax'
            load([param.ruta_modelos '271TSNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
        case 'bn'
            load([param.ruta_modelos '301FGNOM_modelB1_3.mat'],'C_pmj');
            C_pmj_aux(:,:,ii) = C_pmj;
            WV(ii) = ii;
    end
end
C_pmj = C_pmj_aux;
clear C_pmj_aux;

%% Extend the Models to 1/4 Semitone LOG. Freq Resolution
if param.midi_inc>1
    [C_pmj] = BHC_INCmodel(C_pmj,param.midi_inc);
end
fprintf('Done');

%% Inicialize Bases
fprintf('\n\tInitializing the bases ... ');
NMF_inic = [];
NMFparams = nmf_inc_setParamsINC(param,NMFparams);
NMFparams.j_max = length(param.inst);
[NMF_inic.S_fpj] = BHC_generateSwinINC(NMFparams,C_pmj);
fprintf('Done');

%% Synthetic WAV file processing
fprintf('\n\tReading synthetic WAV file ... ');
[sy,fs]=audioread([param.ruta_bbdd param.filename filesep 'Sco_' param.filename '.wav']);
sy=mean(sy,2)';
nsamples = length(sy);
fprintf('Done');

% Compute the spectrogram
fprintf('\n\tComputing the magnitude spectrogram ... ');
[param] = getParametrosMusica(param,nsamples,fs);
Xs_ft = computeCfreqFast(sy,param);
NMFparams.t_max = size(Xs_ft,2);
fprintf('Done');

%% Transcription from MIDI
fprintf('\n\tComputing transcription from MIDI ... ');
[GT] = load([param.ruta_bbdd param.filename filesep 'Sco_' param.filename '.txt']);

nsources = numel(unique(GT(:,3)));
inst = GT(:,3);
notes = GT(:,4);
onsets = round((GT(:,6))*fs / param.salto);
onsets = max(1,onsets);
offsets = round((GT(:,6)+GT(:,7))*fs / param.salto);
offsets = min(param.n_tramas, offsets);
real_clase = zeros(param.midi_max,param.n_tramas,nsources);
for n=1:length(notes)
    real_clase(notes(n),onsets(n):offsets(n),inst(n)) = 1;
end
real_clase = real_clase(param.midi_min:end,:,:);
fprintf('Done');

%% Define States
fprintf('\n\tGenerating states ... ')
[E_kpj,A_tk,param.states_seq,param.states_time,k_max] = getStatesfromPianoRoll(NMFparams,real_clase);
E_kpj(E_kpj>0) = 1;
NMFparams.k_max = k_max;

%% Extend the Models to 1/4 Semitone LOG. Freq Resolution
C_pij = zeros(NMFparams.p_max/4,NMFparams.p_max,NMFparams.j_max);

SV(SV==0) = [];
WV(WV==0) = [];

for jj=1:NMFparams.j_max
    if ismember(jj,SV) % INSTRUMENTS WITH STRONG VIBRATO
        % TAKE 1/2 SEMITONES ABOVE AND BELOW THE CURRENT NOTE
        C_pij(:,:,jj) = computeExtendFreq(param, 2);
    elseif ismember(jj,WV) % INSTRUMENTS WITH WEAK VIBRATO
        % TAKE 1/4 SEMITONES ABOVE AND BELOW THE CURRENT NOTE
        C_pij(:,:,jj) = computeExtendFreq(param, 1);
    else % INSTRUMENTS WITH NO VIBRATO
        % TAKE ONLY THE CURRENT SEMITONE
        C_pij(:,:,jj) = computeExtendFreq(param);
    end
end
E_kpj = tprod(C_pij,[-1 2 3],E_kpj,[1 -1 3]);
E_kpj(E_kpj>0) = 1;
fprintf('Done');

%% Learning States Bases
fprintf('\n\tLearning states bases from synthetic data ... ')
[NMF_inic.E_kpj] = nmf_inc_updateE(Xs_ft,NMFparams,NMF_inic.S_fpj,E_kpj,A_tk');
fprintf('Done');

return