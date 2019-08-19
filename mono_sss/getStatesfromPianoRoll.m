function [E_kpj,S_tmk,states_seq,states_time,k_max] = getStatesfromPianoRoll(NMFparams,R_ptj)
% [E_kp,S_tmk,states_seq,states_time,k_max]=getStatesfromPianoRoll(NMFparams,R_ptj,)
% Generate the states from the MIDI piano roll
%
% Inputs:
% NMFparams - NMF default params
% R_ptj - MIDI transcription
%
% Outputs:
% E_kp              - Matrix relation States x MIDI notes
% S_tmk             - MIDI time x States matrix
% states_seq        - States time sequence
% states_time       - States t_ini:t:fin (from MIDI transcription)
% k_max             - Number of states
%
% Load Parameters
t_max = NMFparams.t_max;
p_max = NMFparams.p_max / 4;
j_max = NMFparams.j_max;
%%%%%%%

% Current capacity (Optimization)
BLOCK_SIZE = 200;
states_Size = BLOCK_SIZE;
seq_Size = BLOCK_SIZE;

% Generate states output parameters
E_kpj = zeros(states_Size,p_max,j_max);
S_tmk = zeros(t_max,states_Size);

% Generate sequence output parameters
states_seq = zeros(1,seq_Size);
states_time = zeros(2,seq_Size);

% Initialize firt state to silence
states_seq(1)=1; % State number
states_time(1,1)=1; % t_ini
states_time(2,1)=1; % t_fin

% Pointers to last free position
states_Ptr = 2;
seq_Ptr = 2;


%% FUNCTION BODY
for t=1:t_max
    if(t == 55)
    end
    cstates_num = find(R_ptj(:,t,:)>0); % Active notes in t
    c_states_bin = zeros(p_max,j_max);    % Init binary notasXstate
    
    % Silence condition
    if ~isempty(cstates_num)
        c_states_bin(cstates_num) = 1; % No silencio, pasa notas a binario
    else
        cstates_num = 0; % Silence -> No notes
    end
    
    % Check if the state already exist
    membership=0;
    for ii=1:states_Ptr-1
        membership(ii) = isequal(squeeze(E_kpj(ii,:,:)),c_states_bin);
    end
    %     membership = ismember(E_kpj(1:states_Ptr-1,:,:),c_states_bin,'rows');
    if ~membership
        % New state
        E_kpj(states_Ptr,:,:) = c_states_bin;   % Save new state
        seq_pos = states_Ptr;      % Add new state to the sequence
        states_Ptr = states_Ptr + 1;         % Increment position pointer
    else
        % If the state already exist, check get the position
        seq_pos = find(membership>0);
    end
    
    % Check states sequence
    if seq_pos~=states_seq(seq_Ptr-1)
        % If change state -> save it in the sequence
        states_seq(seq_Ptr) = seq_pos;
        states_time(:,seq_Ptr) = [t ; t]; % t_ini;t_fin
        
        seq_Ptr = seq_Ptr + 1;  % increment position pointer
    else
        % If no new state -> increase end time
        states_time(2,seq_Ptr-1) = t;  % OJO: t_fin = t_fin+1
    end
    
    % Update MIDItimeXstates
    S_tmk(t,seq_pos) = 1;
    
    % ADD new block of memory if needed
    if( states_Ptr+(BLOCK_SIZE/10) > states_Size )  % less than 10%*BLOCK_SIZE free slots
        states_Size = states_Size + BLOCK_SIZE;     % add new BLOCK_SIZE slots
        E_kpj(states_Ptr+1:states_Size,:,:) = 0;
        S_tmk(:,states_Ptr+1:states_Size) = 0;
    end
    
    % ADD new block of memory if needed
    if( seq_Ptr+(BLOCK_SIZE/10) > seq_Size )  % less than 10%*BLOCK_SIZE free slots
        seq_Size = seq_Size + BLOCK_SIZE;     % add new BLOCK_SIZE slots
        states_seq(seq_Ptr+1:seq_Size) = 0;
        states_time(:,seq_Ptr+1:seq_Size) = 0;
    end
end

k_max = states_Ptr - 1; % No max of states

% REMOVE unused slots
E_kpj(states_Ptr:end,:,:) = [];
S_tmk(:,states_Ptr:end) = [];
states_seq(seq_Ptr:end) = [];
states_time(:,seq_Ptr:end) = [];

return;