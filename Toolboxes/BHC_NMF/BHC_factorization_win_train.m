function [NMFparams,X_ft,Y_ft,A_ptj,C_pmj,S_fpj,distorsion]=BHC_factorization_win_train(NMFparams,NMF_inic,cfreq_amplitudes,param,barra)
% [NMFparams,X_ft,Y_ft,A_ptj,C_pmj,S_fpj,distorsion]=BHC_factorization_win(NMFparams,NMF_inic,cfreq_amplitudes,param,barra)
% Julio Carabias / Francisco Rodriguez -- March 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NMFparams,X_ft]=BHC_setParams(cfreq_amplitudes,param,NMFparams);

NUM_MAX_ITER = NMFparams.NUM_MAX_ITER;   % Num Iter
B = NMFparams.B;

f_max = NMFparams.f_max;
t_max = NMFparams.t_max;
p_max = NMFparams.p_max;
j_max = NMFparams.j_max;
m_max = NMFparams.m_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inicializo A
if isfield(NMF_inic,'A_ptj'),
    A_ptj=NMF_inic.A_ptj;    
else    
    A_ptj=ones(p_max,t_max,j_max);
end;

%% Inicializo C
if isfield(NMF_inic,'C_pmj'),
    C_pmj=NMF_inic.C_pmj;
else
    C_pmj=ones(p_max,m_max,j_max) + rand(p_max,m_max,j_max);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distorsion=zeros(1,NUM_MAX_ITER);

if barra
    % La puta barra
    h = waitbar(0,'Please wait...');
end;

%% Genero S_fpj Inicial    
[S_fpj]=BHC_generateSwin(NMFparams,C_pmj);     
%% Genero Y_ft Inicial
[Y_ft]=BHC_generateY(NMFparams,S_fpj,A_ptj);

%% Inicio el bucle
iter=1;
cond=1;
while cond==1 && iter<NUM_MAX_ITER,       
    if barra
        waitbar(iter/NUM_MAX_ITER,h,['Iteration ' num2str(iter) ' of ' num2str(NUM_MAX_ITER)]);
    end;
    
    %% Actualizo S_pjf
    [C_pmj]=BHC_updateCwin(NMFparams,C_pmj,A_ptj,X_ft,Y_ft);
    % Genero S_fpj    
    [S_fpj]=BHC_generateSwin(NMFparams,C_pmj);    
    % Normalizo S
    [S_fpj, C_pmj, A_ptj]=BHC_normalizeS(NMFparams,S_fpj,C_pmj,A_ptj);
    % Genero Y_ft
    [Y_ft]=BHC_generateY(NMFparams,S_fpj,A_ptj);
    
    %% Actualizo A_pjt   
    [A_ptj]=BHC_updateA(NMFparams,A_ptj,S_fpj,X_ft,Y_ft);       
    % Genero Y_ft
    [Y_ft]=BHC_generateY(NMFparams,S_fpj,A_ptj);
    
    %% Calculo la distorsion
    dmatrix = computedistorsionmatrix(NMFparams.B,X_ft,Y_ft);
    distorsion(iter) = sum(sum(dmatrix));
                               
    %% Actualizo la iteracion actual
    iter=iter+1;
    
end;

% Cierro la puta barra
if barra  
    close(h);
end;

return;        
