function [A_kt]=performDTW_function(distorsionXstate,param)

states_seq = param.states_seq;
states_time = param.states_time;

C = [1,1,1; 1,2,2; 1,3,3; 1,4,4;
    2,1,1;
    3,1,1;
    4,1,1];

nstates = length(states_seq);
r_ntramas = size(distorsionXstate,2);
m_ntramas = states_time(2,end);

matrizTmxTr = zeros(m_ntramas,r_ntramas);
for ss_iter=1:nstates
    current_state = states_seq(ss_iter);
    cstate_t_ini =  states_time(1,ss_iter);
    cstate_t_fin =  states_time(2,ss_iter);
    cstate_miditramas = cstate_t_fin - cstate_t_ini + 1;
    matrizTmxTr(cstate_t_ini:cstate_t_fin,:)=repmat(distorsionXstate(current_state,:),cstate_miditramas,1);
end

matrizTmxTr = matrizTmxTr ./ repmat(sqrt(sum(matrizTmxTr.^2,1)),m_ntramas,1);
matrizTmxTr(~isfinite(matrizTmxTr)) = 0;

p_final = dpfast_var_offline_min(matrizTmxTr,C,0,1);


structMidi = activeKMidi(states_time,states_seq,param);

A_kt = zeros(size(distorsionXstate,1),size(distorsionXstate,2));
for ii=1:length(p_final)
    A_kt(structMidi{p_final(ii)},ii) = 1;
end

end