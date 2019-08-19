function structMidi = activeKMidi(states_time,states_seq,param)

nmidi = max(states_time(2,:));
t_delay = floor(1/param.t_salto);
structMidi = cell(nmidi,1);

for ii = 1:nmidi
    n = 1;
    ini = max(ii-t_delay,1);
    fin = min(ii+t_delay,nmidi);

    for jj=ini:fin
        k=find(states_time(1,:)<=jj);
        
        if(sum(states_seq(k(end)) == structMidi{ii}(:))==0) || isempty(structMidi{ii}(:))
            structMidi{ii}(n) = states_seq(k(end));
            n = n+1;
        end
    end
end
end