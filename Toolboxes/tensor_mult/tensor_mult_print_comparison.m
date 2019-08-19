function print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)
    fprintf('toc_tm   : %f\n', toc_tm );
    fprintf('toc_tp   : %f\n', toc_tp);
    fprintf('toc_for  : %f\n', toc_fo);
    fprintf('fo / tm  : %f\n', toc_fo / toc_tm );
    fprintf('fo / tp  : %f\n', toc_fo / toc_tp );
    fprintf('tp / tm  : %f\n', toc_tp / toc_tm );
    
    sizeC = numel(C_fo);

    fprintf('error tm : %d\n', norm(reshape(C_tm , sizeC, 1) - reshape(C_fo, sizeC, 1)));
    fprintf('error tp : %d\n\n', norm(reshape(C_tp , sizeC, 1) - reshape(C_fo, sizeC, 1)));
end