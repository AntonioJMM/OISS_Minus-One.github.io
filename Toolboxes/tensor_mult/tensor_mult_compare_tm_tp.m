%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Here are some test cases to compare tensor_mult and tprod*
%
% * https://de.mathworks.com/matlabcentral/fileexchange/16275-tprod-arbitary-tensor-products-between-n-d-arrays
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compare_tm_tp
    A = complex(randn(101,101,101), randn(101,101,101));
    B = complex(randn(101,101,101), randn(101,101,101));

    compare(A, B, [1 3], [2 1]);
    compare(real(A), real(B), [1 3], [2 1]);
    compare(A, B, [1 3 2], [2 1 3]);

    A = complex(randn(1001,101,101), randn(1001,101,101));
    B = complex(randn(1001,101), randn(1001,101));

    compare(A, B, 1, 1);
    compare(A, B, [1 2], [1 2]);

    A = complex(randn(51,101), randn(51,101));
    B = complex(randn(51,101), randn(51,101));
    
    compare(A, B, [], []);
    
    A = sparse(diag(randn(1000,1)));
    B = randn(1000,50,10);
    
    compare(A, B, 2, 1);
end

function compare(A, B, sum_idx_A, sum_idx_B)
    [sum_idx_tp_A, sum_idx_tp_B] = tm2tp(A, B, sum_idx_A, sum_idx_B);
    
    tic_tm = tic;
    C_tm = tensor_mult(A, B, sum_idx_A, sum_idx_B);
    toc_tm = toc(tic_tm);
    
    A_full = full(A);
    B_full = full(B);
    
    tic_tp = tic;
    C_tp = tprod(A_full, sum_idx_tp_A, B_full, sum_idx_tp_B);
    toc_tp = toc(tic_tp);
    
    error = norm(reshape(C_tm - C_tp, numel(C_tm), 1));
    
    fprintf('time tm: %f\n', toc_tm);
    fprintf('time tp: %f\n', toc_tp);
    fprintf('tp / tm: %f\n', toc_tp / toc_tm);
    fprintf('error  : %d\n\n', error);
end

function [sum_idx_tp_A, sum_idx_tp_B] = tm2tp(A, B, sum_idx_A, sum_idx_B)

    dim_A = numel(size(A));
    dim_B = numel(size(B));

    size_idx_tm = numel(sum_idx_A);
    
    sum_idx_tp_A = zeros(1, dim_A);
    sum_idx_tp_B = zeros(1, dim_B);
    
    for i = 1:size_idx_tm
        sum_idx_tp_A(sum_idx_A(i)) = -i;
        sum_idx_tp_B(sum_idx_B(i)) = -i;
    end
    
    k = 1;
    for i = 1:dim_A
        if sum_idx_tp_A(i) == 0
            sum_idx_tp_A(i) = k;
            k = k + 1;
        end
    end
    for i = 1:dim_B
        if sum_idx_tp_B(i) == 0
            sum_idx_tp_B(i) = k;
            k = k + 1;
        end
    end
    
end