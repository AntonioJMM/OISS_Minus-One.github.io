
%% Introduction
% Here are some examples of inner and outer products.
% A comparison is made between using a tensor_mult, tprod* and using for 
% loops.
%
% * https://de.mathworks.com/matlabcentral/fileexchange/16275-tprod-arbitary-tensor-products-between-n-d-arrays

%% $C_{j,l} = \sum_{i=1}^I \sum_{k=1}^K \delta_{i,k} A_{i,j} B_{k,l} = \sum_{i=1}^I A_{i,j} B_{i,l}$
clear
I = 700; J = 500;
K = I; L = J;

A = randn(I,J);
B = randn(K,L);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, 1, 1);
toc_tm  = toc(tic_tm );

tic_tp = tic;
C_tp   = tprod(A, [-1 2], B, [-1 3]);
toc_tp = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(J,L);
for j = 1:J
    for l = 1:L
        for i = 1:I
            C_fo(j,l) = C_fo(j,l) + A(i,j) * B(i,l);
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)

%% $C_{j,k,m,n} = \sum_{i=1}^{I} \sum_{l=1}^{L} \delta_{i,l} A_{i,j,k} B_{l,m,n} = \sum_{i=1}^{I} A_{i,j,k} B_{i,m,n}$
clear
I = 20; J = 30; K = 40;
L = I ; M = 50; N = 60;

A = randn(I,J,K);
B = randn(L,M,N);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, 1, 1);
toc_tm  = toc(tic_tm );

tic_tp  = tic;
C_tp    = tprod(A, [-1 2 3], B, [-1 4 5]);
toc_tp  = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(J,K,M,N);
for j = 1:J
    for k = 1:K
        for m = 1:M
            for n = 1:N
                for i = 1:I
                    C_fo(j,k,m,n) = C_fo(j,k,m,n) + A(i,j,k) * B(i,m,n);
                end
            end
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)

%% $C_{i,k,l,n} = \sum_{j=1}^{J} \sum_{m=1}^{M} \delta_{j,m} A_{i,j,k} B_{l,m,n} = \sum_{j=1}^{J} A_{i,j,k} B_{l,j,n}$
clear
I = 20; J = 30; K = 40;
L = 50; M = J ; N = 60;

A = randn(I,J,K);
B = randn(L,M,N);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, 2, 2);
toc_tm  = toc(tic_tm );

tic_tp  = tic;
C_tp    = squeeze(tprod(A, [1 -1 3], B, [4 -1 6]));
toc_tp  = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(I,K,L,N);
for i = 1:I
    for k = 1:K
        for l = 1:L
            for n = 1:N
                for j = 1:J
                    C_fo(i,k,l,n) = C_fo(i,k,l,n) + A(i,j,k) * B(l,j,n);
                end
            end
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)

%% $C_{i,k,l,m,n,p,q} = \sum_{j=1}^{J} \sum_{o=1}^{O} \delta_{j,o} A_{i,j,k,l} B_{m,n,o,p,q} = \sum_{j=1}^{J} A_{i,j,k,l} B_{m,n,j,p,q}$
clear
I = 5; J = 6; K = 7; L = 8;
M = 9; N = 10; O = J ; P = 11; Q = 12;

A = randn(I,J,K,L);
B = randn(M,N,O,P,Q);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, 2, 3);
toc_tm  = toc(tic_tm );

tic_tp  = tic;
C_tp    = tprod(A, [1 -1 2 3 4], B, [5 6 -1 7 8]);
toc_tp  = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(I,K,L,M,N,P,Q);
for i = 1:I
    for k = 1:K
        for l = 1:L
            for m = 1:M
                for n = 1:N
                    for p = 1:P
                        for q = 1:Q
                            for j = 1:J
                                C_fo(i,k,l,m,n,p,q) = C_fo(i,k,l,m,n,p,q) + A(i,j,k,l) * B(m,n,j,p,q);
                            end
                        end
                    end
                end
            end
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)

%% $C_{i,j,k} = A_{i,j} B_{k}$
clear
I = 100; J = 150;
K = 200;

A = randn(I,J);
B = randn(K,1);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, [], []);
toc_tm  = toc(tic_tm );

tic_tp  = tic;
C_tp    = tprod(A, [1 2], B, 3);
toc_tp  = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(I,J,K);
for i = 1:I
    for j = 1:J
        for k = 1:K
            C_fo(i,j,k) = A(i,j) * B(k);
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)

%% $C_{i,j,k,l,m,n} = A_{i,j,k} B_{l,m,n}$
clear
I = 10; J = 15; K = 20;
L = 25; M = 30; N = 35;

A = randn(I,J,K);
B = randn(L,M,N);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, [], []);
toc_tm  = toc(tic_tm );

tic_tp  = tic;
C_tp    = tprod(A, [1 2 3], B, [4 5 6]);
toc_tp  = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(I,J,K,L,M,N);
for i = 1:I
    for j = 1:J
        for k = 1:K
            for l = 1:L
                for m = 1:M
                    for n = 1:N
                        C_fo(i,j,k,l,m,n) = A(i,j,k) * B(l,m,n);
                    end
                end
            end
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)

%% $C_{i,j,k,l,m,n} = A_{i,j,k} B_{l,m,n}$
% Sometimes the for loop version is faster
clear
I = 2; J = 2; K = 2;
L = 2; M = 2; N = 2;

A = randn(I,J,K);
B = randn(L,M,N);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, [], []);
toc_tm  = toc(tic_tm );

tic_tp  = tic;
C_tp    = tprod(A, [1 2 3], B, [4 5 6]);
toc_tp  = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(I,J,K,L,M,N);
for i = 1:I
    for j = 1:J
        for k = 1:K
            for l = 1:L
                for m = 1:M
                    for n = 1:N
                        C_fo(i,j,k,l,m,n) = A(i,j,k) * B(l,m,n);
                    end
                end
            end
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)

%% $C_{k,l,m,n,q} = \sum_{j=1}^{J} \sum_{o=1}^{O} \sum_{i=1}^{I} \sum_{p=1}^{P} \delta_{j,o} \delta_{i,p} A_{i,j,k,l} B_{m,n,o,p,q} = \sum_{j=1}^{J} \sum_{i=1}^{I} A_{i,j,k,l} B_{m,n,j,i,q}$
clear
I = 10; J = 12; K = 14; L = 16;
M = 18; N = 10; O = J ; P = I; Q = 12;

A = randn(I,J,K,L);
B = randn(M,N,O,P,Q);

tic_tm  = tic;
C_tm    = tensor_mult(A, B, [2,1], [3,4]);
toc_tm  = toc(tic_tm );

tic_tp  = tic;
C_tp    = tprod(A, [-1 -2 1 2], B, [3 4 -2 -1 5]);
toc_tp  = toc(tic_tp);

tic_fo = tic;
C_fo = zeros(K,L,M,N,Q);
for k = 1:K
    for l = 1:L
        for m = 1:M
            for n = 1:N
                for q = 1:Q
                    for i = 1:I
                        for j = 1:J
                            C_fo(k,l,m,n,q) = C_fo(k,l,m,n,q) + A(i,j,k,l) * B(m,n,j,i,q);
                        end
                    end
                end
            end
        end
    end
end
toc_fo = toc(tic_fo);

tensor_mult_print_comparison(toc_fo, toc_tm, toc_tp, C_fo, C_tm, C_tp)
