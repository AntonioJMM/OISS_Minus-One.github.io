%  THIS IS WORK IN PROGRESS. FURTHER TESTING IS NEEDED.
%
%  Description
% 
%  Compute the product of two multidimensional arrays in Matlab.
%  This is similar to TPROD* but it doesn't require a compiler.
%  tensor_mult also supports sparse matrices, which TPROD doesn't.
%
%  Sometimes tprod is faster, and sometimes tensor_mult is faster.
%  I haven't figured out why, yet.
% 
%  Note: If you need the result in a different order, use permute, see
%        example 4.
%
%  * https://de.mathworks.com/matlabcentral/fileexchange/16275-tprod-arbitary-tensor-products-between-n-d-arrays


%% ### Example 1 (Matrix-matrix product)
clear
% C_{ac} = A_{ab} B_{bc}
%             12     12
%              ^     ^
%
% sum over second index of A and first index of B
% This is equivalent to C = A * B
A = randn(5,7);
B = randn(7,5);
% tensor_mult command:

C_tm = tensor_mult(A, B, 2, 1);

% equivalent tprod command:

C_tp = tprod(A, [1 -1], B, [-1 2]);

%% ### Example 2 (Transposed matrix-matrix product)
clear
%  C_{bc} = A_{ab} B_{ac}
%              12     12
%              ^      ^
% 
%  sum over first index of A and first index of B
%  This is equivalent to C = A.' * B
A = randn(5,7);
B = randn(5,9);
% tensor_mult command:

C_tm = tensor_mult(A, B, 1, 1);

% equivalent tprod command:

C_tp = tprod(A, [-1 1], B, [-1 2]);

%% ### Example 3 (Tensor product without permutation)
clear
%  C_{cedf} = A_{abce} B_{dfab}
%                1234     1234
%                ab         ab
% 
%  sum over first  index of A and third  index of B and
%  sum over second index of A and fourth index of B
A = randn(5,6,7,8);
B = randn(9,10,5,6);
% tensor_mult command:

C_tm = tensor_mult(A, B, [1 2], [3 4]);

% equivalent tprod command:

C_tp = tprod(A, [-1 -2 1 2], B, [3 4 -1 -2]);
 
%% ### Example 4 (Tensor product with permutation)
clear
%  C_{decf} = A_{abce} B_{dfab}
%                1234     1234
%                ab         ab
% 
%  sum over first  index of A and third  index of B and
%  sum over second index of A and fourth index of B
%  reorder C_{cedf} -> C_{decf}
A = randn(5,6,7,8);
B = randn(9,10,5,6);
% tensor_mult command:

C_tm = permute(tensor_mult(A, B, [1 2], [3 4]), [3,2,1,4]);

% equivalent tprod command:

C_tp = tprod(A, [-1 -2 3 2], B, [1 4, -1 -2]);

%% ### Example 5 (Outer product)
clear
%  C_{abcdefg} = A_{abc} B_{defg}
% 
%  outer product
A = randn(4,5,6);
B = randn(7,8,9,10);
% tensor_mult command:

C_tm = tensor_mult(A, B, [], []);

% equivalent tprod command:

C_tp = tprod(A, [1 2 3], B, [4,5,6,7]);
