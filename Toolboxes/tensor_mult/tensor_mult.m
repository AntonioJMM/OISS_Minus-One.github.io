%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2017 Thomas Leitz (thomas.leitz@web.de, Unruhschuh.com)
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the product of two multidimensional arrays.
% Only tensor_mult.m is needed.
%
% Note: If you need the result in a different order, use permute, see
%       example 4.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1
%
% C_{ac} = A_{ab} B_{bc}
%
% sum over second index of A and first index of B
%
% C = tensor_mult(A, B, 2, 1);
%
% This is equivalent to C = A * B
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2
%
% C_{bc} = A_{ab} B_{ac}
%
% sum over first index of A and first index of B
%
% C = tensor_mult(A, B, 1, 1);
%
% This is equivalent to C = A.' * B
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3
%
% C_{cedf} = A_{abce} B_{dfab}
%               1234     1234
%               ^^         ^^
%
% sum over first  index of A and third  index of B and
% sum over second index of A and fourth index of B
%
% C = tensor_mult(A, B, [1 2], [3 4]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4
%
% C_{decf} = A_{abce} B_{dfab}
%               1234     1234
%               ^^         ^^
%
% sum over first  index of A and third  index of B and
% sum over second index of A and fourth index of B
% reorder C_{cedf} -> C_{decf}
%
% C = permute(tensor_mult(A, B, [1 2], [3 4]), [3,2,1,4]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 5
%
% C_{cedf} = A_{cbea} B_{dfab}
%               1234     1234
%                ^ ^       ^^
% sum over the fourth index of A and third  index of B and
% sum over the second index of A and fourth index of B
%
% C = tensor_mult(A, B, [4 2], [3 4]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6
%
% C_{abcdefg} = A_{abc} B_{defg}
%
% outer product
%
% C = tensor_mult(A, B, [], []);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = tensor_mult(A, B, sum_idx_A, sum_idx_B)

    sum_idx_A = reshape(sum_idx_A, 1, numel(sum_idx_A));
    sum_idx_B = reshape(sum_idx_B, 1, numel(sum_idx_B));

    size_A = size(A);
    size_B = size(B);
    
    num_size_A = numel(size_A);
    num_size_B = numel(size_B);
    
    perm_A = 1:num_size_A;
    perm_B = 1:num_size_B;
    
    num_idx = size(sum_idx_A, 2);
    
    assert(isequal(size(sum_idx_A) , size(sum_idx_B) ));
    assert(isequal(size_A(sum_idx_A), size_B(sum_idx_B)));
    
    sum_dim = prod(size_A(sum_idx_A));

    for i = 1:num_idx
        perm_A = perm_A(perm_A~=sum_idx_A(i));
        perm_B = perm_B(perm_B~=sum_idx_B(i));
    end
    
    perm_A = [perm_A, sum_idx_A];
    perm_B = [sum_idx_B, perm_B];
    
    size_A(sum_idx_A) = [];
    size_B(sum_idx_B) = [];
    
    if isempty(size_A)
        size_A = 1;
    end
    if isempty(size_B)
        size_B = 1;
    end
    
    C = squeeze(reshape(...
        reshape(permute(A, perm_A), [prod(size_A), sum_dim]) * ...
        reshape(permute(B, perm_B), [sum_dim, prod(size_B)]), ...
        [size_A,size_B]));

end
