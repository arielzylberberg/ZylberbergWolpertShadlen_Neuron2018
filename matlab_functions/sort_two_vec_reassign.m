function [Y] = sort_two_vec_reassign(A,B1,B2)
% function [Y] = sort_two_vec_reassign(A,B1,B2)
% sort A and B1; assigns Y the values of B2
% 
% i.e., the position of the largest value in A will have the value of
% B2 for which the value of B1 is maximal; same for 2nd largest, ...

if nargin<3 || isempty(B2)
    B2 = B1;
end

%sort indexes
[~,Ia] = sort(A);
[~,Ib] = sort(B1);

% unsort indexes
[~,J] = sort(Ia);
Y = B2(Ib(J)); 

