% Sherman Morrison
% Inputs:   O = current target, C= current inverse
%           shift, idx: amount of shift, index
% Outputs:  O = new target, C= new inverse

function [O, C] = update_sherman_morrison_diag(O,C,shift,idx)

O(idx,idx) = O(idx,idx) + shift;
c_d = C(idx,idx);
C = C - ((C(:,idx) * shift)*(C(idx,:)))/(1+shift*c_d);

end