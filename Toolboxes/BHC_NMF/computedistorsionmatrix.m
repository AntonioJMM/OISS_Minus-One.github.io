function [dmatrix] = computedistorsionmatrix(B,X_ft,Y_ft)
% [dmatrix] = computedistorsionmatrix(B,X_ft,Y_ft)
% Funcion que calcula la matriz distorsion a partir de la ecuacion de la
% B-divergencia. Evita infinitos y devuelve valor absoluto (No negativos).
% Julio Carabias / Francisco Rodriguez -- December 2011

X_ft=X_ft+eps;
Y_ft=Y_ft+eps;

switch B
    case 0
        dmatrix = X_ft./Y_ft - log( X_ft./Y_ft ) - 1;        
    case 1
        dmatrix = X_ft .* (log(X_ft) - log(Y_ft)) + (Y_ft-X_ft);
    otherwise
        dmatrix = 1/(B*(B-1)) * ( X_ft.^B + (B-1)*Y_ft.^B - B*X_ft.*Y_ft.^(B-1) );        
end;

%pos_infinite=~isfinite(dmatrix);
%dmatrix(pos_infinite)=0;
%dmatrix(pos_infinite)=max(max(dmatrix));
%dmatrix=abs(dmatrix);

dmatrix(~isfinite(dmatrix))=0;

return;

