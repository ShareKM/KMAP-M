function [F, P] = fit_ftest(TAC,C1,n1,C2,n2,w)

% number of frames
numfrm = size(TAC,1);

% weighted residual sum of square
WRSS1 = w'*(C1-TAC).^2;
WRSS2 = w'*(C2-TAC).^2;

% Fisher Ratio
F = (WRSS1-WRSS2)/(n2-n1) ./ (WRSS2/(numfrm-n2));

% p value
if nargout>1
    P = 1 - fcdf(F, n2-n1, numfrm-n2);
end