function [AIC, SIC, MSC] = model_select(craw, cfit, numpar, wgt)

% number of frames
numfrm = size(craw,1);

% number of voxels
numvox = size(craw,2);

% weight
if nargin<4 | isempty(wgt)
    wgt = ones(numfrm,numvox);
end

% sum-of-squares
ss = sum(wgt.*(craw-cfit).^2,1);

% AIC
numpar_new = numpar + 1;
AIC = numfrm*log(ss/numfrm) + 2*numpar_new;
if numfrm/numpar<40
    AIC = AIC + 2*numpar_new*(numpar_new+1)/(numfrm-numpar_new-1);
end

if nargout>1
    % MSC (model selection criterion)
    sm = sum(wgt.*(craw-repmat(mean(craw,1),[numfrm 1])).^2,1);
    MSC = log(sm./ss) - 2*numpar/numfrm;

    % BIC (Schwartz criterion or Bayesian information criterion)
    SIC = numfrm*log(ss/numfrm) + numpar*log(numfrm);
end
