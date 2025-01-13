% function x = myfunction(x1, x2, lambda)
%
% Purpose:
% Simple function performing a linear combination
%
% Inputs:
% value1, value2, weight of the first value \in (0, 1)
%
% Outputs:
% Result of linear combination
% -------------------------------------------------------------------------
function x = myfunction(x1, x2, lambda)

    if (lambda < 0) | (lambda > 1)
        fprintf('ERROR: lambda not in range') 
        x = 'error';      
    else 
        x = x1*(lambda) + x2*(1-lambda);
    end     

end

