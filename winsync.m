function [syncIndex, corrValue] = winsync(data, ref, refLen, corr)

    % in_data & rf_data must be [N * 1] format;
    data2 = abs(reshape(data, [], 1)) .^ 2;
    ref2 = abs(reshape(ref, [], 1)) .^ 2;
    s = ref2(1 : refLen);
    R = [];
    s = s - sum(s) / refLen;
    stda = sqrt(sum(s .^ 2) / refLen);
    
    for k = 1: length(data2) / 2 - refLen + 1
        c = data2(k : k + refLen - 1);
        c = c - sum(c) / refLen;
        stdc = sqrt(sum(c .^ 2) / refLen);
        P(k) = transpose(s) * c ./ refLen;
        R(k) = P(k) / (stda * stdc);
        if abs(R(k)) > corr
            break
        end
    end
    
    corrValue = R(k);
    syncIndex = k;

end

% The method used to jude the state of synchronization is to compare the
% received signal and the refSignal. If the comparison result is satisfied
% the condition 'corr', then break the function and return the corrValue
% and the syncIndex.
