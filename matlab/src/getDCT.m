% compute DCT

function c_mat = getDCT(temp, pp, D)
    nSteps = fix(size(temp, 1) / pp);
    c_mat = zeros(pp, nSteps, size(temp, 2));

    for j = 1:size(temp, 2)
        for i = 1:nSteps
            c_mat(:, i, j) = D' * temp(pp *(i-1) + 1:(pp * i), j);
        end
    end
end