% compute inverse DCT
function d_hat = getiDCT(temp, pp, D)
    nSteps = fix((size(temp, 1) * size(temp, 2)) / pp);
    d_hat = zeros(size(temp, 1) * size(temp, 2), size(temp, 3));

    for j = 1:size(temp, 3)
        for i = 1:nSteps
            d_hat((pp * (i-1)+1):(pp * i), j) = D * temp(:, i, j);
        end
    end
end