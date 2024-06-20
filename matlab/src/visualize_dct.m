% function to visualize DCT matrix
function c = visualize_dct(c)
    for i = 1:size(c, 1)
        c(i,:) = 255 * (c(i,:) - min(c(i,:))) / (max(c(i,:)) - min(c(i,:)));
    end
end