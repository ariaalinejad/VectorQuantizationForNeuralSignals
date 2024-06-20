% function to add normally distributed noise to a compressed signal

%use e.g.: filename = 'vq-lbg/hc1_cr-73_sndr-13.8_c';

function addNormalNoise(fileName, multiplier)
    arguments
        fileName = '';
        multiplier = 5;
    end
    %open file 
    base_folder = 'C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/mat';
    data = load(fullfile(base_folder, [fileName '.mat']));
    x = data.x;
    
    % add noise to signal
    rng(2,"twister"); % set random seed
    noise = randn(size(x))*multiplier;
    s = x + noise;

    base_folder = 'C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/bin';
    % Write the data to a binary file
    fileID = fopen(fullfile(base_folder, [fileName '_noisy' num2str(multiplier) '.bin']), 'wb');
    fwrite(fileID, s, 'double');
    fclose(fileID);

    sprintf('%s saved as Bin', fileName);
end