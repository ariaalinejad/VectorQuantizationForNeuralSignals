% function to save .mat file as .bin file
function matToBin(fileName)

    base_folder = 'C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/mat/';
    data = load(fullfile(base_folder, [fileName '.mat']));
    x = data.x;
    
    base_folder = 'C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/bin/';
    % Write the data to a binary file
    fileID = fopen(fullfile(base_folder, [fileName '.bin']), 'wb');
    fwrite(fileID, x, 'double');
    fclose(fileID);

    sprintf('%s saved as Bin', fileName);
end
