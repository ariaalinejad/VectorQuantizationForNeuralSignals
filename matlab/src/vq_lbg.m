% Started 22.05.24
function vq_lbg(M, N, k, N_ch, pp, n_remove, weighting, saveBin, saveAll)
    arguments
        M = 24000*2;%=48000, num samples
        N = 64;% num quataes
        k = 4; % num samples in time
        N_ch = 1; % num channels in space
        pp = 60; % DCT dimension
        n_remove = 0; % coeffs. to remove after DCT
        weighting = 0; % takes the value 0,1 or 2 (where 0 is SE) 
        saveBin = false; % save output, temp_hat, to bin file?
        saveAll = false; % save all plots and variables?
    end
    close all 
    addpath 'C:\Users\ariaa\OneDrive - NTNU\Masteroppgave\Code & Data\matlab\src'
    addpath 'C:\Users\ariaa\OneDrive - NTNU\Masteroppgave\Code & Data\matlab\notebooks - testing'
     
    tic; % for elapsed time
    
    
    % extract data
    base_folder = 'C:/Users/ariaa/OneDrive - NTNU/Masteroppgave/Code & Data/data/mat/';
    hc1_data = load(fullfile(base_folder, 'hc1_d533101_extracellular.mat'));
    hc1 = hc1_data.x';
    fs = 10000;
    fc_high = 4000;
    fc_low = 300;

    delta_eps = 50; % as defined in thesis
    n_bits = 16; % n bits used in uncompressed signal

    % choose weighting
    if weighting == 1
        w_power = 2;
    elseif weighting == 2
        w_power = round(log2(M/N));
    end
    
    
    % stoping criterion
    it_max = 500; % max num iterations


    %prepare time axis
    t_max = (M/2)/fs;
    stepsize = t_max/( (M/2) - 1 );
    t = 0:stepsize:t_max;

    %% DCT
    
    D = zeros(pp, pp);
    parfor p = 1:pp
        a = zeros(pp, 1);
        a(p) = 1;
        D(:, p) = idct(a);
    end
    
    
    %% Define filters
    
    low_10hz = fc_low / (fs / 2);
    high_10hz = fc_high / (fs / 2);
    [b_10hz, a_10hz] = butter(10, [low_10hz, high_10hz], 'bandpass');
    
    
    %% Filter:
    
    temp = hc1; %use first half as training set
    AP_train = zeros(size(temp));
    parfor i = 1:size(temp, 2)
        AP_train(:,i) =  filtfilt(b_10hz, a_10hz, temp(:,i));
    end
    
    AP_train = AP_train(1:M/2, 1:N_ch); % Use first half as training set
    
    c_mat = getDCT(AP_train, pp, D);
    
    % Visualize resulting coeff.s
    %{
    for i = 1:N_ch
        figure('Position', [0, 0, 1500, 250]);
        imagesc(visualize_dct(abs(c_mat(:,:,i))'));
        colorbar;
    end
    %}
    
    %% Remove coeff.
    
    if n_remove > 0
        c_mat_th = c_mat(1:end-n_remove, :, :);
    else
        c_mat_th = c_mat;
    end
    
    c_th = zeros(size(c_mat_th, 1) * size(c_mat_th, 2), size(c_mat_th, 3));
    for j = 1:size(c_mat_th, 3)
        for i = 1:size(c_mat_th, 2)
            c_th((i-1)*size(c_mat_th, 1) + 1 : i*size(c_mat_th, 1), j) = c_mat_th(:, i, j);
        end
    end
    
    
    %% VQ (code optimized with the help of copilot)
    
    % group coefficients
    temp = zeros(floor(size(c_th, 1) / k), k, N_ch);
    for j = 1:N_ch
        for i = 1:floor(size(c_th, 1) / k)
            temp(i, :, j) = c_th((i - 1) * k + 1 : i * k, j);
        end
    end
    
    % initialize parameters and variables
    l_1 = size(temp,1);
    j_min = zeros(l_1, 1);
    it = 0;
    L_max = N;
    y = randi([int32(min(temp,[], 'all')), int32(max(temp,[], 'all'))], 1, k, N_ch);
    L = 1;
    if weighting == 2 ||  weighting == 1
        w = temp.^w_power;
    else
        w = ones(size(temp,1),k,N_ch);%
    end
    
    % run loop until stopping criterion
    while true
        it = it + 1;
        y_old = y; 
    
        % assign each point the center that gives the lowest error
        parfor i = 1:l_1
            error = sum(w(i,:,:).*(temp(i,:,:) - y_old).^2, [2 3]);
            [~, j_min(i)] = min(error);
        end
        
        % set y as center of points in cluster j_min 
        parfor j = 1:size(y_old,1)
            if any(j_min == j)
                y_old(j,:,:) = sum(w(j_min == j,:,:).*temp(j_min == j,:,:), 1)./sum(w(j_min == j,:,:),1);
            else
                y_old(j,:,:) = temp(randi(l_1), :,:); %take the value of a random point
            end
        end
        
        % display delta y
        dy = sum(abs(y - y_old), 'all');
        disp(['Iteration: ', num2str(it)]);
        disp(['delta y: ', num2str(dy)]);

        % split center, stop, or go to next iteration
        if (L<L_max)
            y = zeros(L*2, k, N_ch);
            y(1:L,:,:) = y_old;
            parfor j = 1:size(y_old,1) % Here a modification has been made to the standard LGB because of computation speed
                eps = randi([-delta_eps, delta_eps], 1, k, N_ch);
                y(L+j, :, :) = y_old(j,:,:) + eps;
            end
            L = L*2;
        elseif( dy == 0 || it>=it_max)% set value of first code vector to 0 vector,  to allow zero value coefficients (making matrix sparse)
            y = y_old;
            y(1,:,:) = zeros(k, N_ch);
            break
        else
            y = y_old;
        end
    
    end
    
    % plot figures
    if k > 1
        figure;
        title([num2str(it), ' iterations centers (across samples)']);
        scatter(temp(:,1,1), temp(:,2,1), 'filled', 'DisplayName', 'x');
        hold on;
        scatter(y(:,1,1), y(:,2,1), 'DisplayName', 'y');
        xlabel('p=1');
        ylabel('p=2');
        legend('show');
        hold off;
    end   
    if N_ch > 1
        figure;
        title([num2str(it), ' iterations centers (across channels)']);
        scatter(temp(:,1,1), temp(:,1,2), 'filled');
        hold on;
        scatter(y(:,1,1), y(:,1,2));
        xlabel('coeff. 1')
        ylabel('coeff. 2')
        legend('show')
        hold off;
    end
    
    % Make a huffman test set
    parfor i = 1:l_1
        error =  mean(((temp(i,:,:) - y).^2), [2 3]);  %MSE
        [~, j_train(i)] = min(error); % argmin 
    end
    
   
    %% Evaluate results
    temp = hc1(M/2+1:M, 1:N_ch); % use second half as test set
    
    temp = reshape(temp, M/2, N_ch);
    
    %filter test set
    AP_test  = zeros(size(temp));
    parfor i = 1:size(temp, 2)
        AP_test(:, i)  = filtfilt(b_10hz, a_10hz, temp(:, i));
    end
    
    
    %% Compress test set
    
    % Perform DCT
    c_test_mat = getDCT(AP_test, pp, D);
    
    % Remove lowest coefficients
    c_test_mat_th = c_test_mat(1:end-n_remove,:,:);
    
    % get dimensions of matrix
    n_cmt_0 = size(c_test_mat_th, 1); % Number of elements in the first dimension of c_mat_th (cmt)
    n_cmt_1 = size(c_test_mat_th, 2);
    n_ct_0 = size(c_th, 1);
    
    % flatten coefficient matrix
    c_test_th = zeros(n_cmt_0 * n_cmt_1, N_ch);
    for j = 1:N_ch
        for i = 1:n_cmt_1
            c_test_th((n_cmt_0 * (i - 1)) + 1 : n_cmt_0 * i, j) = c_test_mat_th(:, i, j);
        end
    end
    
    % Resize to match VQ size
    c_test_vq = zeros(floor(n_ct_0 / k), k, N_ch);
    for j = 1:N_ch
        for i = 1:floor(n_ct_0 / k)
            c_test_vq(i, :, j) = c_test_th((i - 1) * k + 1 : i * k, j);
        end
    end
    
    
    % Find closest code vectors
    parfor i = 1:l_1
        error =  sum(((c_test_vq(i,:,:) - y).^2), [2 3]);  %SE
        [~, j_test(i)] = min(error); % argmin 
    end
    % y(j_test, :, :) are now the values, and j_test are the
    % indexes to be sent
    
    %% Reconstruct

    % Get coeff values
    yj = y(j_test,:,:);
    
    % get dimensions
    n_y_0 = size(yj, 1);
    n_ct_0 = size(c_th, 1);
    
    
    % Get Approximated coefficients
    c_re_th_ap = zeros(n_ct_0, N_ch); 
    for j = 1:N_ch
        for i = 1:n_y_0
            c_re_th_ap((i - 1) * k + 1 : i * k, j) = yj(i, :, j);
        end
    end
    
    % Turn into matrix
    c_re_mat_th_ap = zeros(n_cmt_0, n_cmt_1, N_ch);
    for j = 1:N_ch
        for i = 1:n_cmt_1
            c_re_mat_th_ap(:, i, j) = c_re_th_ap((i - 1) * n_cmt_0 + 1 : i * n_cmt_0, j);
        end
    end
    
    % Pad missing values
    c_re_mat_ap = padarray(c_re_mat_th_ap, [n_remove, 0, 0], 0, 'post'); % Reconstructed coefficients of test set
    
    c_test_mat_ap = padarray(c_test_mat_th, [n_remove, 0, 0], 0, 'post'); % Original coefficients of test set
    
    % Perform huffman coding
    len_huff = huff_VQ(j_test, j_train, N);
    len_original = log2(N)*size(j_test,2);
    cr_huff = len_original/len_huff;
    fprintf('CR from huffmanCoding: %.2f\n', cr_huff);
    
    % compute SNDR for coefficients
    % P_d = norm(temp_mat_ap(:));
    % P_error = norm(temp_mat_ap(:) - c_mat_ap(:));
    % SNDR = 20*log10(P_d/P_error);
    % fprintf('SNDR: %.2f dB\n', SNDR);
    % % fprintf('MSE: %.2f\n', mse);
    % fprintf('CR: %.2f\n', n_bits/log2(N)*k*N_ch*pp/(pp - n_remove)*cr_huff);
    
    % plot frequency of occurance for different coefficient calues
    figure;
    %title('Number of occurances for different coefficient values')
    histogram(reshape(c_re_th_ap, [size(c_re_th_ap,1)*size(c_re_th_ap,2), 1]), 101);
    xlabel('Coefficient Values');
    ylabel('Number of occurances');

    % plot original vs reconstucted coefficients (ch1)
    % figure;
    % plot(reshape(c_test_mat_ap.', 1, []), 'DisplayName', 'Original coeff.');
    % hold on;
    % plot(reshape(c_re_mat_ap.', 1, []), 'DisplayName', 'Reconstructed coeff.');
    % legend('show');
    
    %% Reconstruct
    
    % Compute the iDCT to reconstruct test set signal
    temp_hat = getiDCT(c_re_mat_ap, pp, D);
    
    % Create a copy of AP_test
    temp = AP_test;
    
    % Calculate the power of the original signal
    P_d = norm(temp(:));
    
    % Calculate the power of the error between temp and temp_hat
    P_error = norm(temp(:) - temp_hat(:));
    
    % Compute the Signal-to-Noise and Distortion Ratio (SNDR) in dB
    SNDR = 20 * log10(P_d / P_error);
    
    fprintf('SNDR: %.2f dB\n', SNDR);
    
    % compute MSE of error
    % mse = mean((temp(:) - temp_hat(:)).^2);
    % fprintf('MSE: %.2f\n', mse);
    
    % Compute the Compression Ratio (CR)
    CR = (n_bits*M/2) / len_huff; %
    fprintf('CR: %.2f\n', CR);
    
    % plot
    figure;
    for i=1:N_ch
        subplot(N_ch,1,i);
        plot(t', temp(:, i), 'DisplayName', 'Original');
        hold on;
        plot(t', temp_hat(:, i), 'DisplayName', 'Reconstructed');
        title(['Channel ' num2str(i)]);
        xlabel('seconds [s]');
        legend('show');
        hold off;
    end
    fprintf('Runtime: %.2f sec \n ', toc);
    %% Analyze spike
    spike_ch = 1; 
    
    % choose spike which spike to look at based on signal length, the spike
    % times are set through manual inspection
    if M==480000 
        spike_start = 5750;
        spike_end = 5780;
    elseif M==48000 
        spike_start = 588;
        spike_end = 602;
    elseif M==4800 
        spike_start = 1690;
        spike_end = 1720;
    elseif M==2400000 
        spike_start = 431; %spike at 436 (15 sample window)
        spike_end = 446; 
    elseif M==60000 
        spike_start = 2763;
        spike_end = 2778;
    elseif M==120000 
        spike_start = 259; % 2ms window
        spike_end = 281;
    elseif M==1200000 
        spike_start = 498; % 2ms window
        spike_end = 518;
    else
        spike_start = 1;
        spike_end = 30;
        disp('spike not defined')
    end
    
    % plot
    figure;
    plot(t(spike_start:spike_end)',     temp(spike_start:spike_end, spike_ch), 'DisplayName', 'Original');
    hold on
    plot(t(spike_start:spike_end)', temp_hat(spike_start:spike_end, spike_ch), 'DisplayName', 'Reconstructed');
    legend('show')
    xlabel('seconds [s]')
    title('Spike before and after reconstruction');
    hold off
    
    % Calculate the power of the original signal
    P_d = norm(temp(spike_start:spike_end, spike_ch));
    
    % Calculate the power of the error between temp and temp_hat
    P_error = norm(temp(spike_start:spike_end, spike_ch) - temp_hat(spike_start:spike_end, spike_ch));
    
    % Compute the Signal-to-Noise and Distortion Ratio (SNDR) in dB
    SNDR_spike = 20 * log10(P_d / P_error);
    fprintf('SNDR of spike: %.2f dB\n', SNDR_spike);
    
    
    %% Save reconstructed signal
    if saveBin
        x = temp_hat'; % notice the transpose
        save(fullfile(base_folder,['vq-lbg/hc1_cr-' num2str(round(CR)) '_sndr-' num2str(round(SNDR_spike,1)) '_c.mat']), 'x');
        matToBin(['vq-lbg/hc1_cr-' num2str(round(CR)) '_sndr-' num2str(round(SNDR_spike,1)) '_c']);
    end
    
    %% Save only parameters
    % save(fullfile(base_folder, ['vq-lbg/hc1_cr-' num2str(round(CR)) '_sndr-' num2str(round(SNDR,1)) '_parameters.mat']), "n_bits", "M", "N", "k", "pp", "N_ch", "n_remove", "delta_eps", "w_multiplier", "w_power");
    
    %% Save all variables & Figures
    if saveAll
        FolderName = ['hc1_cr-' num2str(round(CR)) '_sndr-' num2str(round(SNDR_spike,1)) '_M-' num2str(M)]; 
        mkdir([base_folder '/vq-lbg/' FolderName]);
    
        save(fullfile(base_folder, ['vq-lbg/' FolderName '/variables.mat'])); % save all variables
    
    
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        i = 0;
        for iFig = 1:length(FigList)
            i = i + 1;
            FigHandle = FigList(iFig);
            savefig(FigHandle, fullfile(base_folder, ['vq-lbg/' FolderName '/' int2str(i) '.fig']));
        end
        
        fid = fopen(fullfile(base_folder, ['vq-lbg/' FolderName '/output.txt']), 'wt');
        fprintf(fid,'SNDR of spike: %.2f dB\n', SNDR_spike);
        fprintf(fid,'SNDR of signal: %.2f dB\n', SNDR);
        fprintf(fid,'CR: %.2f\n', CR);           
        fprintf(fid,'CR_huff: %.2f\n', cr_huff);
        fprintf(fid,'Time %.2f s\n', toc);
        fprintf(fid,'Converged %d\n', ~isequal(it, it_max));
        fclose(fid);
    end
 
end