function Y_corrected = faces_baseline_correction(Y);
    % Y is a n_electrodes x n_frames matrix
    % Return Y from which the average voltage in the first 100 ms is
    % subtracted
    baseline_end_frame = floor(size(Y,2)/2048*100);
    baseline = mean(Y(:,1:baseline_end_frame,:),2);
    Y_corrected = bsxfun(@minus,Y,baseline);
end

