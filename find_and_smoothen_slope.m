function [slope, smooth_slope] = find_slope(filename)
    % find slope of data in `filename`; also smoothen it
    
    % assumes first row contains header
    M = csvread(filename, 1);
    
    % first column has distances and second column has corresponding
    % elevation
    
    M2 = M(1:end - 1, :);
    M3 = M(2:end, :);
    
    delM = M3 - M2;
    
    % del elevation / del distance
    slope = delM(:, 2) ./ delM(:, 1);
    smooth_slope = smoothdata(delM(:, 2), 'gaussian', 20) ./ smoothdata(delM(:, 1),'gaussian', 20);
    
    % smoothen
    % M(:, 1) = smoothdata(M(:, 1), 'loess');
    % M(:, 2) = smoothdata(M(:, 2), 'loess');
    
    % comment out if necessary
    % elevation vs distance
    plot(smoothdata(M(:, 1), 'loess'), smoothdata(M(:, 2), 'loess'), 'LineWidth', 1.5)
    
    
    
end