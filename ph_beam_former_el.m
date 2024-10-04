function [angles] = ph_beam_former_el(samples)
    cols = 142;
    rows = 37;
    for j = 1:rows
        for k = 1:cols
            X = samples;
            R = (X * X')/size(X,2);
            R_inv = pinv(R); 
            theta_scan = linspace(-pi/12, pi/12, 30); % Region beam is covering
            results = zeros(30,2);

            for i = 1:length(theta_scan)
                theta = theta_scan(i);
                phi = pi * sin(theta); % Using sin not sind because theta_scan in radians
                s = [1; exp(-1j*phi)];
                w = (R_inv * s)/(s' * R_inv * s);
                X_weighted = w' * X;
                power_dB = 10 * log10(var(X_weighted));
                results(i,1) = power_dB;
                results(i,2) = theta;
            end
            results(:,1) = results(:,1) - max(results(:,1));
            sorted = sortrows(results, 1, "descend");
        end
    end
    angles = squeeze(sorted(1:10,2)) * 180/pi;

end