function p = parallel_waitbar_update(n, wbar)
    persistent total count wb
    if nargin == 2
        % Initialization mode
        wb = wbar;
        total = n;
        count = 0;
    else
        % Increment count for afterEach call
        count = count + 1;
        p = count / total;
        percentage = p * 100; % 
        waitbar(p, wb, sprintf('Progress: %.1f%%', percentage)); 
    end
end
