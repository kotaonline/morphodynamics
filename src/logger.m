function logger(logFile, logMode, logData)
    persistent FID
    % Open the file
    if strcmp(logMode, 'open')
        FID = fopen(logFile, 'w');
        if FID < 0
            error('Cannot open file');
        end
    elseif strcmp(logMode, 'close')
        fclose(FID);
        FID = -1;
    elseif strcmp(logMode, 'write')
        fprintf(FID, '%s\n', logData);
    else
        error('Invalid log mode');
    end
end
% Write to the screen at the same time:
% fprintf('%s\n', Data);