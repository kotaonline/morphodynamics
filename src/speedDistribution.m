function speedDistribution(filelist, prefix, plotTitle)
    fh=fopen(filelist,'r');
    nextLine=fgetl(fh);
    linecount=1;
    while(ischar(nextLine))
        [pardir, filename, ext] = fileparts(nextLine);
        if(exist([pardir filesep filename],'dir'))
            mastervelfile = [pardir filesep filename filesep 'PIV' filesep 'masterVels.mat'];
            if(exist(mastervelfile,'file'))
                load(mastervelfile);
                superMag{linecount} = masterMags;
            else
                error('Something wrong! Entered a PIV directory without masterVels.mat\nI''m outta here!');
            end
        else
            error('Could not find a directory with the same name as input tif file. Did you forget to run morphodynamics?');
        end
        nextLine=fgetl(fh);
        linecount = linecount+1;
    end
    fclose(fh);
    
    speeds = cell2mat(superMag);
    speeds(speeds>65)=0;
    reshaped_speeds = reshape(speeds,[numel(speeds),1]);
    reshaped_speeds(isnan(reshaped_speeds))=[];
    
    speedfile = [pardir filesep prefix '-speeds.txt'];
    save(speedfile,'reshaped_speeds','-ascii','-tabs');
    
    hbins = 2.5:5:60;
    counts = hist(reshaped_speeds,hbins);
    total_counts = sum(counts);
    normalized_counts = counts./total_counts;
    
    histfile = [pardir filesep prefix '-hist.txt'];
    normhist = [hbins', normalized_counts'];
    
    save(histfile, 'normhist', '-ascii', '-tabs');
    
    for i=1:numel(superMag)
        nanless_speeds = superMag{i};
        nanless_speeds(isnan(nanless_speeds))=[];
        
        mean_sd{i} = [mean(nanless_speeds) std(nanless_speeds)];% numel(nanless_speeds)]
        
    end
    cell2mat(mean_sd')
        


end