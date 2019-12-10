function getCorrelationLength(filelist,figtitle)
    fh=fopen(filelist,'r');
    nextLine=fgetl(fh);
    linecount=1;
    while(ischar(nextLine))
        [pardir, filename, ext] = fileparts(nextLine);
        if(exist([pardir filesep filename],'dir'))
            mastercorfile = [pardir filesep filename filesep 'Postprocess/DistCorr' filesep 'Spatial_Correlation.mat'];
            if(exist(mastercorfile,'file'))
                load(mastercorfile);
                superCorr{linecount} = mean(spCorr(:,144:168),2);
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
    
    
    corvals=cell2mat(superCorr);
    
    xvals=linspace(1,size(corvals,1),size(corvals,1));
    dvals=all_dists.*(16*1.29);
    plot(dvals,corvals); hold on
    plot(linspace(1,400,400),zeros(400,1),'-','LineWidth',1)
    title(figtitle);
    xlim([0,400]);
    cv=corvals;
    cv(cv<=0.01)=0;
    cv(isnan(cv))=0;
    cvm=cv<=0.01;
    for i=1:size(cvm,2)
        clen(i)=dvals(find(cvm(:,i),1,'first'));
    end
    clen
    median(clen)
end