function rosePlotter(filelist,prefix,rgbtrip,plotTitle)
    fh=fopen(filelist,'r');
    nextLine=fgetl(fh);
    linecount=1;
    while(ischar(nextLine))
        [pardir, filename, ext] = fileparts(nextLine);
        if(exist([pardir filesep filename],'dir'))
            mastervelfile = [pardir filesep filename filesep 'PIV' filesep 'masterVels.mat'];
            if(exist(mastervelfile,'file'))
                load(mastervelfile);
                superTheta{linecount} = masterTheta;
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
    
    supersuperTheta = reshape(superTheta{1},[numel(superTheta{1}),1]);
    for i=2:numel(superTheta)        
        supersuperTheta = cat(1,supersuperTheta,reshape(superTheta{i},[numel(superTheta{i}),1]));
        
    end
    
    [tout,rout]=rose(supersuperTheta);
    p=polar(0,600000);
    set(p,'Visible','off');
    hold on
    p=polar(tout,rout);
    hold off

    set(p, 'Linewidth',1,'Color','k');
    
    xvals = get(p,'XData');
    yvals = get(p,'YData');
    xzeros = xvals==0;
    nxvals=xvals;
    nyvals=yvals;
    nxvals(xzeros)=[];
    nyvals(xzeros)=[];
    set(gca, 'nextplot' ,'add');
    plot(nxvals,nyvals);
    g=patch(nxvals,nyvals,rgbtrip);
    set(g,'facealpha',1);
    plot(xvals,yvals,'Linewidth',1,'Color','k');
    %title(plotTitle);

    set(findall(gca,'String','  50000'),'String',' ')
    set(findall(gca,'String','  100000'),'String',' ')
    set(findall(gca,'String','  150000'),'String',' ')
    set(findall(gca,'String','  200000'),'String',' ')
    set(findall(gca,'String','  300000'),'String',' ')
    set(findall(gca,'String','  400000'),'String',' ')
    set(findall(gca,'String','  500000'),'String',' ')
    set(findall(gca,'String','  600000'),'String',' ')
    set(findall(gca,'String','  700000'),'String',' ')
    set(findall(gca,'String','  800000'),'String',' ')
    set(findall(gca,'String','  900000'),'String',' ')
    set(findall(gca,'String','  1000000'),'String',' ')
    
    figfile = [pardir filesep prefix '-rose.ps'];
    saveas(gcf,figfile)
    
    
    %txt=findall(gcf,'type','text');
    % delete the text objects
    % delete(txt);
    
end