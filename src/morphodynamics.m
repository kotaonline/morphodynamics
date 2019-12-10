function morphodynamics(inpFile,configFile)

    if(nargin~=2)
        error('Usage: morphodynamics([input Tif file],[configuration file])');
    end
    
    %% Add paths to relevant parts of MatPIV.

    addpath('PIV/MatPIV 1.7/src');
    addpath('PIV/MatPIV 1.7/filters');
    addpath('PIV/MatPIV 1.7/postprocessing');

    %% Initialize parameters -- populate params struct, verify input file format and generate output directories
    
    [params, dirs] = init(inpFile,configFile);
    % save final parameter and directory structures
    paramMatFile = [ dirs.expDir 'config.mat' ];
    
    %% Start logging if enabled
    
    if(params.log)
        logger(params.logfile, 'open', '');
    end
    
    %% Check for world coordinates

    if(exist(params.woco,'file'))
        if(params.log)
            logger(params.logfile, 'write', sprintf('Number of frames: %d', params.nFrames));
            logger(params.logfile, 'write', sprintf('Image size: %dx%d', params.width, params.height));
            logger(params.logfile, 'write', sprintf('Experiment name: %s', params.expName));
            logger(params.logfile, 'write', sprintf('Camera pixel size (um): %f', params.pixelSize));
            logger(params.logfile, 'write', sprintf('Pixels per micron: %f', params.pxperum));
            logger(params.logfile, 'write', sprintf('Interrogation window overlap: %f', params.overlap));
            logger(params.logfile, 'write', 'World coordinate file in place? ... yes');
            logger(params.logfile, 'write', ' ');
            logger(params.logfile, 'write', 'More details about the experiment can be found in the params.mat file');
        end
    end
    
    %% Start processing in the following order
     % Detect scratch from each frame
     % Generate velocity fields using PIV
     
     %tic;
     
     if(params.log)
         logger(params.logfile, 'write', ' ');
         logger(params.logfile, 'write', '*** Initializing PIV mask generation ***');
     end
     
     [params, dirs] = generatePIVMasks_par(params,dirs);
     
     if(params.log)
         logger(params.logfile, 'write', '*** Generated all masks for PIV ***');
         logger(params.logfile, 'write', ' ');
         logger(params.logfile, 'write', '*** Computing velocity fields ***');
     end
     
     [params, dirs] = velocityField_par(params,dirs);   
     
     if(params.log)
         logger(params.logfile, 'write', '*** Finished computing velocity fields ***');
         logger(params.logfile, 'write', ' ');
         logger(params.logfile, 'write', '*** Generating velocity field figures and kymographs ***');
     end
     
     %params.vMax = params.vMax*60;
     %save(paramMatFile, 'params', 'dirs', '-mat');
     [params, dirs] = makeFigs(params,dirs);
     
     if(params.log)
         logger(params.logfile, 'write', '*** Finished generating figures ***');
         logger(params.logfile, 'write', ' ');
         if(strcmpi(params.expType,'scatter')==0)
             logger(params.logfile, 'write', 'Looks like we are analyzing scatter assays today!');
             logger(params.logfile, 'write', 'There isn''t much else for me to do here then');
         else
            logger(params.logfile, 'write', '*** Starting spatial correlation analysis ***');
     
         end
     end
     
     if(strcmpi(params.expType,'scratch'))
        [params, dirs] = spatialCorrelation(params,dirs);
     
        pivmatfile = [dirs.pivDir filesep 'masterVels.mat'];
     
        if(params.log)
            logger(params.logfile, 'write', '*** Generated spatial correlation profiles and images ***');
            logger(params.logfile, 'write', ' ');
            logger(params.logfile, 'write', '*** Starting BOD analysis ***');
        end
        
     
        BOD(pivmatfile,[dirs.BOD filesep params.expName '_BOD'],1,144);
     
        if(params.log)
            logger(params.logfile, 'write', '*** Done with BOD analysis ***');
            logger(params.logfile, 'write', ' ');
            logger(params.logfile, 'write', '*** All done! Cleaning up... ***');
            logger(params.logfile, 'write', 'Good bye!');
            logger(params.logfile, 'close', ' ');
        end
     end
     save(paramMatFile, 'params','dirs','-mat');

     %toc;
     
     %clear all;
     %close all;
end


%{
[imgFile, imgPath] = uigetfile({'*.tif;*.tiff','Tif files (*.tif, *.tiff)'},'Select a file');
outDir = uigetdir('','Choose output directory');

min_interrog_win = 32;

if(~isequal(imgFile,0))
    
    stackFullFile = fullfile(imgPath, imgFile);
    stkInfo = imfinfo(stackFullFile);
    nFrames = numel(stkInfo);
    mImage=stkInfo(1).Width;
    nImage=stkInfo(1).Height;
    inpStack=zeros(nImage,mImage,nFrames,'uint8');
 
    TifLink = Tiff(stackFullFile, 'r');
    for i=1:nFrames
        TifLink.setDirectory(i);
        inpStack(:,:,i)=TifLink.read();
    end
    TifLink.close();
    
    desW = 512;
    desH = 512;

    nWindows = (desW*2)/min_interrog_win;
    kymo = zeros(nFrames-1,nWindows-1);
    pivm = zeros(nWindows-1,nWindows-1,nFrames-1);
    corr = zeros(nFrames-1,nWindows-1);
    s_pol = zeros(nFrames-1);
    for i=1:1:nFrames-1
       
        I1=inpStack(:,:,i);
        I2=inpStack(:,:,i+1);
        
        %desC = int16([(imgW-desW)/2 (imgH-desH)/2]);
        desC = int16([1 (nImage-desH)/2]);
        cropI1 = imcrop(I1,[desC(1) desC(2) desW-1 desH-1]);
        cropI2 = imcrop(I2,[desC(1) desC(2) desW-1 desH-1]);
        imwrite(cropI1,'tp1.tif');
        imwrite(cropI2,'tp2.tif');
        if i==1
            defineWorld(cropI1,0.625,'woco.mat');
        end
        olay=findWound(cropI1,cropI2,'woco.mat');%return;
        outFileName = sprintf('%s_%d','Edge',i);
        fullOutPath = fullfile(outDir,outFileName);
        saveas(olay,fullOutPath,'tif');
        close(olay);
 
        %[x,y,u,v,snr,pkh] = matpiv('tp1.tif','tp2.tif',[64 64;64 64;32 32;32 32;16 16;16 16;8 8;8 8],300,0.5,'multin','woco.mat','polymask.mat');
        [x,y,u,v,snr,pkh] = matpiv('tp1.tif','tp2.tif',[64 64;64 64;32 32;32 32],300,0.5,'multin','woco.mat','polymask.mat');
        [su,sv] = snrfilt(x,y,u,v,snr,1.3);
        [pu,pv] = peakfilt(x,y,su,sv,pkh,0.2);
        [gu,gv] = globfilt(x,y,pu,pv,2);
        [mu,mv] = localfilt(x,y,gu,gv,2,'median',3);
        load('polymask.mat');
        [fu,fv] = naninterp2(mu,mv,maske,x,y);
        %[x,y,u,v]=matpiv('tp1.tif','tp2.tif',64,600,0.5,'multi','woco.mat');
        %[x,y,u,v]=matpiv('MatPIV 1.7/Demo3/mpim1b.bmp','MatPIV 1.7/Demo3/mpim1c.bmp',64,0.0012,0.5,'multi','~/Desktop/PIV/MatPIV 1.7/Demo3/worldco.mat');
        %quiver(x,y,u,v);
        %Mp = x(:)';
        %Np = y(:)';
        %[M, N] = meshgrid(linspace(min(Mp),max(Mp),256), linspace(min(Np),max(Np),256));
        %P = griddata(x,y,u,M,N);
        %P(isnan(P))=0;
        %colormap('jet');
        %imagesc(P,colorbar)
        
        
        % The following generates a pseudocolored figure displaying PIV
        % 
        pivfig=figure('visible','off');
        w = magnitude(x,y,fu,fv);
        w(isnan(w))=0;
        colormap('jet');
        pcolor(flipud(x),flipud(y),w), axis square, shading flat, colorbar
        
        %
        
        %pivm(:,:,i)=w;
        %imshow(imread('tp1.tif'));
        %hold on
        %q=quiver(flipud(x),flipud(y),fu,fv);
        %q.LineWidth=2;
        %q.Color='blue';
        %hold off
        outFileName = sprintf('%s_%d','PIV',i);
        fullOutPath = fullfile(outDir,outFileName);
        saveas(pivfig,fullOutPath,'tif');
        close(pivfig);
        kymo(i,:)=mean(w);
        corr(i,:)=spCorr(fu);
        s_pol(i)=mean(mean(cos(atan2(fv./fu))))
        
    end
else
    return;
end
%}



