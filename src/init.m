function [params, dirs] = init(inpFile,paramFile)
    %% Check user input file. Currently I check whether the extension is .tif or .tiff and assume that it is a tif file if the extension matches.

    [ parentDir, expName, ext ] = fileparts(inpFile);

    if(~strcmp(ext,'.tif') && strcmp(ext,'.tiff'))
        error(strcat('File format - ',ext,' - not supported. Please choose a .tif or .tiff file'));
    end
    
    %% Read parameter file and populate params struct
    
    paramFileHdl=fopen(paramFile);
    nextLine = fgetl(paramFileHdl);
    while ischar(nextLine)
        lineParts = strsplit(nextLine,'=');
        params.(lineParts{1}) = lineParts{2};
        nextLine = fgetl(paramFileHdl);
    end
    fclose(paramFileHdl);
    
    % expName is the same as the tif file name without extension
    params.expName = expName;
    
    % Camera pixel size and frame rate are essential parameters. Units for
    % pixel size are um and frame rate is expressed in minutes
    
    if(~isfield(params,'pixelSize') || ~isfield(params,'timePerFrame'))
        error('pixelSize and timePerFrame are required parameters.');
    end
    
    % params.overlap determines extent of overlap between interrogration
    % windows in PIV analysis. If not set in the input text file, this
    % value is set to 0.5 (50% overlap)
    
    if(~isfield(params,'overlap'))
        params.overlap = 0.5;
    else
        params.overlap = str2double(params.overlap);
    end
    
    % params.overwrite determines whether computation must overwrite
    % existing files at each step. Default behavior is to skip
    % (re)computation if output files exist. One can get around this by
    % 'touch'ing necessary output files. However, I strongly recommend
    % against this since some analysis depends not only on the 'existence'
    % of other files but their contents as well.
    
    if(~isfield(params,'overwrite'))
        params.overwrite = 0;
    else
        params.overwrite = str2num(params.overwrite);
    end
    
    % This is a dummy parameter. Intended for future expansion of code to
    % include migration in both directions. Currently migration is assumed
    % to proceed only in the positive x-direction.
    
    if(~isfield(params,'edges'))
        params.edges = 2;
    else
        params.edges = str2num(params.edges);
    end
    
    % This parameter is not fully supported either. The code was developed
    % with scratch and scatter assays in mind. Scratch assays are tested
    % more rigorously than scatter assays.
    
    if(~isfield(params,'expType'))
        params.expType = 'scratch';
    end
    
    % arbirtrary cutoff - throwing away small areas not filled with cells.
    params.bound_cutoff = 500; 
    
    %arbitrary - hope that clusters are not smaller than this size.
    if(strcmp(params.expType,'scatter'))
        params.bound_cutoff = 200; 
    end
    
    % Smoothing kernel size - used for smoothing 2d matrices using smooth2a
    % while generating figures
    
    if(~isfield(params,'smoothX'))
        params.smoothX      = 5;
    end
    
    if(~isfield(params,'smoothY'))
        params.smoothY      = 5;
    end
    
    params.pixelSize    = str2double(params.pixelSize);
    params.timePerFrame = str2double(params.timePerFrame);
    params.pxperum      = 1/params.pixelSize;
    params.sMax         = 0.0;
    
    params.inpFile      = inpFile;
    
    % Here I read the input tif info to set basic parameters such as width
    % and height of the frame and number of frames in the movie.
    
    stkInfo = imfinfo(inpFile);
    params.nFrames = numel(stkInfo);
    params.width = stkInfo(1).Width;
    params.height = stkInfo(1).Height;
    
    params.ext = ext;
    
    % woco stands for 'world coordinates'. It is the mapping between camera
    % pixel coordinates (pixels) and real world values (um). The mat file
    % contains a vector with x-y scaling info
    
    params.woco = [parentDir filesep expName filesep 'woco.mat'];
    
    % This parameter is not fully supported either. It is meant for
    % starting analysis at an arbitrary value other than the first frame.
    % It defaults to the first frame, so no harm done if not set.
    
    if(~isfield(params,'sframe'))
        params.sframe=1;
    else
        params.sframe = str2num(params.sframe);
    end
    
    % Similar to sframe
    
    if(~isfield(params,'eframe'))
        params.eframe = params.nFrames-1;
    else
        params.eframe = str2num(params.eframe);
    end

    % This parameter is set for determining the cutoffs for spatial
    % correlation analysis. 
    
    if(~isfield(params,'radAnnulus'))
        params.radAnnulus = 40;
    end
    
    %% Generate output directories with input filename as the parent directory
    % Output directory structure:
    %
    % parentDir - directory in which input tif file is present
    %  |
    %  ---> expDir - has the same name as input tif file without ext
    %  |
    %  ---> Preprocess
    %        |
    %        ---> iframes - individual frames in the multiframe tiff
    %        |
    %        ---> overlay - tifs overlaid with detected scratch edge
    %        |
    %        ---> maskmat - mat file masks removing areas with no cells
    %        |
    %        ---> masktif - same mask files as in mat, but in tif format
    %  |
    %  ---> PIV - contains all PIV analysis
    %        |
    %        ---> velFields - contains images with quiver plots
    %        |
    %        ---> OrderParam - contains images with quiver plots
    %  |
    %  ---> Postprocess
    %        |
    %        ---> Kymograph - contains velocity and orderparam kymographs
    %        |
    %        ---> DistCorr  - contains spatial correlation data and plots
    %        |
    %        ---> BOD - contains biorthogonal decomposition data and plots
    %  |
    %  ---> log - temporary directory to hold all log files in parallel
    %  loops. Note that this directory is deleted once log files are
    %  consolidated.
    
    dirs.parentDir = [parentDir filesep];
    dirs.expDir = [dirs.parentDir filesep expName filesep];
    
    % Directory for preprocessing images
    
    dirs.preproc = [dirs.expDir filesep 'Preprocess' filesep];
    
    dirs.iframes = [dirs.preproc filesep 'iframes' filesep];
    dirs.overlay = [dirs.preproc filesep 'overlay' filesep];
    dirs.masktif = [dirs.preproc filesep 'masktif' filesep];
    dirs.maskmat = [dirs.preproc filesep 'maskmat' filesep];
    
    % Directory for PIV analysis
    
    dirs.pivDir   = [dirs.expDir filesep 'PIV' filesep];
    dirs.vfImages = [dirs.pivDir filesep 'velFields' filesep];
    dirs.SImages  = [dirs.pivDir filesep 'OrderParam' filesep];
    
    % Directory for Postprocessing
    
    dirs.postproc = [dirs.expDir filesep 'PostProcess' filesep];
    
    dirs.kymo     = [dirs.postproc filesep 'Kymographs' filesep];
    dirs.spCorr   = [dirs.postproc filesep 'DistCorr' filesep];
    dirs.BOD      = [dirs.postproc filesep 'BOD' filesep];
    
    if(~exist(dirs.expDir,'dir'))
        mkdir(dirs.expDir);
    end
    
    % Create preprocess directories
    
    if(~exist(dirs.preproc, 'dir'))
        mkdir(dirs.preproc);
    end
    
    if(~exist(dirs.iframes, 'dir'))
        mkdir(dirs.iframes);
    end
    
    if(~exist(dirs.overlay,'dir'))
        mkdir(dirs.overlay);
    end
    
    if(~exist(dirs.maskmat,'dir'))
        mkdir(dirs.maskmat);
    end
    
    if(~exist(dirs.masktif,'dir'))
        mkdir(dirs.masktif);
    end
    
    % Directory for PIV analysis
    
    if(~exist(dirs.pivDir,'dir'))
        mkdir(dirs.pivDir);
    end
    
    if(~exist(dirs.vfImages, 'dir'))
        mkdir (dirs.vfImages);
    end
    
    if(~exist(dirs.SImages, 'dir'))
        mkdir (dirs.SImages);
    end
    
    % Directories for postprocessing
    
    if(~exist(dirs.postproc, 'dir'))
        mkdir(dirs.postproc);
    end
    
    if(~exist(dirs.kymo, 'dir'))
        mkdir(dirs.kymo);
    end
    
    if(~exist(dirs.spCorr, 'dir'))
        mkdir(dirs.spCorr);
    end
    
    if(~exist(dirs.BOD, 'dir'))
        mkdir (dirs.BOD);
    end
    
    % Directory for logging is created only if logging is enabled. By
    % default, I write output to log files
    
    if(~isfield(params,'log'))
        params.log = 1;
    end
    
    if(params.log)
        %dirs.log = [dirs.expDir filesep 'log'];
        params.logfile = [dirs.expDir filesep params.expName '.log'];
        %if(~exist(dirs.log,'dir'))
        %    mkdir(dirs.log);
        %end
    end

    %% Convert input stack to images

    fname = [parentDir filesep expName ext];

    if ~exist(fname,'file')    
        t = 1;
        while exist([dirs.iframes sprintf('%03d',t) '.tif'],'file')
            t = t + 1;
        end
        nFrames = t - 1;
    
        if nFrames < 2    
            error('File %s nor images exist',fname);
        end
    end
    
    TifLink = Tiff(inpFile, 'r');
    for i = 1 : params.nFrames
        
        TifLink.setDirectory(i);
        
        if(~exist(params.woco,'file'))
            %fprintf('Define world coordinates - ');
            % Define world coordinates -- Needed for MatPIV to scale images properly
            I=TifLink.read();
            defineWorld(I, params.pxperum, params.woco);
            %fprintf(' Done\n');
        end
        
        if(~exist([dirs.iframes filesep sprintf('%03d',i) params.ext],'file'))
            I=TifLink.read();
    
            eval(['imwrite(I,''' [dirs.iframes sprintf('%03d',i) '.tif'''] ',''tif'')']);
        end
    end
    TifLink.close();

end