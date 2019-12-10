function imageReg_wrapper(fileList,outDir,x,y,w,h,s,e)
%% Input Parameters
% Input is provided as a newline separated list of filenames that need to
% be registered. Absolute paths to files are preferable while relative
% paths are OK. If using bare filenames, user is responsible for setting
% paths to code.

% Registration appears to work better if an invariant background field
% devoid of cells is provided. In order to provide this input, please input
% x,y co-ordinates of the top left of the ROI, width (w) and height (h) of
% the rectangular ROI. At this time, only rectangular ROIs are processed.
% Optionally start (s) and end (e) frames can also be provided as input to
% register a subset of images.
    
    fileH = fopen(fileList);
    nextLine = fgetl(fileH);
    while(ischar(nextLine))
        fprintf('Current File: %s ',nextLine);
        
        registerInput(nextLine,outDir,x,y,w,h,s,e);
        fprintf(' - Done\n');
        nextLine = fgetl(fileH);
    end
    fclose(fileH);
end

function registerInput(imgFile,oDir,x,y,w,h,s,e)
    
    [fPar, fName, fExt] = fileparts(imgFile);
    stkInfo = imfinfo(imgFile);
    nFrames = numel(stkInfo);
    imgW    = stkInfo(1).Width;
    imgH    = stkInfo(1).Height;
    inpStack = uint8(zeros(imgH,imgW,nFrames));
    regStack = uint8(zeros(imgH,imgW,nFrames));
    outTif   = [oDir filesep fName '-reg' fExt];
    
    TifLink = Tiff(imgFile, 'r');

    for i = 1 : nFrames
        TifLink.setDirectory(i);
        I=TifLink.read();
        inpStack(:,:,i) = I;
    end
    
    TifLink.close();
    if(x==-1)
        x=1;
    end
    if(y==-1)
        y=1;
    end
    if(w==-1)
        w=imgW;
    end
    if(h==-1)
        h=imgH;
    end
    if(s==-1)
        s=1;
    end
    if(e==-1)
        e=nFrames;
    end
   
    %regStack(:,:,1)  = inpStack(:,:,1);
    regStack = inpStack(:,:,s:e);
    fftIm1 = fft2(regStack(:,:,1));
    fftreg1 = fft2(regStack(y:y+h-1,x:x+w-1,1));
    imwrite(regStack(:,:,1),outTif);
    [nr,nc]=size(fftIm1);
    Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
    Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
    [Nc,Nr] = meshgrid(Nc,Nr);
    for i=2:e
        fprintf(['Processing frame - ' sprintf('%03d\n',i)]);

        fftIm2  = fft2(regStack(:,:,i));
        fftreg2 = fft2(regStack(y:y+h-1,x:x+w-1,i));
        [output, fout] = dftregistration(fftreg1,fftreg2,100);
        

        Greg = fftIm2.*exp(1i*2*pi*(-output(3)*Nr/nr-output(4)*Nc/nc));
        Greg = Greg*exp(1i*output(2));

        regStack(:,:,i) = uint8(abs(ifft2(Greg)));
        imwrite(regStack(:,:,i),outTif,'WriteMode','append');
    
    end
    
end


