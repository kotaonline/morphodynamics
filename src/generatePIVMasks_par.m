function [params, dirs] = generatePIVMasks_par(params, dirs)

    %% Process one frame at a time
    local_ifrdir = dirs.iframes;
    local_ext = params.ext;
    
    parfor i = 1 : params.nFrames-1
        ifname1 = [ local_ifrdir filesep sprintf('%03d',i) local_ext ];
        ifname2 = [ local_ifrdir filesep sprintf('%03d',i+1) local_ext ];

        Img1 = imread(ifname1);
        Img2 = imread(ifname2);
        
        Img = imsubtract(Img1,Img2);
        
        [params, dirs] = generateMask(Img,Img1,i,params,dirs);
    end

end