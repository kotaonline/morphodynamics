function defineWorld(Img,pscl,matfile)
    if(ischar(Img))
        wI = imread(Img);
    else
        wI = Img;
    end
    [wIx, wIy] = size(wI);
    
    if (mod(wIx,2)==0)
        rc = wIx/2;
    else
        rc = (wIx+1)/2;
    end
    
    if (mod(wIy,2)==0)
        hc = wIy/2;
    else
        hc = (wIy+1)/2;
    end
    
    wIc = [20 20;...
           20 hc;...
           20 wIy-20;...
           rc 20;...
           rc hc;...
           rc wIy-20;...
           wIx-20 20;...
           wIx-20 hc;...
           wIx-20 wIy-20];
       
    iIc = [ones(size(wIc,1),1) wIc*pscl];
    
    comap = iIc\wIc(:,:);
    
    if (ischar(matfile))
        save(matfile,'comap');
    else
        printf('Error! <matfile> must be a valid file name');
    end
end