
if psixSL > psixPL % SFD-Plus
    
    if size(spRZLP,1) == 3
        spRZLP = spRZLP([1 3],:);
    elseif size(spRZLP,1) == 4
        spRZLP = spRZLP([1 4],:);
    end
    if size(spRZLS,1) == 4
        spRZLS = spRZLS([2 3],:);
    end
    
else % SFD-Minus
    
    if rxPL < rxSL % LFS
        
        if size(spRZLP,1) == 3
            spRZLP = spRZLP([1 2],:);
        end
        if size(spRZLP,1) == 4
            spRZLP = spRZLP([1 2],:);
        end
        if size(spRZLS,1) == 4
            spRZLS = spRZLS([3 4],:);
        end
           
    else % HFS        
        
        if size(spRZLP,1) == 3
            spRZLP = spRZLP([2 3],:);
        end
        if size(spRZLP,1) == 4
            spRZLP = spRZLP([3 4],:);
        end
        if size(spRZLP,1) == 5
            spRZLP = spRZLP([4 5],:);
        end
        if (size(spRZLS,1) == 3) || (size(spRZLS,1) == 4)
            spRZLS = spRZLS([1 2],:);
        end
       
    end
   
end
