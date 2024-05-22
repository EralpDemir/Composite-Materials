function ANG_perm = LP_permute_blnc(ANG, nPly)

k = 0;
for i = 1:size(ANG,1)
    if nPly == 2
        if ((ANG(i,1)>0 && ANG(i,2)<0) || (ANG(i,1)<0 && ANG(i,2)>0))&& ...      
                (ANG(i,1)==-ANG(i,2))
            k = k+1;
            ANG_perm(k,:) = ANG(i,:);
        elseif (ANG(i,1)==0 && ANG(i,2)==-90)
            k = k+1;
            ANG_perm(k,:) = [-90 0];
        elseif (ANG(i,1)==-90 && ANG(i,2)==0)
            k = k+1;
            ANG_perm(k,:) = [0 -90];
        end
        
    elseif nPly == 4
        if ((ANG(i,1)>0 && ANG(i,2)>0 && ANG(i,3)<0 && ANG(i,4)<0) || ...
                (ANG(i,1)>0 && ANG(i,2)<0 && ANG(i,3)>0 && ANG(i,4)<0) || ...
                (ANG(i,1)>0 && ANG(i,2)<0 && ANG(i,3)<0 && ANG(i,4)>0) || ...
                (ANG(i,1)<0 && ANG(i,2)<0 && ANG(i,3)>0 && ANG(i,4)>0) || ...
                (ANG(i,1)<0 && ANG(i,2)>0 && ANG(i,3)<0 && ANG(i,4)>0) || ...
                (ANG(i,1)<0 && ANG(i,2)>0 && ANG(i,3)>0 && ANG(i,4)<0)) && ...      
                ((ANG(i,1)==-ANG(i,3) && ANG(i,2)==-ANG(i,4)) || ...
                (ANG(i,1)==-ANG(i,4) && ANG(i,2)==-ANG(i,3)))
            k = k+1;
            ANG_perm(k,:) = ANG(i,:);
        elseif (ANG(i,1)==0 && ANG(i,2)==0 && ANG(i,3)==-90 && ANG(i,4)==-90)
            k = k+1;
            ANG_perm(k,:) = [-90 -90 0 0];
        elseif (ANG(i,1)==0 && ANG(i,2)==-90 && ANG(i,3)==0 && ANG(i,4)==-90)
            k = k+1;
            ANG_perm(k,:) = [-90 0 -90 0];
        elseif (ANG(i,1)==0 && ANG(i,2)==-90 && ANG(i,3)==-90 && ANG(i,4)==0)
            k = k+1;
            ANG_perm(k,:) = [-90 0 0 -90];
        elseif (ANG(i,1)==0 && ANG(i,2)==-90 && ANG(i,3)==-90 && ANG(i,4)==0)
            k = k+1;
            ANG_perm(k,:) = [-90 0 0 -90];
        elseif (ANG(i,1)==-90 && ANG(i,2)==-90 && ANG(i,3)==0 && ANG(i,4)==0)
            k = k+1;
            ANG_perm(k,:) = [0 0 -90 -90];
        elseif (ANG(i,1)==-90 && ANG(i,2)==0 && ANG(i,3)==-90 && ANG(i,4)==0)
            k = k+1;
            ANG_perm(k,:) = [0 -90 0 -90];
        elseif (ANG(i,1)==-90 && ANG(i,2)==0 && ANG(i,3)==0 && ANG(i,4)==-90)
            k = k+1;
            ANG_perm(k,:) = [0 -90 -90 0];

        end
    end
end
