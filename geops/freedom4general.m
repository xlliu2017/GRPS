    function [fidx,bidx] = freedom4general(DinT,idx,Dtype)
        
        switch Dtype
            case 'P'
                flag = sum(DinT(:,idx),2);
                fidx = find(flag == 6);
                bidx = find((flag>=1).*(flag <=3));
            case 'E'
                flag = sum(DinT(:,idx),2);
                fidx = find(flag == 2);
                bidx = find(flag == 1);
            case 'T'
                fidx = idx;
                bidx = [];
        end        
    end