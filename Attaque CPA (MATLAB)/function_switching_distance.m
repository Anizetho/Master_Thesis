function [sd] = function_switching_distance(V1, V2, delta, mode)
    for byte=1:16
        dist(:,[((256*mode)-(256*((mode+1)-byte)))+1:((256*mode)-(256*(mode-byte)))]) =  bitxor(V1(:,i),V2(:,[((256*mode)-(256*((mode+1)-byte)))+1:((256*mode)-(256*(mode-byte)))]));
        C = [((256*mode)-(256*((mode+1)-byte)))+1:((256*mode)-(256*(mode-byte)))];
        
        hd(:,C) = reshape(sum(de2bi(dist(:,C)),2), size(dist(:,C)));
        tmp1(:,C) = bsxfun(@bitand, V1(:,i), hd(:,C));
    
        hw_tmp1(:,C) = sum(de2bi(tmp1(:,C)), 2);
        hw_tmp1(:,C) = reshape(hw_tmp1(:,C),size(tmp1(:,C)));
        sd(:,C) = hd(:,C) - delta*hw_tmp1(:,C);
    end
    
    %final = hd(:,C) - delta*hw_tmp1(:,C);
            
            
%     dist =  bitxor(V1,V2);
%     hd = reshape(sum(de2bi(dist),2), size(dist));
%     
%     tmp1 = bsxfun(@bitand, V1, hd);
%     
%     hw_tmp1 = sum(de2bi(tmp1), 2);
%     hw_tmp1 = reshape(hw_tmp1,size(tmp1));
%     %hw_hd = sum(de2bi(hd), 2);
%     sd = hd - delta*hw_tmp1;
    
end

