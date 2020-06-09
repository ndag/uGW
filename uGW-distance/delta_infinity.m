function [result] = delta_infinity(x,y)
% This function calculates dtelat_infinity between (x and y)
  
   if abs(x-y)<10^(-15)
       result=0;
   else
       result=max(x,y);
   end
   
end


