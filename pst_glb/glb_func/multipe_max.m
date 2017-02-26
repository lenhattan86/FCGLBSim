function [ maxes idxes ] = multipe_max( array, top_num)
len = length(array);
    
   max_prev = inf;
   for i=1:top_num
      max_val = -inf;  
      for j=1:len
         if max_val < array(j) && array(j)<max_prev
            max_val = array(j);
            maxes(i) = max_val;
            idxes(i) = j;
         end
      end
      max_prev = max_val;
   end
end

