function ret=my_prctile(data, v)
% ret=my_prctile(data, v)

  if all(size(data)>1)%apply percentile on columns
    assert(length(size(data))==2)
    ret = [];
    for k=1:size(data,2)
      x = sort(data(:,k));
      vec = [];
      for j=1:length(v)
        vec = [vec percentile(x, v(j))];
      end
      ret = [ret vec'];
    end
  elseif any(size(data)==0)
	ret = [];
  else
    x = sort(data);
    vec = [];
    for j=1:length(v)
      vec = [vec percentile(x, v(j))];
    end
    ret = vec;
  end
return


function perc = percentile(x, k_percentile)

  %x = sort(data);
  num_x = length(x);
  k = k_percentile/100 * (num_x-1) + 1;
  
  [k_int, D] = rat(k);
  if isequal(D,1),  % k is an integer, take percentile element
    perc = x(k_int);
  else              % take midpoint between two elements
    perc = mean(x(floor(k):ceil(k)));
  end
return
