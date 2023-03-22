% qselect - select the k-th smallest out of n numbers
% Implements Hoare's Quickselect algorithm (https://en.wikipedia.org/wiki/Quickselect)
% with the median-of-3 pivot strategy. Operates in-place, avoiding sorting and recursion
%
%    kth = qselect(a, k)
%
%    a      - array of n elements
%    k      - specifies the desired k-th smallest element
%             the k-th *largest* element can be found by passing length(a)+1-k
%
% Returns
%    kth    - the sought element

% sample use:
%
% mat=[3, 2, 7, 4, 5, 1, 4, -1];
% fprintf('second smallest element is %g\n', qselect(mat, 2));                % 1
% fprintf('fourth largest element is %g\n', qselect(mat, length(mat)+1-4));   % 4

% Manolis Lourakis 2007-18
% Institute of Computer Science, Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece

% Sep  2018  - Original version. (v. 1.0)

% Adapted by Alexandra Gallyas-Sanhueza (ag753@cornell.edu) to return the
% number of operations (OPS), and to return the two values around the
% median if given k1 = floor((D+1)/2), k2 = ceil((D+1)/2)

function [k2th,k1th,OPS] = qselect_modified(a, k1, k2)
    
  k = k2;

  if(~isvector(a))
    error('qselect: expecting a vector input');
  end

  l=1; r=length(a);
  %if(k<0), kth=qselect(a, r+1+k); return; end  % return k-th largest element for negative k

  OPS = 0;
  
  while(true)
    % median of three, i.e. a(l), a(s), a(r), as pivot
    % this optimization generally improves performance and
    % eliminates worst-case behavior for sorted/reverse-sorted data
    s=fix((r + l)/2); % 1 
    if(a(s)<a(r)), temp=a(s); a(s)=a(r); a(r)=temp; end % 1
    if(a(s)<a(l)), temp=a(s); a(s)=a(l); a(l)=temp; end % 1
    if(a(r)<a(l)), temp=a(r); a(r)=a(l); a(l)=temp; end % 1
    pivot=a(r); % median
    
    OPS = OPS + 4;

    i=l;
    for j=l:r-1 
      if(a(j)<=pivot) 
        temp=a(i); a(i)=a(j); a(j)=temp; 
        i=i+1; 
        OPS = OPS + 1;
      end
      OPS = OPS + 1;
    end

    temp=a(r); a(r)=a(i); a(i)=temp; 

    s=i-l+1;
    OPS = OPS + 1;
    if(k<s)
      r=i-1;
      OPS = OPS + 1;
    elseif(k>s)
      l=i+1; k=k-s;
      OPS = OPS + 2;
    else
      k2th=a(i);
      k1th=max(a(1:k1));
      OPS = OPS + k1;
      return;
    end
    OPS = OPS + 2;
  end
end

