% function x=alexic(l)
% computes the coordinates of a D dimensional vector, given the 
% lexicographical index l, on a lattice of size L x L x L x ... x L
% The dimension D and lattice size L are global variables.
%
% alexic computes the inverse of
% l = x(1) + (x(2)-1)*L + (x(3)-1)*L*L + ...
function x=alexic(l)
   global L D
   l=l-1;
   x=zeros(D,1);
   for d=1:D
       x(d) = mod(l,L)+1;
       l= floor(l/L);
   end
end