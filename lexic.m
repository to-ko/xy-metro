% function l=lexic(x)
%
% computes the lexicographical index for a D dimensional vector x
% The lattice size is L x L x L x ... x L, where L is a global variable
% D is arbitrary
%
% l = x(1) + (x(2)-1)*L + (x(3)-1)*L*L + ...
function l=lexic(x)
   global L
   l=0;
   for d=length(x):-1:1
       l=l*L;
       l=l+x(d)-1;
   end
   l=l+1;
end