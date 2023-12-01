% H=energy(theta)
%
% computes the energy of a spin configuration given by theta.
% spin on site with lexic. index l is 
% [cos(theta(l)), sin(theta(l))]
%
% the dot product between sites l and k is
%
% s(l)*s(k) = cos(theta(l))*cos(theta(k)) + sin(theta(l))*sin(theta(k))
%           = cos(theta(l)-theta(k))
function H=energy(theta)
   global h L D
   H = 0;
   for l=1:L^D
      for mu=1:D
         k = h(l,mu);
         H = H - cos(theta(l)-theta(k));
      end
   end
end