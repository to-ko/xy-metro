% chi=susceptibility(theta)
%
% computes the susceptibility of a spin configuration given by theta.
% spin on site with lexic. index l is 
% [cos(theta(l)), sin(theta(l))]
%
% s(l)*s(k) = cos(theta(l))*cos(theta(k)) + sin(theta(l))*sin(theta(k))
%           = cos(theta(l)-theta(k))
function chi=susceptibility(theta)
   global L D
   chi = 0;
   for l=1:L^D
       for k=1:L^D
          chi = chi + cos(theta(l)-theta(k));
       end
   end
   chi = chi / L^D;
end