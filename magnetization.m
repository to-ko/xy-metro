% M=magnetization(theta)
%
% computes the magnetization of a spin configuration given by theta.
% spin on site with lexic. index l is 
% [cos(theta(l)), sin(theta(l))]
%
function M=magnetization(theta)
   global L D
   M = [0, 0];
   for l=1:L^D
      M = M + [cos(theta(l)), sin(theta(l))];
   end
   M = norm(M) / L^D;
end