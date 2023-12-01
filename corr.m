% G=corr(theta)
%
% computes the correlation function
% G(t) = S(0)*S(t)
% where
% S is a slice-spin
% S(t) = 1/L sum_x s(x,t)
% spin on site with lexic. index l is 
% [cos(theta(l)), sin(theta(l))]
%
function G=corr(theta)
   global L D
   G = zeros(L,1);
   
   % compute time-slice spins
   for x1=1:L
      s{x1} = [0,0];
      for x2=1:L
         l = lexic([x1,x2]);
         s{x1} = s{x1} + [cos(theta(l)), sin(theta(l))];
      end
      x{x1} = s{x1} / L;
   end
   
       
   for t=0:L-1
      for x1=1:L
         G(t+1) = G(t+1) + s{x1} * s{mod(x1+t-1,L)+1}.';
      end
      G(t+1) = G(t+1)/L;
   end
end
