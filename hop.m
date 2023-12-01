% function h=hop()
% generate the neighborhood relations on a D dimensional Torus of size L
% in each direction
%
% returns a L^D x 2*D dimensional matrix, such that
% h(l,mu) is the lexicographical index of a nearest neighbor of the point with
% index l. mu=1..D give the positive neighbors, mu=D+1..2*D gives the 
% negative neighbors
function h=hop()
global L D
   h = zeros(L^D,2*D);
   for l=1:L^D
      x = alexic(l);
      for mu=1:D
         xp = x;
         % positive neighbor in mu direction
         xp(mu) = mod(x(mu),L)+1;
         h(l,mu) = lexic(xp);
         % negative neighbor in mu direction
         xp(mu) = mod(x(mu)-2+L,L)+1;
         h(l,D+mu) = lexic(xp);
      end
   end
end