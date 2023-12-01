function [acc, Iam1] = sweep()
   global L D beta h theta
   Delta = pi/4; % interval for proposals
   
   acc = zeros(L^D,1);
   DH  = zeros(L^D,1);
   for l=1:L^D
      % compute local magnetization
      b = zeros(D,1);
      for mu=1:D
         % contribution from positive neighbors
         k = h(l,mu);
         b = b + [cos(theta(k)); sin(theta(k))];
         % contribution from negative neighbors
         k = h(l,D+mu);
         b = b + [cos(theta(k)); sin(theta(k))];
      end
      
      % propose change
      phi = (rand-0.5)*2*Delta; % uniformly distributed in [-Delta,+Delta]
      tlp = theta(l)+phi;
      s   = [cos(theta(l)); sin(theta(l))]; % old value of spin on site l
      sp  = [cos(tlp); sin(tlp)];          % proposed new value
      DH(l) = (s'-sp')*b;
      
      % accept-reject step
      A = exp(-beta*DH(l));
      if rand < A
         % proposal accepted
         acc(l) = 1;
         theta(l) = tlp; % update the field
      else
         % proposal rejected
         acc(l) = 0; 
      end
   end
   Iam1 = mean(exp(-beta*DH));
   acc  = sum(acc); % number of accepted proposals
end

