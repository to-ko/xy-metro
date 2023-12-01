function plot_cnfg( theta )
   global D L h
   
   if D~=2
      error('plotting only possible with D=2'); 
   end
   
   figure('Color',[1 1 1]);
   axis([0 L+1 0 L+1]); hold on
   caxis([-1 1])
   colorbar
   col = colormap;
   % plot links: color according to energy contribution
   for l=1:L^D
       x = alexic(l);
       for mu=1:D
          y = alexic(h(l,mu));
          e = cos(theta(l) - theta(h(l,mu)));
          % compute index in colormap 'jet'
          Ci=fix((e*0.9999+1)/2*size(col,1))+1;
          if x(1)==L && y(1)==1, y(1) = L+0.5; end % treat boundary points 
          if x(2)==L && y(2)==1, y(2) = L+0.5; end
          line([x(1) y(1)],[x(2) y(2)],'Color',col(Ci,:),'LineWidth',2);
       end
   end
   
   % plot sites
   for l=1:L
       for k=1:L
           plot(l,k,'k.')
       end
   end
   
       
   sx = reshape(cos(theta),L,L)';
   sy = reshape(sin(theta),L,L)'; 
   quiver([1:L],[1:L],sx,sy,0.4,'Color',[0 0 0]);
   xlabel('x_1')
   ylabel('x_2')
   title(sprintf('spin cnfg on a %d x %d torus',L,L));
end

