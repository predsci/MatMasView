function [xmin,xmax]=get_nice_axis_lim(x,p)

  xmax = max(x(:));
  xmin = min(x(:));
  xrange = xmax-xmin;
  xbuf = p*xrange;
  
  xmax=xmax+xbuf;
  xmin=xmin-xbuf;

return
end