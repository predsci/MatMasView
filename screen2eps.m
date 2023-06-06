function screen2eps(h,filename,res)
%Sean P. McCarthy
if nargin < 2
  error('Not enough input arguments!')
end

dpi=96;

oldscreenunits = get(h,'Units');
oldpaperunits  = get(h,'PaperUnits');
oldpaperpos    = get(h,'PaperPosition');

set(h,'Units','pixels');
scrpos = get(h,'Position');
newpos = scrpos/dpi;
set(h,'PaperUnits','inches','PaperPosition',newpos)

print(h,'-depsc2', res, filename);

set(h,'Units',oldscreenunits,'PaperUnits',oldpaperunits,...
      'PaperPosition',oldpaperpos)
end
