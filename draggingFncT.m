function draggingFncT(hfig,events)
    
   global tidx tval;

   N = evalin('base','length(tvec_cut)');      

   tidx = tidx - events.VerticalScrollCount;    
      
   if(tidx > N)
      tidx = tidx-N;
   end
   if(tidx < 1)
      tidx = tidx+N;
   end
      
   tval = evalin('base','tvec_cut (tidx)');
      
   evalin('base','set(h_tcut,''CData'',squeeze(MAS_DATA_CUT(:,tidx,:)))');
   evalin('base','set(h_tcut,''XData'',squeeze(X(:,tidx,:)))');
   evalin('base','set(h_tcut,''YData'',squeeze(Y(:,tidx,:)))');
   evalin('base','set(h_tcut,''ZData'',squeeze(Z(:,tidx,:)))');
      
   %Set 2D cut plots CDATA (plot selection done in startdragfunc)
   evalin('base','[X_cut Y_cut]=pol2cart(squeeze(P(:,tidx,:)),squeeze(R(:,tidx,:)));');
   evalin('base','set(plot_slice,''CData'',squeeze(MAS_DATA_CUT(:,tidx,:)),''XData'',X_cut,''YData'',Y_cut);');
   evalin('base','title(axis_slice,[''\theta='',num2str(tval)],''Color'',''w'');');
   evalin('base','set(fig_2Dslice, ''Name'',[field_name, ''R-\phi \theta='',num2str(tval)]);');

   %Set YDATA in 1D plots (plot selection done elsewhere)
   evalin('base','set(plot_1D,''YData'',squeeze(MAS_DATA_CUT(pidx,tidx,:)))');
   evalin('base','set(axis_1Dplot,''YLim'',[min(squeeze(MAS_DATA_CUT(pidx,tidx,:))) max(squeeze(MAS_DATA_CUT(pidx,tidx,:)))] );');
   evalin('base','tval = tvec_cut (tidx);');
   evalin('base','title(axis_1Dplot,[''\phi='',num2str(pval),''      \theta='',num2str(tval)],''Color'',''k'');');
   evalin('base','set(fig_1Dplot,''Name'',[field_name,'' T='',num2str(tval),'' P='',num2str(pval),''  Frame: '',num2str(id,''%03d'')]);');
   
         
end
