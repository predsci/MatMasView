function draggingFncP(hfig,events)
    
   global pidx pval;

   N = evalin('base','length(pvec_cut )');  
           
   pidx = pidx - events.VerticalScrollCount;    
      
   if(pidx > N)
      pidx = pidx-N;
   end
   if(pidx < 1)
      pidx = pidx+N;
   end
      
   pval = evalin('base','pvec_cut (pidx);');
   
   %Set 3D cut frame:
   evalin('base','set(h_pcut,''CData'',squeeze(MAS_DATA_CUT(pidx,:,:)))');
   evalin('base','set(h_pcut,''XData'',squeeze(X(pidx,:,:)))');
   evalin('base','set(h_pcut,''YData'',squeeze(Y(pidx,:,:)))');
   evalin('base','set(h_pcut,''ZData'',squeeze(Z(pidx,:,:)))');
      
   %Set 2D cut plots CDATA (plot selection done in startdragfunc)
   evalin('base','[X_cut Y_cut]=pol2cart(squeeze(T(pidx,:,:)),squeeze(R(pidx,:,:)));');
   evalin('base','set(plot_slice,''CData'',fliplr(squeeze(MAS_DATA_CUT(pidx,:,:))''),''XData'',fliplr(Y_cut''),''YData'',fliplr(X_cut''))');
   evalin('base','title(axis_slice,[''\phi='',num2str(pval)],''Color'',''w'');');
   evalin('base','set(fig_2Dslice, ''Name'',[field_name,'' R-\theta \phi='',num2str(pval)]);');

   %Set YDATA in 1D plots (plot selection done elsewhere)
   evalin('base','set(plot_1D,''YData'',squeeze(MAS_DATA_CUT(pidx,tidx,:)))');
%   evalin('base','set(axis_1Dplot,''YLim'',[min(squeeze(MAS_DATA_CUT(pidx,tidx,:))) max(squeeze(MAS_DATA_CUT(pidx,tidx,:)))] );');
   evalin('base','pval = pvec_cut (pidx);');
   evalin('base','title(axis_1Dplot,[''\phi='',num2str(pval),''      \theta='',num2str(tval)],''Color'',''k'');');
   evalin('base','set(fig_1Dplot,''Name'',[field_name,'' T='',num2str(tval),'' P='',num2str(pval)]);');

   
         
end
