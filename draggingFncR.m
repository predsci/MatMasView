function draggingFncR(hfig,events)
    
   global ridx rval;

   N = evalin('base','length(rvec_cut)');  
           
   ridx = ridx - events.VerticalScrollCount;    
      
   if(ridx > N)
      ridx = ridx-N;
   end
   if(ridx < 1)
      ridx = ridx+N;
   end
      
   rval = evalin('base','rvec_cut (ridx)');
     
   evalin('base','set(h_rcut,''CData'',squeeze(MAS_DATA_CUT(:,:,ridx)))');
   evalin('base','set(h_rcut,''XData'',squeeze(X(:,:,ridx)))');
   evalin('base','set(h_rcut,''YData'',squeeze(Y(:,:,ridx)))');
   evalin('base','set(h_rcut,''ZData'',squeeze(Z(:,:,ridx)))');      
      
   %Set 2D cut plots CDATA (plot selection done in startdragfunc)
   evalin('base','set(plot_slice,''CData'',squeeze(MAS_DATA_CUT(:,:,ridx)),''XData'',P(:,:,ridx),''YData'',T(:,:,ridx))');
   evalin('base','title(axis_slice,[''r='',num2str(rval)],''Color'',''w'');');
   evalin('base','set(fig_2Dslice, ''Name'',[field_name,'' \theta-\phi R='',num2str(rval)]);');
      
   %Set YDATA in 1D plots (plot selection done elsewhere) 
   %Can put horizontal versus vertical line graph options here.
   evalin('base','set(plot_1D,''YData'',squeeze(MAS_DATA_CUT(pidx,:,ridx)))');
   evalin('base','set(axis_1Dplot,''YLim'',[min(squeeze(MAS_DATA_CUT(pidx,:,ridx))) max(squeeze(MAS_DATA_CUT(pidx,:,ridx)))] );');
   evalin('base','rval = rvec_cut (ridx);');
   evalin('base','title(axis_1Dplot,[''r='',num2str(rval),''      \phi='',num2str(pval)],''Color'',''k'');');
   evalin('base','set(fig_1Dplot,''Name'',[field_name,'' R='',num2str(rval),'' P='',num2str(pval)]);');
         
end
