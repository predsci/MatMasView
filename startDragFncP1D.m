function startDragFncP1D(hfig,events)
    
   %Recreate 1D plots:
   evalin('base','set(fig_1Dplot,''Name'',[field_name,'' T='',num2str(tval),'' P='',num2str(pval),''  Frame: '',num2str(id,''%03d'')]);');
   evalin('base','plot_1D = plot(axis_1Dplot,rvec_cut,squeeze(MAS_DATA_CUT(pidx,tidx,:)),''ko-'',''LineWidth'',1.5);');

   %evalin('base','axis_1Dplot = get(fig_1Dplot,''CurrentAxes'');');   
   evalin('base','xlabel(axis_1Dplot,''r'',''FontSize'',fsize); ');
   evalin('base','ylabel(axis_1Dplot,[field_name,'' ('',field_units,'')''],''FontSize'',fsize);');
   evalin('base','plot1d_title=title(axis_1Dplot,[''\phi='',num2str(pval),''      \theta='',num2str(tval)],''Color'',''k'');');
   evalin('base','set(axis_1Dplot,''FontSize'',fsize);');
   evalin('base','axis(axis_1Dplot,''tight'');');
   evalin('base','set(axis_1Dplot,''YLim'',[cmin cmax]);');  
%   evalin('base','set(axis_1Dplot,''YLim'',[min(squeeze(MAS_DATA_CUT(pidx,tidx,:))) max(squeeze(MAS_DATA_CUT(pidx,tidx,:)))] );');       
   
   evalin('base','grid(axis_1Dplot,''on'');');
  
   evalin('base', 'set(fig_1Dplot, ''WindowScrollWheelFcn'', {@draggingFncP});');
        
end
   
   