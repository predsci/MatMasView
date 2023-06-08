function startDragFncR1D(hfig,events)
    
   %Recreate 1D plots:
   evalin('base','set(fig_1Dplot,''Name'',[field_name,'' R='',num2str(rval),'' P='',num2str(pval),''  Frame: '',num2str(id,''%03d'')]);');
   evalin('base','plot_1D = plot(axis_1Dplot,tvec_cut,squeeze(MAS_DATA_CUT(pidx,:,ridx)),''ko-'',''LineWidth'',1.5);');

   evalin('base','xlabel(axis_1Dplot,''\theta'',''FontSize'',fsize); ');
   evalin('base','ylabel(axis_1Dplot,[field_name,'' ('',field_units,'')''],''FontSize'',fsize);');
   evalin('base','plot1d_title=title(axis_1Dplot,[''r='',num2str(rval),''      \phi='',num2str(pval)],''Color'',''k'');');
   evalin('base','set(axis_1Dplot,''FontSize'',fsize);');
   evalin('base','axis(axis_1Dplot,''tight'');');
   evalin('base','set(axis_1Dplot,''YLim'',[cmin cmax]);');
   evalin('base','grid(axis_1Dplot,''on'');');

   evalin('base', 'set(fig_1Dplot, ''WindowScrollWheelFcn'', {@draggingFncR});');
        
end
   
   
