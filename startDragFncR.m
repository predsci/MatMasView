function startDragFncR(hfig,events)

global sliceid draw_cut_lines;

sliceid = 'tp';

%Set 2D plots:

evalin('base','plot_slice = pcolor(axis_slice,squeeze(P(:,:,ridx)),squeeze(T(:,:,ridx)),squeeze(MAS_DATA_CUT(:,:,ridx)));');

evalin('base','set(axis_slice,''YDir'',''Reverse'');');

evalin('base','set(fig_2Dslice, ''Name'',[field_name,'' P-T R='',num2str(rval),''  Frame: '',num2str(id,''%03d'')],''NumberTitle'',''off'',''Color'',[0.0 0.0 0.0]);');

evalin('base','set(axis_slice,''Color'',[0 0 0],''XColor'',''w'',''YColor'',''w'',''ZColor'',''w'');');

evalin('base','if(showmesh==1);shading(axis_slice, ''faceted'');else;shading(axis_slice,shading_str);end;');

evalin('base','colormap(axis_slice,cmap);');
evalin('base','cb = colorbar(''peer'',axis_slice,''FontSize'',fsize);');
evalin('base','cblabel(axis_slice,[field_name '' ('' field_units '')''],''FontSize'',fsize);');
evalin('base','caxis(axis_slice,[cmin cmax]);');
evalin('base','set(axis_slice,''FontSize'',fsize);');
evalin('base','axis(axis_slice, ''equal'');');
evalin('base','axis(axis_slice, ''tight'');');
evalin('base','xlabel(axis_slice,''\phi'',''FontSize'',fsize);');
evalin('base','ylabel(axis_slice,''\theta'',''FontSize'',fsize);');
evalin('base','plot2d_title=title(axis_slice,[''r='',num2str(rval)],''Color'',''w'');');

if(draw_cut_lines==1)
  evalin('base','hold(axis_slice,''on'');');
  evalin('base','X_line = squeeze(squeeze(P(pidx,:,ridx)));');
  evalin('base','Y_line = squeeze(squeeze(T(pidx,:,ridx)));');
  evalin('base','axes(axis_slice);');
  evalin('base','line(X_line,Y_line,''Color'',''w'',''LineWidth'',1.5);');
  evalin('base','plot(axis_slice,X_line(tidx),Y_line(tidx),''wo'',''LineWidth'',1,''MarkerSize'',5,''MarkerFaceColor'',''c'');');
  evalin('base','hold(axis_slice,''off'');');
end

evalin('base', 'set(fig_cube,   ''WindowScrollWheelFcn'', {@draggingFncR});');
evalin('base', 'set(h_rcut,   ''FaceAlpha'',1.0);');   
evalin('base', 'set(h_tcut,   ''FaceAlpha'',alphaval);'); 
evalin('base', 'set(h_pcut,   ''FaceAlpha'',alphaval);');     
evalin('base', 'set(fig_2Dslice,   ''WindowScrollWheelFcn'', {@draggingFncR});');

%Set 1D plots:   
startDragFncR1D;

evalin('base','drawnow;');

end
   
   
