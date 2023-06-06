function startDragFncT(hfig,events)
     
global sliceid draw_cut_lines;
   
sliceid = 'rp';
    
%Set 2D plot:    
evalin('base','[X_cut Y_cut] = pol2cart(squeeze(P(:,tidx,:)),squeeze(R(:,tidx,:)));');
evalin('base','plot_slice = pcolor(axis_slice,X_cut,Y_cut,squeeze(MAS_DATA_CUT(:,tidx,:)));');

evalin('base','set(fig_2Dslice, ''Name'',[field_name,'' R-P T='',num2str(tval)],''NumberTitle'',''off'',''Color'',[0.20 0.20 0.25]);');
%evalin('base','set(fig_2Dslice, ''Name'',[field_name,'' R-P T='',num2str(tval)],''NumberTitle'',''off'',''Color'',''w'');');


evalin('base','set(axis_slice,''Color'',[0.20 0.20 0.25],''XColor'',''w'',''YColor'',''w'',''ZColor'',''w'');');
%evalin('base','set(axis_slice,''Color'',''w'',''XColor'',''k'',''YColor'',''k'',''ZColor'',''k'');');

evalin('base','if(showmesh==1);shading(axis_slice, ''faceted'');else;shading(axis_slice,shading_str);end;');
evalin('base','colormap(axis_slice,cmap);');
evalin('base','cblabel(axis_slice,field_units,''FontSize'',fsize);');
evalin('base','caxis(axis_slice,[cmin cmax]);');
evalin('base','set(axis_slice,''FontSize'',fsize);');
evalin('base','axis(axis_slice, ''equal'');');
evalin('base','axis(axis_slice, ''tight'');');

evalin('base','cb = colorbar(''peer'',axis_slice,''FontSize'',fsize);');

%evalin('base','axis(axis_slice, [-30 30 -30 30]);');

evalin('base','xlabel(axis_slice,''r cos(\phi)'',''FontSize'',fsize);');
evalin('base','ylabel(axis_slice,''r sin(\phi)'',''FontSize'',fsize);');

evalin('base','plot2d_title=title(axis_slice,[''\theta='',num2str(tval) ''   Frame:'' num2str(id,''%03d'')],''Color'',''w'');');  
%evalin('base','plot2d_title=title(axis_slice,[''\theta='',num2str(tval)],''Color'',''k'');');  

%evalin('base','set(axis_slice,''Visible'',''off'');');

if(draw_cut_lines==1)
evalin('base','hold(axis_slice,''on'');');
evalin('base','X_line = squeeze(X_cut(pidx,:));');
evalin('base','Y_line = squeeze(Y_cut(pidx,:));');
evalin('base','axes(axis_slice);');
evalin('base','line(X_line,Y_line,''Color'',''w'',''LineWidth'',1.5);');
evalin('base','plot(axis_slice,X_line(ridx),Y_line(ridx),''wo'',''LineWidth'',1,''MarkerSize'',5,''MarkerFaceColor'',''c'');');
evalin('base','hold(axis_slice,''off'');');
end


evalin('base', 'set(fig_cube,    ''WindowScrollWheelFcn'', {@draggingFncT});');
evalin('base', 'set(h_rcut,   ''FaceAlpha'',alphaval);');    
evalin('base', 'set(h_tcut,   ''FaceAlpha'',1.0);'); 
evalin('base', 'set(h_pcut,   ''FaceAlpha'',alphaval);'); 
evalin('base', 'set(fig_2Dslice, ''WindowScrollWheelFcn'', {@draggingFncT});');

%Set 1D plots:    
startDragFncT1D;

evalin('base','drawnow;');

end
   
   