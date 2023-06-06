   function startDragFncP(hfig,events)
     
    global sliceid draw_cut_lines;
    
    sliceid = 'rt';
    
    %Recreate 2D plots:
    evalin('base','[X_cut Y_cut] = pol2cart(squeeze(T(pidx,:,:)),squeeze(R(pidx,:,:)));');
    %evalin('base','plot_slice = pcolor(axis_slice,X_cut,Y_cut,squeeze(MAS_DATA_CUT(pidx,:,:)));');
    evalin('base','plot_slice = pcolor(axis_slice,fliplr(Y_cut''),fliplr(X_cut''),fliplr(squeeze(MAS_DATA_CUT(pidx,:,:))''));');    
    evalin('base','set(fig_2Dslice, ''Name'',[field_name,'' R-T P='',num2str(pval),''  Frame: '',num2str(id,''%03d'')],''NumberTitle'',''off'',''Color'',[0.20 0.20 0.25]);');
    evalin('base','set(axis_slice,''Color'',[0.20 0.20 0.25],''XColor'',''w'',''YColor'',''w'',''ZColor'',''w'');');
    evalin('base','if(showmesh==1);shading(axis_slice, ''faceted'');else;shading(axis_slice,shading_str);end;');
    evalin('base','colormap(axis_slice,cmap);');
    evalin('base','cb = colorbar(''peer'',axis_slice,''FontSize'',fsize);');
    evalin('base','cblabel(axis_slice,field_units,''FontSize'',fsize);');
    evalin('base','caxis(axis_slice,[cmin cmax]);');
    evalin('base','set(axis_slice,''FontSize'',fsize);');
    evalin('base','axis(axis_slice, ''equal'');');
    evalin('base','axis(axis_slice, ''tight'');');
    evalin('base','ylabel(axis_slice,''r cos(\theta)'',''FontSize'',fsize);');
    evalin('base','xlabel(axis_slice,''r sin(\theta)'',''FontSize'',fsize);');
    evalin('base','plot2d_title=title(axis_slice,[''\phi='',num2str(pval)],''Color'',''w'');');

    if(draw_cut_lines==1)
      evalin('base','hold(axis_slice,''on'');');
      evalin('base','X_line = squeeze(X_cut(tidx,:));');
      evalin('base','Y_line = squeeze(Y_cut(tidx,:));');
      evalin('base','axes(axis_slice);');
      evalin('base','line(X_line,Y_line,''Color'',''w'',''LineWidth'',1.5);');
      evalin('base','plot(axis_slice,X_line(ridx),Y_line(ridx),''wo'',''LineWidth'',1,''MarkerSize'',5,''MarkerFaceColor'',''c'');');
      evalin('base','hold(axis_slice,''off'');');
    end
        
    evalin('base', 'set(fig_cube,   ''WindowScrollWheelFcn'', {@draggingFncP});');
    evalin('base', 'set(h_rcut,   ''FaceAlpha'',alphaval);');    
    evalin('base', 'set(h_tcut,   ''FaceAlpha'',alphaval);'); 
    evalin('base', 'set(h_pcut,   ''FaceAlpha'',1.0);'); 
    evalin('base', 'set(fig_2Dslice,   ''WindowScrollWheelFcn'', {@draggingFncP});');
    
    %Set 1D plots:    
    startDragFncP1D;

    evalin('base','drawnow;');
    
   end
   
   