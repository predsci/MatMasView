%MAT_MAS_VIEW_GEN:
%Visualize a 3D MAS output (hdf4)
%Predictive Science Inc.

close all;
clear all;

global ridx pidx tidx;
global rval pval tval r1dval t1dval p1dval;
global sliceid draw_cut_lines shading_str id set_clim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   INPUT PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_tag     = 'run_name';  

datasetdir = 'PATH-TO-DATA';

field_vec   = {'vr' 'rho' 't' 'vt' 'vp' 'br' 'bt' 'bp' 'jr' 'jt' 'jp'};
idx_start   = 3;  %Starting file index.
idx_end     = 3;  %Ending file index.

show_movie = 1;   %Show movie while saving.
save_movie = 1;   %Save pngs of movie.
pause_first_frame=0; %Pause after showing first frame to allow positioning b4 mov.

showmesh   = 0;     %Show grid on surfaces
fsize      = 18;    %Font size of figures (sometimes only shown at saves/prints)
alphaval   = 0.6;   %Transparacy value for non-selected 3D slices
skp        = 1;     %Data skip for slices in 3D plot.
show_init_rot  = 0; %Show initial rotation movie for 3D view.
draw_cut_lines = 1; %Draw cut lines on higher dimenion plots.
init_view='P';
shading_str='interp';   % or 'flat';
helio_scaling=0;   %Scale fields by thier solar wind trends (R/R0, etc).
R0 = 214.9395;  %R0 for the helio scaling (1AU in MAS units).
log_scale=0;

tdview=[-20 20];       
tdzoom=0.75;

rval=1.03;
tval=1.20;
pval=5.9062;

%Set domain to plot (set to 0 to do full domain):
r_cut_domain_phys = [1 1.1];
t_cut_domain_phys = [1.15 1.45];
p_cut_domain_phys = [5.85 5.95];

%Set clim (overwritten if set_clim=1):
set_clim = 0;     %Set color range to be fixed at min/max (1) or manual (0)
cmin=-4;
cmax=4;

%Select if want plot values to be stretched to uniform (1) or not (0):
r_uniform = 0;
t_uniform = 0;
p_uniform = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_n = idx_end-idx_start+1;
id=idx_start;

scrsz = get(0,'ScreenSize');
loadpsipals;

%Start field loop:
for k=1:length(field_vec)

close all;

field = field_vec{k};

file_name = [field num2str(idx_start,'%03d') '.hdf']; %Name of first data file [file_name].hdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap_val_n=256;
%Set field name and units based on filename:
if(file_name(1) == 't') 
   field_name  = 'Temperature (MK)';
   field_units = 'MK';
   field_scale = 28.0706672;
   cmap = hot(cmap_val_n);
elseif(strcmp(file_name(1:2),'br')) 
   field_name  = 'Magnetic Field (r)';
   if(log_scale==1)
     field_units = 'log_1_0 Gauss';
   else
     field_units = 'Gauss';
   end
   field_scale = 2.20689140;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1:2),'bt')) 
   field_name  = 'Magnetic Field (\theta)';
   field_units = 'Gauss';
   field_scale = 2.20689140;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1:2),'bp')) 
   field_name  = 'Magnetic Field (\phi)';
   field_units = 'Gauss';
   field_scale = 2.20689140;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1:2),'ar')) 
   field_name  = 'Vector Potential (r)';
   field_units = 'um...';
   field_scale = 1;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1:2),'at')) 
   field_name  = 'Vector Potential (\theta)';
   field_units = 'um...';
   field_scale = 1;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1:2),'ap')) 
   field_name  = 'Vector Potential (\phi)';
   field_units = 'um...';
   field_scale = 1;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1:2),'vr')) 
   field_name  = 'Velocity (r)';
   field_units = 'km/s';
   field_scale = 481.37106736;
   cmap = psi_pal_bluered;%jet(cmap_val_n);
elseif(strcmp(file_name(1:2),'vt')) 
   field_name  = 'Velocity (\theta)';
   field_units = 'km/s';
   field_scale = 481.37106736;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1:2),'vp')) 
   field_name  = 'Velocity (\phi)';
   field_units = 'km/s';
   field_scale = 481.37106736;
   cmap = psi_pal_bluered;%psi_pal_bluered;
elseif(strcmp(file_name(1:3),'rho')) 
   field_name  = 'Density';
   if(log_scale==1)
     field_units = 'log_1_0 g/cm^3';
   else
     field_units = 'g/cm^3';
   end
   field_scale = 1.67e-16;
   cmap = (bone(cmap_val_n));
elseif(strcmp(file_name(1),'p')) 
   field_name  = 'Pressure';
   field_units = 'dyn/cm^2';   
   field_scale = 0.38757170;
   cmap = psi_pal_banded;
elseif(strcmp(file_name(1),'j')) 
   field_name  = 'Current Density';
   field_units = 'stat-Amp/cm^2';   
   field_scale = 0.07564541;
   cmap = psi_pal_bluered;
elseif(strcmp(file_name(1),'A')) 
   field_name  = '|A|';
   field_units = 'Gauss-cm';   
   field_scale = 2.20689140;
   cmap = jet(cmap_val_n);
else
   field_name  = 'Unknown';
   field_units = 'Unknown';  
   field_scale = 1;
   cmap = jet(cmap_val_n);
end

%Read MAS data:
MAS_DATA  = double(hdfread([datasetdir file_name], '/Data-Set-2'));

%Get grid data for half-mesh:
MAS_grid_info = hdfinfo([datasetdir file_name]);
pvec = double(MAS_grid_info.SDS.Dims(1).Scale);
tvec = double(MAS_grid_info.SDS.Dims(2).Scale);
rvec = double(MAS_grid_info.SDS.Dims(3).Scale);

%Thin out data by skp factor:
MAS_DATA = MAS_DATA(1:skp:end,1:skp:end,1:skp:end);
rvec=rvec(1:skp:end);
tvec=tvec(1:skp:end);
pvec=pvec(1:skp:end);

%Get resolution of mesh:
nr = length(rvec);
nt = length(tvec);
np = length(pvec);

%Set internal index vectors:
nrin = 2:nr-1;
ntin = 2:nt-1;
npin = 2:np-1;   

if(length(r_cut_domain_phys)==1)
   r_cut_domain_phys = [rvec(1) rvec(end)];
end
if(length(t_cut_domain_phys)==1)
   t_cut_domain_phys = [tvec(1) tvec(end)];
end
if(length(p_cut_domain_phys)==1)
   p_cut_domain_phys = [pvec(1) pvec(end)];
end

%Set initial cut locations:
if(rval==-9999)
  rval = r_cut_domain_phys(1) + 0.5*(r_cut_domain_phys(2)-r_cut_domain_phys(1));
end
if(tval==-9999)
  tval = t_cut_domain_phys(1) + 0.5*(t_cut_domain_phys(2)-t_cut_domain_phys(1));
end
if(pval==-9999)
  pval = p_cut_domain_phys(1) + 0.5*(p_cut_domain_phys(2)-p_cut_domain_phys(1));%2.8377%2.5;%
end

%Set initial r,t,p values for cut locations:
r1dval = rval;
t1dval = tval;
p1dval = pval;
                
%Find index range of desired physical range:
r_cut_range = find(rvec>=r_cut_domain_phys(1),1,'first'):...
              find(rvec<=r_cut_domain_phys(2),1,'last');
t_cut_range = find(tvec>=t_cut_domain_phys(1),1,'first'):...
              find(tvec<=t_cut_domain_phys(2),1,'last');
p_cut_range = find(pvec>=p_cut_domain_phys(1),1,'first'):...
              find(pvec<=p_cut_domain_phys(2),1,'last');
          
%Get cooridinate value arrays for cuts:
rvec_cut = rvec(r_cut_range);
tvec_cut = tvec(t_cut_range);
pvec_cut = pvec(p_cut_range);             

%Display uniform spacing (stretch) if wanted:
if(r_uniform == 1)
    rvec_label_vals = rvec_cut;
    rvec_cut = linspace(rvec_cut(1),rvec_cut(end),length(rvec_cut));
end
if(t_uniform == 1)
    tvec_label_vals = tvec_cut;
    tvec_cut = linspace(tvec_cut(1),tvec_cut(end),length(tvec_cut));
end
if(p_uniform == 1)
    pvec_label_vals = pvec_cut;
    pvec_cut = linspace(pvec_cut(1),pvec_cut(end),length(pvec_cut));
end

%Get indices into cut arrays:
ridx = find(rvec_cut>=rval,1,'first');
tidx = find(tvec_cut>=tval,1,'first');
pidx = find(pvec_cut>=pval,1,'first');

rval = rvec_cut(ridx);
tval = tvec_cut(tidx);
pval = pvec_cut(pidx);

%Get resolution of cut cube:
nr_cut = length(rvec_cut);
nt_cut = length(tvec_cut);
np_cut = length(pvec_cut);

%Get 3D cut cube coordinate matrices:                    
[P,T,R] = ndgrid(pvec_cut,tvec_cut,rvec_cut);

%Convert coordinate matrices to Cartesean:
X = R.*sin(T).*cos(P);
Y = R.*sin(T).*sin(P);
Z = R.*cos(T);

%Extract portion wanted and get into proper units:
MAS_DATA_CUT = field_scale.*MAS_DATA(p_cut_range,t_cut_range,r_cut_range);

if(helio_scaling==1)
  if(strcmp(file_name(1:3),'rho') || strcmp(file_name(1),'b'))
    MAS_DATA_CUT=MAS_DATA_CUT.*(R./R0).^2;
  elseif(strcmp(file_name(1),'p') || strcmp(file_name(1),'j'))
    MAS_DATA_CUT=MAS_DATA_CUT.*(R./R0).^3;
  elseif(strcmp(file_name(1),'t'))
    MAS_DATA_CUT=MAS_DATA_CUT.*(R./R0);
  end
elseif(log_scale==1)
  MAS_DATA_CUT2=MAS_DATA_CUT;
  MAS_DATA_CUT(MAS_DATA_CUT2<0)=log10(-MAS_DATA_CUT2(MAS_DATA_CUT2<0));
  MAS_DATA_CUT(MAS_DATA_CUT2>0)=log10(MAS_DATA_CUT2(MAS_DATA_CUT2>0));   
  clear MAS_DATA_CUT2;
end

if(set_clim == 1)
   cmin = min(MAS_DATA_CUT(:));
   cmax = max(MAS_DATA_CUT(:));
   
   if((strcmp(file_name(1),'b') || strcmp(file_name(1),'v') || strcmp(file_name(1),'a'))) % || strcmp(file_name(1),'j'))
     cabsmax = max(abs(cmin),abs(cmax));
     cmin = -cabsmax;
     cmax = cabsmax;
   end

   if(strcmp(file_name(1:2),'vr'))
        cmin = 0;
        cmax = cabsmax;
   end
   
end

%%%%%%%%%%%%%%%%%%%PLOT_INIT_FRAME%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%PLOT CUTS IN 3D:
fig_cube   = figure;
set(fig_cube,'Position',[scrsz(4)/8 scrsz(4)/2 scrsz(4)/2 scrsz(4)/2],'NumberTitle','off','Color','k');          
set(fig_cube,'InvertHardCopy','off');


h_rcut = surface(squeeze(X(:,:,ridx)),squeeze(Y(:,:,ridx)),squeeze(Z(:,:,ridx)),...
         squeeze(MAS_DATA_CUT(:,:,ridx)),'FaceAlpha',1.0);
set(h_rcut,'ButtonDownFcn', {@startDragFncR});    

hold on;

 h_pcut = surface(squeeze(X(pidx,:,:)),squeeze(Y(pidx,:,:)),squeeze(Z(pidx,:,:)),...
         squeeze(MAS_DATA_CUT(pidx,:,:)),'FaceAlpha',alphaval);
 set(h_pcut,'ButtonDownFcn', {@startDragFncP});    
 
 
 h_tcut =surface(squeeze(X(:,tidx,:)),squeeze(Y(:,tidx,:)),squeeze(Z(:,tidx,:)),...
         squeeze(MAS_DATA_CUT(:,tidx,:)),'FaceAlpha',alphaval);
 set(h_tcut,'ButtonDownFcn', {@startDragFncT});  


hold off;

set(gca,'FontSize',fsize);  
xlabel('x','FontSize',fsize); ylabel('y','FontSize',fsize); zlabel('z','FontSize',fsize);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w');
if(showmesh==1)
    shading faceted;
else
    shading(gca,shading_str);
    set(h_rcut,'LineStyle','none');    
    set(h_tcut,'LineStyle','none');    
    set(h_pcut,'LineStyle','none');    
end
axis equal;
axis tight;
axis vis3d;

caxis([cmin cmax]);
cb = colorbar;
colormap(cmap);
cblabel([field_name ' (' field_units ')'],'FontSize',fsize);
set(cb,'FontSize',fsize);
grid on;
view(tdview);
camzoom(tdzoom);

if (idx_n>1)
  dim = [0.1 0.8 0.1 0.1];
  str = ['Frame: ' num2str(idx_start,'%03d')];
  plot3d_annot=annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','w','FontSize',20,'LineWidth',2,'EdgeColor','w');
end

if(show_init_rot==1)
%Rotation animation
idx=0;
for i=1:5:360
   camorbit(5,0,'camera')
   drawnow
   fnamemovie3d = ['movies/' run_tag '_' field '_' num2str(idx_start) '_3D_rotation2_',num2str(idx,'%06d'),'.png'];
   screen2png(fig_cube, fnamemovie3d,'-r96');
   idx=idx+1;
   pause(0.01)
end
end

%%%%%%%%%%%TWO-DIMENSIONAL PLOTS%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_2Dslice = figure;
%customPos(3, 1, 2);
%set(fig_2Dslice, 'Name',[field_name,' R-T P=',num2str(pval)],'NumberTitle','off','Color',[0.20 0.20 0.25],...
% 'Position',[(scrsz(4)/2)+scrsz(4)/7 scrsz(4)/2 scrsz(4)/2 scrsz(4)/2]);
set(fig_2Dslice,'NumberTitle','off','Color','k',...
 'Position',[(scrsz(4)/2)+scrsz(4)/7 scrsz(4)/2 scrsz(4)/2 scrsz(4)/2]);
set(fig_2Dslice,'InvertHardCopy','off');
axis_slice=axes();

%%%%%%%%%%%ONE-DIMENSIONAL PLOTS%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_1Dplot = figure;
set(fig_1Dplot,'NumberTitle','off','Color',[0.90 0.90 0.95],'Position',...
    [(scrsz(4)/2)+(scrsz(4)/2)+scrsz(4)/6 scrsz(4)/2 scrsz(4)/2 scrsz(4)/2]);
set(fig_1Dplot,'InvertHardCopy','off');
axis_1Dplot=axes();

%Start drag function on initial view to plot first frame:
if(init_view=='R')
  startDragFncR;
elseif(init_view=='T')
  startDragFncT;
elseif(init_view=='P')
  startDragFncP;
end

if (idx_n>1)
    set(fig_cube,    'Name',[field_name, '  Frame: ' num2str(idx_start,'%03d')])
else
    set(fig_cube,    'Name',field_name)
end
str_tmp=get(plot1d_title,'String');
set(plot1d_title,'String',str_tmp);
set(fig_1Dplot,  'Name',str_tmp)

str_tmp=get(plot2d_title,'String');
set(plot2d_title,'String',str_tmp);
set(fig_2Dslice, 'Name',str_tmp);

%%%%%%%%%%%%%%%%%%%%%%%%START ANIMATION%%%%%%%%%%%%%%%%%%%%

if(show_movie==0)
  set(fig_1Dplot, 'Visible','off');
  set(fig_2Dslice, 'Visible','off');
  set(fig_cube, 'Visible','off');
end

disp('Starting animation...');

if(pause_first_frame==1)
  disp('Press enter to start animation');
  pause
end

%Loop over frames
for id = idx_start:idx_end
   
   disp(['Loading frame ',num2str(id)]);
   
   fnamemovie3d = ['movies/' run_tag '_' field '_' num2str(idx_start) '-' num2str(idx_end)  '_3D_',num2str(id,'%06d'),'.png'];
   fnamemovie2d = ['movies/' run_tag '_' field '_' num2str(idx_start) '-' num2str(idx_end)  '_2D_',num2str(id,'%06d'),'.png'];
   fnamemovie1d = ['movies/' run_tag '_' field '_' num2str(idx_start) '-' num2str(idx_end)  '_1D_',num2str(id,'%06d'),'.png'];

   %Load new data:
   file_name = [field num2str(id,'%03d') '.hdf'];
   MAS_DATA  = double(hdfread([datasetdir file_name], '/Data-Set-2'));

   %Thin out data by skp factor:
   MAS_DATA = MAS_DATA(1:skp:end,1:skp:end,1:skp:end);

  %Extract portion wanted and get into proper units:
   MAS_DATA_CUT = field_scale.*MAS_DATA(p_cut_range,t_cut_range,r_cut_range);
   
   clear MAS_DATA;
   
if(helio_scaling==1)
  if(strcmp(file_name(1:3),'rho') || strcmp(file_name(1),'b'))
    MAS_DATA_CUT=MAS_DATA_CUT.*(R./R0).^2;
  elseif(strcmp(file_name(1),'p') || strcmp(file_name(1),'j'))
    MAS_DATA_CUT=MAS_DATA_CUT.*(R./R0).^3;
  elseif(strcmp(file_name(1),'t'))
    MAS_DATA_CUT=MAS_DATA_CUT.*(R./R0);
  end
elseif(log_scale==1)
  MAS_DATA_CUT2=MAS_DATA_CUT;
  MAS_DATA_CUT(MAS_DATA_CUT2<0)=log10(-MAS_DATA_CUT2(MAS_DATA_CUT2<0));
  MAS_DATA_CUT(MAS_DATA_CUT2>0)=log10(MAS_DATA_CUT2(MAS_DATA_CUT2>0));   
  clear MAS_DATA_CUT2;
end
   
   
   %Set 3D cut frame:
   set(h_rcut,'CData',squeeze(MAS_DATA_CUT(:,:,ridx)));
   set(h_tcut,'CData',squeeze(MAS_DATA_CUT(:,tidx,:)));
   set(h_pcut,'CData',squeeze(MAS_DATA_CUT(pidx,:,:)));   
   
   if(strcmp(sliceid,'rt'))    
      %Set 2D cut plots CDATA (plot selection done in startdragfunc)
      set(plot_slice,'CData',fliplr(squeeze(MAS_DATA_CUT(pidx,:,:))'));           
      %Set YDATA in 1D plots (plot selection done elsewhere)
      set(plot_1D,'YData',squeeze(MAS_DATA_CUT(pidx,tidx,:)));
     % set(axis_1Dplot,'YLim',[min(squeeze(MAS_DATA_CUT(pidx,tidx,:))) max(squeeze(MAS_DATA_CUT(pidx,tidx,:)))] );
   elseif(strcmp(sliceid,'rp'))      
      %Set 2D cut plots CDATA (plot selection done in startdragfunc)
      set(plot_slice,'CData',squeeze(MAS_DATA_CUT(:,tidx,:)));           
      %Set YDATA in 1D plots (plot selection done elsewhere)
      set(plot_1D,'YData',squeeze(MAS_DATA_CUT(pidx,tidx,:)));
     % set(axis_1Dplot,'YLim',[min(squeeze(MAS_DATA_CUT(pidx,tidx,:))) max(squeeze(MAS_DATA_CUT(pidx,tidx,:)))] );
   elseif(strcmp(sliceid,'tp'))      
      %Set 2D cut plots CDATA (plot selection done in startdragfunc)
      set(plot_slice,'CData',squeeze(MAS_DATA_CUT(:,:,ridx)));           
      %Set YDATA in 1D plots (plot selection done elsewhere)
      set(plot_1D,'YData',squeeze(MAS_DATA_CUT(pidx,:,ridx)));
     % set(axis_1Dplot,'YLim',[min(squeeze(MAS_DATA_CUT(pidx,:,ridx))) max(squeeze(MAS_DATA_CUT(pidx,:,ridx)))] );
   end
     
   set(fig_cube,    'Name',[field_name, '  Frame: ' num2str(id,'%03d')])
   if(idx_n>1)
     set(plot3d_annot,'String',['Frame: ' num2str(id,'%03d')]);   
   end

   str_tmp=get(plot1d_title,'String');
   str_tmp(end-2:end)=num2str(id,'%03d');
   set(plot1d_title,'String',str_tmp);
   set(fig_1Dplot,  'Name',str_tmp)
   
   str_tmp=get(plot2d_title,'String');
   str_tmp(end-2:end)=num2str(id,'%03d');
   set(plot2d_title,'String',str_tmp);
   set(fig_2Dslice, 'Name',str_tmp);

   drawnow;

   if(save_movie==1)
     screen2png(fig_cube,   fnamemovie3d,'-r96');
     screen2png(fig_2Dslice,fnamemovie2d,'-r96');
     screen2png(fig_1Dplot, fnamemovie1d,'-r96');
   end

end

disp('Done!');

end

return


