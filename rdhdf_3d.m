function [data,xscales,yscales,zscales]=rdhdf_3d(filename)
  data = hdfread(filename,'/Data-Set-2');
  grid_info = hdfinfo(filename);
  xscales = grid_info.SDS.Dims(1).Scale;
  yscales = grid_info.SDS.Dims(2).Scale;
  zscales = grid_info.SDS.Dims(3).Scale;
end
