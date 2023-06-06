function [data,xscales,yscales]=rdhdf_2d(filename)
  data = hdfread(filename,'/Data-Set-2');
  grid_info = hdfinfo(filename);
  xscales = grid_info.SDS.Dims(1).Scale;
  yscales = grid_info.SDS.Dims(2).Scale;
end
