function [data,xscales]=rdhdf_1d(filename)
  data = hdfread(filename,'/Data-Set-2');
  grid_info = hdfinfo(filename);
  xscales = grid_info.SDS.Dims(1).Scale;
end