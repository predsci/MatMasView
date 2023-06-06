function wrhdf_3d(filename,data,xscales,yscales,zscales,data_type)

  nx = length(xscales);
  ny = length(yscales);
  nz = length(zscales);

  import matlab.io.hdf4.*
  sdID = sd.start(filename,'create');
  sd.setFillMode(sdID,'nofill');
  sdsID = sd.create(sdID,'Data-Set-2',data_type,[nz ny nx]);
  dimID = sd.getDimID(sdsID,0);
  sd.setDimScale(dimID,zscales);
  dimID = sd.getDimID(sdsID,1);
  sd.setDimScale(dimID,yscales);
  dimID = sd.getDimID(sdsID,2);
  sd.setDimScale(dimID,xscales);
  sd.writeData(sdsID,[0 0 0],permute(data,[3 2 1]));
  sd.endAccess(sdsID);
  sd.close(sdID);
   if(strcmp(data_type,'single') || strcmp(data_type,'double'))
%     %Now convert hdf4 to be compatible with old DFSD hdf4 tools:
     syscallstr=['export LD_LIBRARY_PATH="/usr/local/psi_tools/lib:${LD_LIBRARY_PATH}:/usr/local/psi_external_libraries/hdf5/lib:/usr/local/psi_external_libraries/hdf4/lib"; export PATH=${PATH}:/usr/local/psi_tools/gen_tools/bin; hdfsd2hdf ',filename,' ',filename];
     system(syscallstr);
   end
end
