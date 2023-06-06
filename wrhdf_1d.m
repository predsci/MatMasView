function wrhdf_1d(filename,data,xscales,data_type)

  nx = length(xscales);

  import matlab.io.hdf4.*
  sdID = sd.start(filename,'create');
  sd.setFillMode(sdID,'nofill');
  sdsID = sd.create(sdID,'Data-Set-2',data_type,nx);
  dimID = sd.getDimID(sdsID,0);
  sd.setDimScale(dimID,xscales);
  sd.writeData(sdsID,[0],data);
  sd.endAccess(sdsID);
  sd.close(sdID);
   if(strcmp(data_type,'single') || strcmp(data_type,'double'))
%     %Now convert hdf4 to be compatible with old DFSD hdf4 tools:
     syscallstr=['export LD_LIBRARY_PATH="/usr/local/psi_tools/lib:${LD_LIBRARY_PATH}:/usr/local/psi_external_libraries/hdf5/lib:/usr/local/psi_external_libraries/hdf4/lib"; export PATH=${PATH}:/usr/local/psi_tools/gen_tools/bin; hdfsd2hdf ',filename,' ',filename];
     system(syscallstr);
   end

end
