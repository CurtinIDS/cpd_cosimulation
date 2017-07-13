function values=cpd(inp)
	filename='cpdheat_fresh.inp'; %CPD input file
	myFile=fileread(filename);
      VAL1=inp(1,:);
      VAL2=inp(2,:);
      VAL3=inp(3,:);
	  fileID=fopen('cpdheat.inp','w'); 
	  newFile=strrep(strrep(strrep(myFile,'VAL1',num2str(VAL1,'%.1e%15d')),'VAL2',num2str(VAL2)),'VAL3',num2str(VAL3));  
	  fprintf(fileID,'%s',newFile);
	  fclose(fileID);
	  system('./cpd');  %Path to the executable
	  k=dlmread('cpdheat.out1',' ',1);
      values=zeros(size(k,1),2);
      values(:,1)=k(:,6);
      values(:,2)=k(:,19)*1000;
end