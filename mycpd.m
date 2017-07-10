function [out,in]=mycpd(samples,timesteps)
	%samples - how many samples do I generate? 
	%timestep - and for how many timesteps? 
	%ex. mycpd(1000,100);
	filename='cpdheat_fresh.inp'; %CPD input file
	myFile=fileread(filename);
	Nexp=samples;
	Out=[];
	In=[];
	for i=1:Nexp
	  VAL1=3.1*1E15*myrand(0.3); %pre-exponential factor
	  VAL2=61200*myrand(0.5);    %activation energy
	  VAL3=8100*myrand(0.5);     %standard deviation
	  fileID=fopen('cpdheat.inp','w'); 
	  newFile=strrep(strrep(strrep(myFile,'VAL1',num2str(VAL1,'%.1e%15d')),'VAL2',num2str(VAL2)),'VAL3',num2str(VAL3));  
	  fprintf(fileID,'%s',newFile);
	  fclose(fileID);
	  system('./cpd');  %Path to the executable
	  k=dlmread('cpdheat.out1',' ',1);
	  size(k)
	  if size(k,1)>=4000  %Discards results where the convergence failed
	      In(size(In,1)+1,:)=k(1:timesteps,19);
	      Out(size(Out,1)+1,:)=log([VAL1,VAL2,VAL3]);
	  end
	  
	end
	out=Out;
	in=In;
end
