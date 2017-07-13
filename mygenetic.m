function fitness=mygenetic(INP)
global exp1;
global output;
global myvars;
fitness=msd_cpd(NearestPoints(cpd(INP'),exp1),exp1);
%if fitness<0.
%output(size(output,1)+1,:)=fitness;
%myvars(size(myvars,1)+1,:)=INP;
end