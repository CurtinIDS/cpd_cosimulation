function [values,gradients]=Optimize(inp)
global exp1;
global output;
global vars_cpd;
values=msd_cpd(NearestPoints(cpd(inp),exp1),exp1);

gradients=zeros(size(inp,1),1);

for i=1:size(inp,1)
    %right derivative
    tmpinp=inp;
    if i==1
       tmpinp(i,:)=log(inp(1,:));
       val_right=tmpinp(1,:)*1.01;
       tmpinp(1,:)=exp(val_right);
    else
        val_right=tmpinp(i,:)*1.3;
        tmpinp(i,:)=val_right;
        
    end
    gradient_right=msd_cpd(NearestPoints(cpd(tmpinp),exp1),exp1);
    %left derivative
    tmpinp=inp;    
    if i==1
       tmpinp(i,:)=log(inp(1,:));
       val_left=tmpinp(1,:)*0.99;    
       tmpinp(1,:)=exp(val_right);
    else
        val_left=tmpinp(i,:)*0.7;   
        tmpinp(1,:)=val_right;
        
    end
    gradient_left=msd_cpd(NearestPoints(cpd(tmpinp),exp1),exp1);
    gradients(i,:)=((gradient_right-gradient_left))/(val_right-val_left);
   
end

output(size(output,1)+1,:)=values;
vars_cpd(size(vars_cpd,1)+1,:)=inp';
end