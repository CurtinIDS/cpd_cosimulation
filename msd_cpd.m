function out=msd_cpd(inp,exper)
     out=sqrt(sum((inp-exper(:,2)).^2));
end