function out=NearestPoints(inp,exper)
    out=zeros(size(exper,1),1);
    for i=1:size(exper,1)
        [val,pos]=min(abs(inp(:,1)-exper(i,1)));
        out(i,:)=inp(pos,2);
    end
end