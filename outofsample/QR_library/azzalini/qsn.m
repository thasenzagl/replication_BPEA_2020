function q=qsn(p,mu,sl,sr)
%quantiles of a split normal distribution

for jp = 1:length(p)
    if p(jp)<=.5
        temp = p(jp)*(sl+sr)/(2*sl);
        q(jp) = mu+sl*norminv(temp);    
    else
        temp = (1-p(jp))*(sl+sr)/(2*sr);
        q(jp) = mu+sr*norminv(1-temp); 
    end
        
end;
