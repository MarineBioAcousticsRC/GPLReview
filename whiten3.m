function [spc,qs,mu] = whiten3(sp,fac)

if nargin==1
 fac=1;end

[sz1,sz2]=size(sp);
qs=zeros(sz1,ceil(sz2/2)+2)+nan;
for i=1:sz1
    [ks]=base3(sp(i,:));
    qs(i,1:length(ks))=sp(i,ks);
end
sqs=sum(qs);
k=find(isfinite(sqs));
qs=qs(:,1:k(end));
[sz3,sz4]=size(qs);
mu = mean(qs')';
spc=sp-fac*mu*ones(1,sz2);
qs=qs-fac*mean(qs')'*ones(1,sz4);