function sgates=shuffle(signal,N)

signal=signal(:);
M=length(signal);
sgates=zeros(M,N);

for j=1:N
    s1=randperm(M);
    sgates(:,j)=signal(s1);
end

