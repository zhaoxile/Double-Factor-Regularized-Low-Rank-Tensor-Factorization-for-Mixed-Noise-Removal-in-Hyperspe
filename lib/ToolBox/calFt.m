function X=calFt(X,n1,n2) %%%calculate F^{*}X
        [N,m]=size(X);
        for ii=1:m
            temp2=ifft2(reshape(X(:,ii),n1,n2));
            X(:,ii)=reshape(temp2,N,1);
        end
end