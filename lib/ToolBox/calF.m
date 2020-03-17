function X=calF(X,n1,n2) %%%calculate FX
        [N,m]=size(X);
        for ii=1:m
            temp1=fft2(reshape(X(:,ii),n1,n2));
            X(:,ii)=reshape(temp1,N,1);
        end
end
