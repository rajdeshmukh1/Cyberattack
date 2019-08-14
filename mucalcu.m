function mu = mucalcu(resid,P0,H,R,cbar,nmod)
for j = 1:nmod
    sigma = H*P0(:,:,j)*H'+R;
    alpa(j,1) = cbar(j)*((2*pi*det(sigma))^-(0.5));
%     if(~isreal(alpa(j,1)))
%         msg = 'Aiyyo.';
%         error(msg);
%     end
    bet(j,1) = -0.5*resid(:,j)'*pinv(sigma)*resid(:,j);    
end

clear dummy;

for j = 1:nmod    
    dummy = 0;
    for i = 1:nmod
        if(i~=j)
            dummy = dummy + (alpa(i)/alpa(j))*exp(bet(i)-bet(j));
        end
    end
    mu(j) = 1/(1+dummy);
end
