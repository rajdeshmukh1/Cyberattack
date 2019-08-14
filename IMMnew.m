function [xhat, P, mu] = IMMnew(adj, xhat, P, z, A, B, H, Q, R, PI, mu,~)

nmod = size(PI,1);
nag = size(adj,1);

for i=1:nag
    for j=1:nmod
        for k=1:nmod
            MU(i,j,k) = PI(j,k)*mu(j,i);
        end
    end
end

for i = 1:nag
    for j =1:nmod
        temp(i,j)=0;
        for k =1:nmod
            temp(i,j) = temp(i,j) + MU(i,k,j);
        end
        for k=1:nmod
            MU(i,k,j) = MU(i,k,j)/temp(i,j);
        end
    end
end
% clear temp
cbar = temp;
clear temp;
xhat0(:,:,:) = zeros(size(xhat));

for i = 1:nag
    for j = 1:nmod
        for k = 1:nmod
            xhat0(:,i,j) = xhat0(:,i,j)+ (xhat(:,i,k)*MU(i,k,j));
        end
    end
end

P0(:,:,:,:,:) = zeros(size(P));
for i = 1:nag
    for j = 1:nag
%         if(adj(i,j))
          if(true)
            for k =1:nmod
                for s = 1:nmod
                    if i==j
                        P0(:,:,i,j,k) = P0(:,:,i,j,k)+ MU(i,s,k)*...
                            (P(:,:,i,j,s)+((xhat(:,i,s)-xhat0(:,i,k))*(xhat(:,j,s)-xhat0(:,j,k))'));
                    else
                        P0(:,:,i,j,k) = P0(:,:,i,j,k)+ MU(i,s,k)*MU(j,s,k)*...
                            (P(:,:,i,j,s)+((xhat(:,i,s)-xhat0(:,i,k))*(xhat(:,j,s)-xhat0(:,j,k))'));
                    end
                end
            end
        end
    end
end

%call KCF within loop for modes (xhat_i P_ij residuei) for mode l
for i = 1:nmod
    [xhat(:,:,i), P(:,:,:,:,i), resid(:,:,i),residP(:,:,:,:,i)] = KCFnew(adj,A(:,:,i),B,H,xhat0(:,:,i),P0(:,:,:,:,i),z,Q{i},R);
end

for i = 1:nag
    %         for j = 1:nmod
    %             lambdaz(i,j) = mvnpdf(resid(:,i,j),zeros(size(resid(:,i,j))),H*residP(:,:,i,i,j)*H'+R);
    %         lambdaz(i,j) = gausss(resid(:,i,j),zeros(size(resid(:,i,j))),H*P0(:,:,i,i,j)*H'+R);
    %         end
    %     mu(:,i) = mucalcu(resid(:,i,:),P(:,:,i,i,:),H,R,cbar(i,:),nmod);
    mu(:,i) = mucalcu(resid(:,i,:),residP(:,:,i,i,:),H,R,cbar(i,:),nmod);
end

% clear dummy;
% for i = 1:nag
%     for j = 1:nmod
%         dummy = 0;
%         for s = 1:nmod
%             dummy = dummy + PI(s,j)*mu(s,i);
%         end
%         mu(j,i) = lambdaz(i,j)*dummy;
%     end
% end
% 
% for i=1:nag
%     dum = sum(mu(:,i));
%     for j=1:nmod
%         mu(j,i) = mu(j,i)/dum;
%     end
%     clear dum;
% end