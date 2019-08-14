function [xhat, M , resid,residP] = KCFnew(adj,A,B,H,xhat,P,z,Q,R)

nag = size(adj,1);

% for i=1:nag
%     xhat(:,i) = A*xhat(:,i);
%     for j = 1:nag
%         P(:,:,i,j) = (A*P(:,:,i,j)*A' + B*Q*B');
%     end
% end

K(:,:,:) = zeros(size(xhat,1),size(z,1),nag);
C(:,:,:) = zeros(size(xhat,1),size(xhat,1),nag);

for i = 1:nag
    PriPi = zeros(size(P,1),size(P,1));
    for r=1:nag
        if(adj(i,r))
            PriPi = PriPi + (P(:,:,r,i) - P(:,:,i,i));
        end
    end
    clear E
    E = zeros(size(xhat,1),size(xhat,1));
    for r=1:nag
        for s=1:nag
            if(adj(r,i) && adj(s,i))
                E = E + P(:,:,i,i) -P(:,:,r,i) - P(:,:,i,s);
            end
            if(adj(r,i) && adj(s,i)) % && adj(r,s))        % (from i-j update)
                E = E + P(:,:,r,s);
            end
        end
    end
    
    clear G del
    del = R+H*P(:,:,i,i)*H';
    G = E - PriPi*H'*pinv(del)*H*PriPi';
    I = eye(size(xhat,1),size(xhat,1));
    if(det(G)~=0)
        C(:,:,i) = (P(:,:,i,i)*H'*pinv(del)*H - I)*PriPi'*pinv(G);
        %         C(:,:,i) = 0.1.*P(:,:,i,i)/(1+norm(P(:,:,i,i),'fro'));
    end
    K(:,:,i) = (P(:,:,i,i) + C(:,:,i)*PriPi)*H'*pinv(del);
    F(:,:,i) = I - K(:,:,i)*H;
    clear PriPi
end

for i=1:nag
    xshat(:,i) = xhat(:,i) + K(:,:,i)*(z(:,i) - H*xhat(:,i));
    for j=1:nag
        if (adj(i,j))
            xshat(:,i) = xshat(:,i)  + C(:,:,i)*(xhat(:,j) - xhat(:,i));
        end
    end
end

D = zeros(size(P,1),size(P,1),nag,nag);

for i = 1:nag
    for j = 1:nag
        if (true)       % i-j update
            for r = 1:nag
                for s = 1:nag
                    %                     if (adj(r,s) && adj(r,i)  && adj(s,j))
                    if(adj(i,r) && adj(s,j))
                        D(:,:,i,j) = D(:,:,i,j) + P(:,:,r,s);
                    end
                    %                     if (adj(r,i) && adj(s,j))
                    if(adj(i,r) && adj(s,j))
                        D(:,:,i,j) = D(:,:,i,j) + P(:,:,i,j);
                    end
                    %                     if (adj(r,i) && adj(s,j)  && adj(r,j))
                    if(adj(i,r) && adj(s,j))
                        D(:,:,i,j) = D(:,:,i,j) - P(:,:,r,j);
                    end
                    %                     if (adj(r,i)  && adj(s,j) && adj(i,s))
                    if(adj(i,r) && adj(s,j))
                        D(:,:,i,j) = D(:,:,i,j) - P(:,:,i,s);
                    end
                end
            end
        end
    end
end

M = zeros(size(P));
for i=1:nag
    for j=1:nag
        if (true)       % i-j update
            if i==j
                M(:,:,i,j) = F(:,:,i)*P(:,:,i,j)*F(:,:,j)' + C(:,:,i)*D(:,:,i,j)*C(:,:,j)' + K(:,:,i)*R*K(:,:,j)';
            else
                M(:,:,i,j) = F(:,:,i)*P(:,:,i,j)*F(:,:,j)' + C(:,:,i)*D(:,:,i,j)*C(:,:,j)';
            end
            for r=1:nag
                %                 if (adj(i,r) && adj(r,j))
                if (adj(i,r))
                    M(:,:,i,j) = M(:,:,i,j) + C(:,:,i)*(P(:,:,r,j)-P(:,:,i,j))*F(:,:,j)';
                end
            end
            for s=1:nag
                %                 if (adj(j,s) && adj(i,s))
                if (adj(j,s))
                    M(:,:,i,j) = M(:,:,i,j) + F(:,:,i)*(P(:,:,i,s)-P(:,:,i,j))*C(:,:,j)';
                end
            end
        end
    end
end

residP = P;
for i=1:nag
    resid(:,i) = z(:,i) - H*xhat(:,i);
end
for i=1:nag
    xhat(:,i) = A*xshat(:,i);
    for j = 1:nag
        if(true)        % i-j update
            M(:,:,i,j) = A*M(:,:,i,j)*A' + B*Q*B';
        end
    end
end
%
% for i=1:nag
%     resid(:,i) = z(:,i) - H*xhhat(:,i);
% end