clear; clc; 
%close all;

% Setup MC
wantfig = 0;
turntime = 30;
returntime = 90;
returntime2 = 150;
mcRUN = 10;

TIME = 230;
nag = 4;
nmod = 3;
adj = [1    0   0   1;
    0   1   0   1;
    0   0   1   1;
    1   1   1   1];
% adj = eye(nag);
% adj = ones(nag,nag);

ts = 5;
w = deg2rad(1);
A(:,:,1) = [1 ts 0 0; 0 1 0 0; 0 0 1 ts; 0 0 0 1];
% Aestim(:,:,1) = A(:,:,1);
% Aestim(:,:,2) = A(:,:,1);
% Aestim(:,:,3) = A(:,:,1);
A(:,:,2) = [1 sin(w*ts)/w 0 -(1-cos(w*ts))/w; 0 cos(w*ts) 0 -sin(w*ts); 0 (1-cos(w*ts))/w 1 sin(w*ts)/w; 0 sin(w*ts) 0 cos(w*ts)];
w = -w;
A(:,:,3) = [1 sin(w*ts)/w 0 -(1-cos(w*ts))/w; 0 cos(w*ts) 0 -sin(w*ts); 0 (1-cos(w*ts))/w 1 sin(w*ts)/w; 0 sin(w*ts) 0 cos(w*ts)];

q1 = 0.01;

B = [0.5*ts*ts 0; ts 0 ;0 0.5*ts*ts; 0 ts];
H = [1 0 0 0; 0 0 1 0];

Q{1} = q1.*eye(2);     %process noise covar
Q{2} = 10.*eye(2);     %process noise covar
Q{3} = 10.*eye(2);     %process noise covar

R = 10000*eye(2);

a = 0.9;
b = 0.95;
c = 0.95;
PI=[a (1-a)/2 (1-a)/2;...
    (1-b)/2 b (1-b)/2
    (1-c)/2 (1-c)/2 c];  % Markov transition matrix

tic
for mc = 1:mcRUN%:2e1
%     seed = randi(40000,1);
%     rng(seed);
    seed = randi(40000,1);
    % seed = 2236;%errorpartcon
    % seed = 38563;%goodfullcon
    rng(seed);
    clear P PP
    for i =1:nag
        mu(:,i) = [0.9 0.1 0]';
    end
    
    x(:,1) = [25000 120 12000 0]';
    
    XHAT = zeros(size(x,1),nag,TIME+1);
    PP = zeros(size(x,1),size(x,1),nag,nag,TIME+1);
    P = zeros(size(x,1),size(x,1),nag,nag,nmod);
    
    for i = 1:nag
        for j = 1:nag
            for k= 1:nmod
                xhat(:,i,k) = mvnrnd(x(:,1),P(:,:,i,i,k))';
                XHAT(:,i,1) = xhat(:,i,1);                      %mode fused result
            end
        end
    end
    
    for i = 1:nag
        for j = 1:nag
            for k= 1:nmod
                if(adj(i,j))
                    P(:,:,i,j,k) = P(:,:,i,j,k) + B*Q{k}*B' + 1e-12.*eye(4);
                end
%                 xhat(:,i,k) = A(:,:,k)*xhat(:,i,k);
                xhat(:,i,k) = xhat(:,i,k);
%                 XHAT(:,i,1) = xhat(:,i,1);                      %mode fused result
            end
        end
    end
    
    [a, b] = max(mu(:,1));
    sigma(1) = b;
    for i = 1:nag
        MHAT(i,1) = sigma(1);
    end
    
    % for i = 1:nag
    %     z(:,i) = H*x(:,1) + mvnrnd(zeros(2,1),R)';
    % end
    for t = 1:TIME
%         t
        x(:,t+1) = A(:,:,sigma(t))*x(:,t) + B*mvnrnd(zeros(2,1),Q{sigma(t)})';
        for i = 1:nag
            z(:,i) = H*x(:,t+1) + mvnrnd(zeros(2,1),R)';
        end
        
        [xhat, P, mu] = IMMnew(adj, xhat, P, z, A, B, H, Q, R, PI, mu);
        if(nnz(~isreal(mu)~=0))
            msg = 'Aiyyo.';
            error(msg);
        end
        for i = 1:nag
            for j = 1:nmod
                XHAT(:,i,t+1) = XHAT(:,i,t+1)+ (xhat(:,i,j)*mu(j,i));             %for cyberattacking
            end
        end
        
        
        for i=1:nag
            for j=1:nag
                dummy = zeros(size(P,1),size(P,1));
                if adj(i,j)
                    for s = 1:nmod
                        for vari = 1:nmod
                            dummy = dummy + mu(s,i)*mu(vari,j)*(P(:,:,i,j,s)+...
                                (xhat(:,i,s)-XHAT(:,i,t))*(xhat(:,j,vari)-XHAT(:,j,t))');%+...
                            %(xhat(:,i,s)-XHAT(:,i,t))*(xhat(:,i,s)-xhat(:,j,s))');
                        end
                    end
                end
                PP(:,:,i,j,t+1) = dummy;                            %for cyberattacking
            end
        end
        
        for i =1:nag
            [cyberattack, MHAT(i,t)] = max(mu(:,i));                 %for cyberattacking
            mucalc(:,t) = mu(:,1);
        end
        
        
        
        if t<turntime
            sigma(t+1)=1;
        elseif t<returntime
            sigma(t+1)=1;
        elseif t<returntime+30
            sigma(t+1)=2;
        elseif t<returntime2
            sigma(t+1)=2;
        else
            sigma(t+1)=1;
        end
        %     x(:,t+1) = A(:,:,sigma(t))*x(:,t) + B*mvnrnd(zeros(2,1),Q{sigma(t)})';
        %     x(:,t+1) = A(:,:,sigma(t))*x(:,t) + B*mvnrnd(zeros(2,1),Q{sigma(t)})';
        zplot(t,:) = z(:,1);
        %     for i = 1:nag
        %         z(:,i) = H*x(:,t+1) + mvnrnd(zeros(2,1),R)';
        %     end
    end
    seed
    this(1,:) = XHAT(1,1,:);
    this(2,:) = XHAT(3,1,:);
    for cases=2:TIME
        err(mc,cases) = ((this(1,cases-1)-x(1,cases)).^2+(this(2,cases-1)-x(3,cases)).^2).^.5;
    end
%     for cases = 1:TIME
%         val(cases,mc) = sum(diag(PP(:,:,1,1,cases)));
%     end
    mc
end
err2 = mean(err,1);
% val2 = mean(val,2);
% % %% plotting results
% fontname = 'Helvetica';
% set(0,'defaultaxesfontname',fontname);

% set(0,'defaulttextfontname',fontname);
% set(0,'defaulttextinterpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% plot(err2,'b')
% save dataset1 val1
%
if wantfig==1
    figure(1)
    this(1,:) = XHAT(1,1,:);
    this(2,:) = XHAT(3,1,:);
    hold on; plot((this(1,:)),(this(2,:)),'ro','MarkerSize',6)
    plot(x(1,:),x(3,:),'k','LineWidth',2);
    %     plot(zplot(:,1),zplot(:,2),'>')
    % figure(2)
    % for cases = 1:TIME
    % val(cases) = sum(diag(PP(:,:,1,1,cases)));
    % end
    % plot(val)
    axis equal
    xlabel({'East-position of aircraft'},'interpreter','latex')
    ylabel({'North-position of aircraft'},'interpreter','latex')
    plot(x(1,turntime),x(3,turntime),'b*','MarkerSize',10);
    legend({'Predicted aircraft location','Actual aircraft location','Waypoint'},'interpreter','latex')
    plot(x(1,returntime),x(3,returntime),'b*','MarkerSize',10);
    plot(x(1,returntime+30),x(3,returntime+30),'b*','MarkerSize',10);
    plot(x(1,returntime2),x(3,returntime2),'b*','MarkerSize',10);
    hold off;
    
    fontname = 'Helvetica';
    set(0,'defaultaxesfontname',fontname);
    set(0,'defaulttextfontname',fontname);
    set(0,'defaulttextinterpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    %
    figure(2)
    plot((mucalc(1,:)),'b*-')
    hold on
    plot((mucalc(2,:)),'s-','color',[0.3 0.7 0])
    plot((mucalc(3,:)),'mo-')
    %     plot(sigma-1,'r')
    axis([0 TIME -0.05 1.15])
    xlabel({'Time'},'interpreter','latex')
    ylabel({'Mode probability'},'interpreter','latex')
    h = fill([0 turntime turntime 0],[0 0 1 1],'g');
    set(h,'facealpha',.1)
    set(h,'EdgeColor','none')
    set(h,'FaceColor',[0 1 0])
    h = fill([turntime returntime returntime turntime ],[0 0 1 1],'b');
    set(h,'facealpha',.1)
    set(h,'EdgeColor','none')
    h = fill([ returntime returntime+30 returntime+30 returntime ],[0 0 1 1],'m');
    set(h,'facealpha',.1)
    set(h,'EdgeColor','none')
    legend({'$$\mu^4_1$$','$$\mu^4_2$$','$$\mu^4_3$$','Mode 2','Mode 1','Mode 3'},'interpreter','latex','Location','northwest','Orientation','horizontal')
    h = fill([returntime+30 returntime2 returntime2 returntime+30 ],[0 0 1 1],'g');
    set(h,'facealpha',.1)
    set(h,'EdgeColor','none')
    set(h,'FaceColor',[0 1 0])
    h = fill([returntime+30 TIME TIME returntime+30],[0 0 1 1],'b');
    set(h,'facealpha',.1)
    set(h,'EdgeColor','none')
    fontname = 'Helvetica';
    set(0,'defaultaxesfontname',fontname);
    set(0,'defaulttextfontname',fontname);
    set(0,'defaulttextinterpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
end
toc

figure(3)
plot(err2, 'r-')