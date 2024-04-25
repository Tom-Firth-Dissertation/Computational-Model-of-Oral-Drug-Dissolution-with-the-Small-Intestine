
niter = 5000;
pniter = 167; %parabolic iterations
rniter = pniter - 1;
nx = 400; % nodify these to a 1:4 ratio to speed up processing time, only takes a few 1000s of iterations then
ny = 100;
mx = 400;
my = 2*ny;
R = ny/2;
delt = 1;
delx = 1;
tau = 1*delt;
p = 22*pniter/1000; %inversely proportional to the rate of contraction
conc_total = zeros(101,5);
conc_level = 1;

% D2Q9 velocity set parameters
ndir = 9;
cssq = 1/3;
cx = [0, 1, -1, 0, 0, 1, -1, -1, 1];
cy = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
cssqinv = 1/cssq;
omega = 1/tau;
omomega = 1 - omega;
D = cssq*(tau - 0.5);

% Initialisation of the concentration and velocity field at time t = 0
rho = ones(mx, my); %density lattice formed
ux = zeros(mx,my);
uy = zeros(mx,my);

u_max = 0.05;
for y = 0:ny
    ux(:,y+ny/2) = u_max*(1-((y-R)/R)^2); %this subsitutes our bound going from 0 to 2R as oppose to -R to R
end
%Concentration
x_0 = 1*nx/4; %100
y_0 = my/2; %100
alpha = ny/4;
alpha2 = floor(ny/8);
C = zeros(mx, my); %density lattice formed

%%%%%%%%% CIRCLE CONC SHAPE %%%%%%%%%%

for a = 1:mx
    for b = 1:my
        distance = (a-x_0)^2 + (b-y_0)^2; 
        if distance < (alpha + 1)^2 %should encompass circle of radius alpha
            C(a,b) = 1;
        end
    end
end

%%%%%%%% DIAMOND CONC SHAPE %%%%%%%%%

% for x = 1:nx
%     for y = 1:my
%         if y>-x+100-y_0/4+x_0 && y>x+100-y_0/4-x_0  && y<x+100+y_0/4-x_0 && y<-x+100+y_0/4+x_0
%             C(x,y) = 1;
%         end
%     end
% end

%%%%%%%%%% OVAL CONC SHAPE %%%%%%%%%%

% for a = 1:nx
%     for b = 1:my
%         distance = (a-x_0)^2 + (b-y_0)^2; 
%         if distance < (alpha2 + 1)^2 %should encompass circle of radius alpha
%             if a>x_0
%                 C(a+21,b) = 1;
%             else
%                 C(a-21,b) = 1;
%             end
%         end
%         if x_0-21<a && a<x_0+22 && y_0-13<b && b<y_0+13
%             C(a,b) = 1;
%         end
%     end
% end

% Initialisation of Boundary Peristalsis
Boundaries = zeros(nx,my);
Boundaries(:,150) = 1;
Boundaries(:,50) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reynolds = u_max*2*alpha/(cssq^2*(tau-0.5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation of the particle distribution functions for fluid and
% concentration
feq = zeros(mx, my, ndir); 
geq = zeros(mx, my, ndir); 
for k = 1:ndir
    cdotu = cx(k).*ux + cy(k).*uy; %changing ux, uy to nxm matrices caused NaN errors which had to be resolved by.* and ./
    udotu = ux.^2 + uy.^2;
    feq(:, :, k) = w(k).*rho.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);
    geq(:, :, k) = w(k).*C.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);
    
end
f = feq;
fcol = zeros(mx, my, ndir);
g = geq;
gcol = zeros(mx, my, ndir);




% Simulation loop
fprintf('Starting simulation \n');
tic
for t = 1:niter
    if mod(t,1000)==0
        for y = 0:ny
            ux(:,y+ny/2) = u_max*(1-((y-R)/R)^2); %this subsitutes our bound going from 0 to 2R as oppose to -R to R
        end
    end
    % Collision
    fcol = omomega.*f + omega.*feq;
    %this is the BGK collision equation 8.28 in book
    %page 306

    % Streaming - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for j = 1:my
            for i = 1:mx
                xstreamed = mod(i + cx(k), mx);
                if xstreamed == 0
                    xstreamed = mx;
                end
                ystreamed = mod(j + cy(k), my);
                if ystreamed == 0
                    ystreamed = my;
                end
                f(xstreamed, ystreamed, k) = fcol(i, j, k);
            end
        end
    end

% Boundary conditions - this is where we implement
% g(X_b,t+delt) = -g*(X_b,t) the BB for fluids (ABB for chemicals)
% Info about moving walls pg. 180/200 LBM P&P
    uw = [0,50/pniter]; %speed of movement, constant but may depend of t!
    opp = [1,3,2,5,4,8,9,6,7];
    for k = 1:ndir
        f(nx,50:150,k) = f(nx-1,50:150,k); %apply Neumann to outlet x = nx
    end
    %Applied Zou&He Boundary Conditions for inlet (west side conditions
    %from thesis)
    rho_w = ones(nx,my);
    for j = 1:ny-1
        rho_w(:,j+ny/2) = (f(:,j+ny/2,1)+f(:,j+ny/2,4)+f(:,j+ny/2,5)+2*(f(:,j+ny/2,3)+f(:,j+ny/2,7)+f(:,j+ny/2,8)))/(1-ux(1,j+ny/2));
        f(1,j+ny/2,6) = f(1,j+ny/2,8)-(f(1,j+ny/2,4)-f(1,j+ny/2,5))/2 + (rho_w(1,j+ny/2)*ux(1,j+ny/2))/6;
        f(1,j+ny/2,9) = f(1,j+ny/2,7)+(f(1,j+ny/2,4)-f(1,j+ny/2,5))/2 + (rho_w(1,j+ny/2)*ux(1,j+ny/2))/6;
        f(1,j+ny/2,2) = f(1,j+ny/2,3)+2*rho_w(1,j+ny/2)*ux(1,j+ny/2)/3;
    end

    
    

    % Macroscopic variables
     %sum of density (rho)

    rho = f(:, :, 1) + f(:, :, 2) + f(:, :, 3) + f(:, :, 4) + f(:, :, 5)...
       + f(:, :, 6) + f(:, :, 7) + f(:, :, 8) + f(:, :, 9);
    %new functions for velocity!
    ux = (f(:,:,2)+f(:,:,6)+f(:,:,9)-f(:,:,3)-f(:,:,7)-f(:,:,8))./rho; 
    %note uy velocity should be <10^(-6) - nelgible
    uy = (f(:,:,4)+f(:,:,6)+f(:,:,7)-f(:,:,5)-f(:,:,8)-f(:,:,9))./rho;

    % Equilibrium distribution function
    for k = 1:ndir
        cdotu = cx(k).*ux + cy(k).*uy; %changing ux, uy to nxm matrices caused NaN errors which had to be resolved by.* and ./
        udotu = ux.^2 + uy.^2;
        feq(:, :, k) = w(k).*rho.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);
    end
    


    % Reassign Boundary


    %New Method
    Boundaries = zeros(nx,my);
    for i = 1:nx
        if mod(t,2*pniter) == 2*pniter-1
            Boundaries(:,150) = 1;
            Boundaries(:,50) = 1;
        elseif mod(floor(t/pniter),2) == 0 %even so Expand!
            Boundaries(i,floor(abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))-(t-floor(t/pniter)*pniter)/p)+150) = 1; %t value here must be between 0 and 1000
            Boundaries(i,ceil(abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))-(t-floor(t/pniter)*pniter)/p)+150) = 1;
            Boundaries(i,floor(-abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))+(t-floor(t/pniter)*pniter)/p)+50) = 1;
            Boundaries(i,ceil(-abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))+(t-floor(t/pniter)*pniter)/p)+50) = 1;
        elseif mod(floor(t/pniter),2) == 1 %odd so contract
            Boundaries(i,floor(abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))-(ceil(t/rniter)*rniter-t)/p)+150) = 1; %again t val between 0 and 1000 but inversed
            Boundaries(i,ceil(abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))-(ceil(t/rniter)*rniter-t)/p)+150) = 1;
            Boundaries(i,floor(-abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))+(ceil(t/rniter)*rniter-t)/p)+50) = 1;
            Boundaries(i,ceil(-abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))+(ceil(t/rniter)*rniter-t)/p)+50) = 1;
        end
    end

    x = 0:400;
    % Post Sim Code
    if mod(t, 250) == 0 || t == 1
        ConcMax = max(C,[],"all");
        fprintf('Iteration: %d, Time: %f \n', t, toc);
        % figure;imagesc(C.');colorbar;xlim([0,nx]);ylim([0,ny]);axis('equal');
        figure
        tiledlayout(2,1)
        nexttile;
        hold on
        imagesc(ux.');colorbar;colormap("hot");clim([0 u_max+0.001]);xlim([0,nx]);ylim([0,my]);axis('equal');title('Poiseuille Flow Model')
        if mod(t,2*pniter) == 2*pniter-1
            plot(50,'b');
            plot(150,'b');
        elseif mod(floor(t/pniter),2) == 0 %even so Expand!
            plot(x,abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))-(t-floor(t/pniter)*pniter)/p+150,('b'));
            plot(x,-abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))+(t-floor(t/pniter)*pniter)/p+50,('b'));
        elseif mod(floor(t/pniter),2) == 1 %odd so contract
            plot(x,abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))-(ceil(t/rniter)*rniter-t)/p+150,('b'));
            plot(x,-abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))+(ceil(t/rniter)*rniter-t)/p+50,('b'));
        end

        hold off

        nexttile;
        hold on
        imagesc(C.');colorbar;colormap("hot");axis('equal');clim([0 ConcMax]);xlim([0,nx]);ylim([0,my]);
        % title('Drug Concentration Diffusion Model')

        if mod(t,2*pniter) == 2*pniter-1
            plot(50,'b');
            plot(150,'b');
        elseif mod(floor(t/pniter),2) == 0 %even so Expand!
            plot(x,abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))-(t-floor(t/pniter)*pniter)/p+150,('b'));
            plot(x,-abs((t-floor(t/pniter)*pniter)/p*sin(pi*x/100))+(t-floor(t/pniter)*pniter)/p+50,('b'));
        elseif mod(floor(t/pniter),2) == 1 %odd so contract
            plot(x,abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))-(ceil(t/rniter)*rniter-t)/p+150,('b'));
            plot(x,-abs((ceil(t/rniter)*rniter-t)/p*sin(pi*x/100))+(ceil(t/rniter)*rniter-t)/p+50,('b'));
        end

        hold off
        % save(sprintf('Conc%d.mat',t),'C');
        % save(sprintf('Velocity%d.mat',t),'ux');
    end

    if mod(t,50) == 0 || t == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
        %figure;histogram(C.',10,"BinLimits",[0,0.5]);title('Concentration Histogram')
        for x = 1:nx
            for y = 1:my
                if C(x,y) < 0.05 && C(x,y) > 0.01
                    conc_total(1 + t/50,1) = conc_total(1 + t/50,1) + 1;
                elseif C(x,y) < 0.10 && C(x,y) > 0.05
                    conc_total(t/50 + 1,2) = conc_total(1 + t/50,2) + 1;
                elseif C(x,y) < 0.15 && C(x,y) > 0.10
                    conc_total(1 + t/50,3) = conc_total(1 + t/50,3) + 1;
                elseif C(x,y) < 0.20 && C(x,y) > 0.15
                    conc_total(1 + t/50,4) = conc_total(1 + t/50,4) + 1;
                elseif C(x,y) < 0.25 && C(x,y) > 0.20
                    conc_total(1 + t/50,5) = conc_total(1 + t/50,5) + 1;
                end
            end
        end
        ConcMax = max(C,[],"all");
        conc_level = [conc_level,ConcMax];
    end
%Of Flow Simulation Loop Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Simulation loop
        % Collision - remains the same
    if t > 0
        gcol = omomega.*g + omega.*geq; 
        %omega =1 so so indepent of g ^
        %this is the BGK collision equation 8.28 in book
        %page 306
    
        % Streaming - Explicit version -  remains the same
        % This streaming implementation automatically applies periodic boundary
        % conditions in all edges of the computational domain.
        for k = 1:ndir
            for j = 1:my
                for i = 1:mx
                    xstreamed = mod(i + cx(k), mx);
                    if xstreamed == 0
                        xstreamed = mx;
                    end
                    ystreamed = mod(j + cy(k), my);
                    if ystreamed == 0
                        ystreamed = my;
                    end
                    g(xstreamed, ystreamed, k) = gcol(i, j, k);
                end
            end
        end
    
    % Boundary conditions - this is where we implement
    % g(X_b,t+delt) = -g*(X_b,t) the BB for fluids (ABB for chemicals)
    % Info about moving walls pg. 180/200 LBM P&P
        for a = 1:nx
            for b = 1:my
                for k = 1:ndir
                    if Boundaries(a,b) == 1 && b >= 100
                        %g(a,b,k) = g(a,b-1,k); %apply Neumann to top wall
                        g(a,b,opp(k)) = gcol(a,b,k)-2*w(k)*C(a,b)*(uw(1)*cx(k)+uw(2)*cy(k))/cssq;
                        % g(a,b,k) = 0;
                    end
                    if Boundaries(a,b) == 1 && b <= 100
                        %g(a,b,k) = g(a,b+1,k); %apply Neumann to bottom wall 
                        g(a,b,opp(k)) = gcol(a,b,k)-2*w(k)*C(a,b)*(uw(1)*cx(k)+uw(2)*cy(k))/cssq;
                        % g(a,b,k) = 0; %dirichlet
                    end
                end
            end 
        end
        for k = 1:ndir
            g(mx,1+ny/2:3*ny/2,k) = g(mx-1,1+ny/2:3*ny/2,k); %apply Neumann to outlet x = nx
            g(1,1+ny/2:3*ny/2,k) = g(2,1+ny/2:3*ny/2,k); %apply Neumann to inlet x = 1
        end
        
    
        % Macroscopic variables
         %sum of concentration
        C = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5) + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);
    
    
        % Equilibrium distribution function
        for k = 1:ndir
            cdotu = cx(k).*ux + cy(k).*uy;
            udotu = ux.^2 + uy.^2;
            feq(:, :, k) = w(k).*rho.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);
            geq(:, :, k) = w(k).*C.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);
        end
    end
end

Max_Oval_SegmentationRe60 = conc_level;
Total_Oval_SegmentationRe60 = conc_total;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on
plot(Max_Oval_SegmentationRe75(1:20),'-.b');
plot(Max_Oval_SegmentationRe30(1:20),'-.g');
plot(Max_Oval_SegmentationRe60(1:20),'-.r');
xlabel('Iterations (x50)');
ylabel('Max Concentration');
% title('Re 30 - All Shapes');
legend('Reynolds 7.5','Reynolds 30','Reynolds 60','Location','northeast');
%title(legend,'Drug Shapes');
axis('tight')
hold off