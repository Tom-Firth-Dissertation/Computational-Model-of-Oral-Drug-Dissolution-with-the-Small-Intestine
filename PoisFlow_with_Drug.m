
niter = 5001;                                                              
nx = 400; % nodify these to a 1:4 ratio to speed up processing time, only takes a few 1000s of iterations then
ny = 100;
R = ny/2;
delt = 1;
delx = 1;
tau = 1*delt;


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
rho = ones(nx, ny); %density lattice formed
ux = zeros(nx,ny);
uy = zeros(nx,ny);
u_max = 0.05;% how do we actually calculate this? Require the change in pressure from a segment in the intestine
for y = 2:(ny-1)
    ux(1,y) = u_max*(1-((y-R)/R)^2); %this subsitutes our bound going from 0 to 2R as oppose to -R to R
end
%Concentration
x_0 = 3*nx/4;
y_0 = ny/2;
alpha = ny/4; 
C = zeros(nx, ny); %density lattice formed
for a = 1:nx
    for b = 1:ny
        distance = (a-x_0)^2 + (b-y_0)^2; 
        if distance < (alpha + 1)^2 %should encompass circle of radius alpha
            C(a,b) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reynolds = u_max*2*alpha/(cssq^2*(tau-0.5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation of the particle distribution functions for fluid and
% concentration
feq = zeros(nx, ny, ndir); 
geq = zeros(nx, ny, ndir); 
for k = 1:ndir
    cdotu = cx(k).*ux + cy(k).*uy; %changing ux, uy to nxm matrices caused NaN errors which had to be resolved by.* and ./
    udotu = ux.^2 + uy.^2;
    feq(:, :, k) = w(k).*rho.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);
    geq(:, :, k) = w(k).*C.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);

end
f = feq;
fcol = zeros(nx, ny, ndir);
g = geq;
gcol = zeros(nx, ny, ndir);




% Simulation loop
fprintf('Starting simulation \n');
tic
for t = 1:niter
    % Collision
    fcol = omomega.*f + omega.*feq;
    %this is the BGK collision equation 8.28 in book
    %page 306

    % Streaming - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for j = 1:ny
            for i = 1:nx
                xstreamed = mod(i + cx(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(j + cy(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                f(xstreamed, ystreamed, k) = fcol(i, j, k);
            end
        end
    end

% Boundary conditions - this is where we implement
% g(X_b,t+delt) = -g*(X_b,t) the BB for fluids (ABB for chemicals)
% Info about moving walls pg. 180/200 LBM P&P
    opp = [1,3,2,5,4,8,9,6,7];
    for a = 1:nx 
        for k = 1:ndir
            f(a,ny,k) = fcol(a,ny,opp(k)); %BB for top wall
            f(a,1,k) = fcol(a,1,opp(k)); %BB for bottom wall
        end
    end
    for k = 1:ndir
        f(nx,1:ny,k) = f(nx-1,1:ny,k); %apply Neumann to outlet x = nx
    end
    %Applied Zou&He Boundary Conditions for inlet (west side conditions
    %from thesis)
    for j = 2:ny-1
        rho_w = (f(:,:,1)+f(:,:,4)+f(:,:,5)+2*(f(:,:,3)+f(:,:,7)+f(:,:,8)))/(1-ux(1,j));
        f(1,j,6) = f(1,j,8)-(f(1,j,4)-f(1,j,5))/2 + (rho_w(1,j)*ux(1,j))/6;
        f(1,j,9) = f(1,j,7)+(f(1,j,4)-f(1,j,5))/2 + (rho_w(1,j)*ux(1,j))/6;
        f(1,j,2) = f(1,j,3)+2*rho_w(1,j)*ux(1,j)/3;
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
        cdotu = cx(k).*ux + cy(k).*uy;
        udotu = ux.^2 + uy.^2;
        feq(:, :, k) = w(k).*rho.*(1 + cssqinv*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5.*cssqinv.*udotu);
    end
    
    % Processing Point %

    if mod(t, 1250) == 1
        fprintf('Iteration: %d, Time: %f \n', t, toc);
        figure
        tiledlayout(2,1)

        nexttile
        imagesc(ux.'); 
        colorbar;
        colormap("hot");
        xlim([0,nx]);
        ylim([0,ny]);
        axis('equal');
        title('Poiseuille Flow Model')
        
        nexttile
        imagesc(C.');
        colorbar;
        axis('equal');
        xlim([0,nx])
        ylim([0,ny])
        title('Drug Concentration Diffusion Model')

        % save(sprintf('Conc%d.mat',t),'C');
        % save(sprintf('Velocity%d.mat',t),'ux');        


        


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
            for j = 1:ny
                for i = 1:nx
                    xstreamed = mod(i + cx(k), nx);
                    if xstreamed == 0
                        xstreamed = nx;
                    end
                    ystreamed = mod(j + cy(k), ny);
                    if ystreamed == 0
                        ystreamed = ny;
                    end
                    g(xstreamed, ystreamed, k) = gcol(i, j, k);
                end
            end
        end
    
    % Boundary conditions - this is where we implement
    % g(X_b,t+delt) = -g*(X_b,t) the BB for fluids (ABB for chemicals)
    % Info about moving walls pg. 180/200 LBM P&P
        for a = 1:nx 
            for k = 1:ndir
                g(a,ny,k) = g(a,ny-1,k); %apply Neumann to top wall y = ny
                g(a,1,k) = g(a,2,k); %apply Neumann to bottom wall y = 1
                % g(a,ny,k) = 0; %apply Dirichlet to top wall y = ny
                % g(a,1,k) = 0; %apply Dirichlet to bottom wall y = 1
            end
        end
        for k = 1:ndir
            g(nx,1:ny,k) = g(nx-1,1:ny,k); %apply Neumann to outlet x = nx
            g(1,1:ny,k) = g(2,1:ny,k); %apply Neumann to inlet x = 1
        end
        
    
        % Macroscopic variables
         %sum of concentration
        C = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5) + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);
    
    
        % Equilibrium distribution function
        for k = 1:ndir
            cdotu = cx(k).*ux + cy(k).*uy;
            udotu = ux.^2 + uy.^2;
            feq(:, :, k) = w(k).*rho.*(1 + cssqinv.*cdotu + 0.5*cssqinv^2.*cdotu.^2 - 0.5*cssqinv.*udotu);
            geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);

        end
    end
end
