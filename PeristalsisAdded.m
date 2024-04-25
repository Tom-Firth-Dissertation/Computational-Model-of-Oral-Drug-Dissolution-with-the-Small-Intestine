niter = 4000;
pniter = 1000; %parabolic iterations
nx = 400; % nodify these to a 1:4 ratio to speed up processing time, only takes a few 1000s of iterations then
ny = 100;
mx = nx;
my = 2*ny;
R = ny/2;
delt = 1;
delx = 1;
tau = 1*delt;
conc_total = zeros(101,5);
conc_level = 1;
expansion = 25; %scales diameter at maximum expansion to 150%

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

u_max = 0.05;% how do we actually calculate this? Require the change in pressure from a segment in the intestine
for y = 0:ny
    ux(:,y+ny/2) = u_max*(1-((y-R)/R)^2); %this subsitutes our bound going from 0 to 2R as oppose to -R to R
end


%Concentration
x_0 = 100; %100
y_0 = 100; %100
alpha = ny/4; %25

C = zeros(mx, my); %density lattice formed

%%%%%%%%% CIRCLE CONC SHAPE %%%%%%%%%%
%area = 2109

for a = 1:mx
    for b = 1:my
        distance = (a-x_0)^2 + (b-y_0)^2; 
        if distance < (alpha + 1)^2 %should encompass circle of radius alpha
            C(a,b) = 1;
        end
    end
end

%%%%%%%% DIAMOND CONC SHAPE %%%%%%%%%
%area = 2113

% l = 3.1; %similar area to circle!
% for x = 1:nx
%     for y = 1:my
%         if y>-x+100-y_0/l+x_0 && y>x+100-y_0/l-x_0  && y<x+100+y_0/l-x_0 && y<-x+100+y_0/l+x_0
%             C(x,y) = 1;
%         end
%     end
% end

%%%%%%%%%% OVAL CONC SHAPE %%%%%%%%%%
% area = 2095

% alpha2 = floor(ny/6.5); 
% h = 16;
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
%         if x_0-21<a && a<x_0+22 && y_0-h<b && b<y_0+h
%             C(a,b) = 1;
%         end
%     end
% end

%%%%%%%%%% TRIANGLE CONC SHAPE %%%%%%%%
%area = 2113


% a = 0.5;
% b = 40;
% h = 104;
% g = 80;
% for x = 1:nx
%     for y = 1:my
%         if y <= a*(x-b)+g && y >= -a*(x-b)+200-g && x <= h+b
%             C(x,y) = 1;
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% S T A R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s = 52; %Current Scalar of magnitude
% f = -y_0; %positions in y axis
% k = x_0 - 2*s/3; %positions x axis
% for x = 1:nx
%     for y = 1:my
%         if y+f <= sqrt(1/3)*(x-k) && y+f >= -sqrt(1/3)*(x-k) && x-k <= s || (x-k) >= s/3 && y+f >= sqrt(1/3)*((x-k)-4*s/3) && y+f <= -sqrt(1/3)*((x-k)-4*s/3)
%             C(x,y) = 1;
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  T E A R D R O P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r = 16;
% a = 6.5;
% 
% for x = 1:nx
%     for y = 1:my
%         if x <= 100 && (x-100)^2+(y-100)^2 <= r^2
%             C(x,y) = 1;
%         elseif x >= 100 && y-y_0 >= (1/a)*(x-x_0)-r && y-y_0 <= (-1/a)*(x-x_0)+r
%             C(x,y) =1;
%         end
%     end
% end



% Volume Calculation
Volume = sum(C,"all");
disp(Volume)


% Initialisation of Boundary Peristalsis
Boundaries = zeros(nx,my);
Boundaries(:,150) = 1;
Boundaries(:,50) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reynolds = u_max*ny/D;
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
    uw = [0,50/niter]; %speed of movement, constant but may depend of t!
    opp = [1,3,2,5,4,8,9,6,7];
    for a = 1:nx
        for b = 1:my
            for k = 1:ndir
                if Boundaries(a,b) == 1
                    f(a,b,opp(k)) = fcol(a,b,k)-2*w(k)*rho(a,b)*(uw(1)*cx(k)+uw(2)*cy(k))/cssq; %Moving BB for top wall y = 150
                    %f(a,50,k) = fcol(a,50,opp(k))-2*w(opp(k))*rho*(uw(0)*cx(opp(k))+uw(1)*cy(opp(k)))/cssq; %Moving BB for bottom wall y = 50
                end
            end
        end 
    end
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

    Boundaries = zeros(nx,my);
    for i = 1:nx
        if mod(floor(t/500),2) == 0 %even so Expand!
            Boundaries(i,floor(expansion/pniter*(t-floor(t/1000)*1000).*sin(pi*i/400))+150) = 1; %t value here must be between 0 and 1000
            Boundaries(i,ceil(expansion/pniter*(t-floor(t/1000)*1000).*sin(pi*i/400))+150) = 1;
            Boundaries(i,floor(-expansion*(t-floor(t/1000)*1000)/pniter.*sin(pi*i/400))+50) = 1;
            Boundaries(i,ceil(-expansion/pniter*(t-floor(t/1000)*1000).*sin(pi*i/400))+50) = 1;
        elseif mod(floor(t/500),2) == 1 %odd so contract
            Boundaries(i,floor(expansion/pniter*(ceil(t/1000)*1000-t).*sin(pi*i/400))+150) = 1; %again t val between 0 and 1000 but inversed
            Boundaries(i,ceil(expansion/pniter*(ceil(t/1000)*1000-t).*sin(pi*i/400))+150) = 1;
            Boundaries(i,floor(-expansion*(ceil(t/1000)*1000-t)/pniter.*sin(pi*i/400))+50) = 1;
            Boundaries(i,ceil(-expansion/pniter*(ceil(t/1000)*1000-t).*sin(pi*i/400))+50) = 1;
        end
    end

    
    x = 0:400;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%  POST SIM CODING / PLOT THIS BUGGER %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if mod(t, 250) == 0 || t == 1
        ConcMax = max(C,[],"all");
        fprintf('Iteration: %d, Time: %f \n', t, toc);
        figure
        tiledlayout(2,1)
        nexttile;
        hold on
        imagesc(ux.');
        colorbar;
        colormap("hot");
        clim([0 u_max]);
        xlim([0,nx]);
        ylim([0,my]);
        axis('equal');
        title('Poiseuille Flow Model')
        if mod(floor(t/500),2) == 0 %even so Expand!
            plot(x,expansion/pniter*(t-floor(t/1000)*1000).*sin(pi*x/400)+150,('b'));
            plot(x,-expansion*(t-floor(t/1000)*1000)/pniter.*sin(pi*x/400)+50,('b'));
        elseif mod(floor(t/500),2) == 1 %odd so contract
            plot(x,expansion/pniter*(ceil(t/1000)*1000-t).*sin(pi*x/400)+150,('b'));
            plot(x,-expansion*(ceil(t/1000)*1000-t)/pniter.*sin(pi*x/400)+50,('b'));
        end
        hold off

        nexttile;
        hold on
        imagesc(C.');
        colorbar;colormap("hot");
        axis('equal');
        clim([0 ConcMax]);
        xlim([0,nx]);
        ylim([0,my]);
        %title('Drug Concentration Diffusion Model')

        if mod(floor(t/500),2) == 0 %even so Expand!
            plot(x,expansion/pniter*(t-floor(t/1000)*1000).*sin(pi*x/400)+150,('b'));
            plot(x,-expansion*(t-floor(t/1000)*1000)/pniter.*sin(pi*x/400)+50,('b'));
        elseif mod(floor(t/500),2) == 1 %odd so contract
            plot(x,expansion/pniter*(ceil(t/1000)*1000-t).*sin(pi*x/400)+150,('b'));
            plot(x,-expansion*(ceil(t/1000)*1000-t)/pniter.*sin(pi*x/400)+50,('b'));
        end
        hold off

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
                    if Boundaries(a,b) == 1
                        % g(a,b,k) = 0; %dirichlet
                        g(a,b,opp(k)) = g(a,b,k); %dirichlet
                        % g(a,b,opp(k)) = gcol(a,b,k)-2*w(k)*C(a,b)*(uw(1)*cx(k)+uw(2)*cy(k))/cssq; %Moving Bounce Back
                        % if b > 100
                        %     g(a,b,k) = g(a,b-1,k);
                        % elseif b < 100
                        %     g(a,b,k) = g(a,g(a,b,k)b+1,k);
                        % end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Processing Code %
% Find maximum values
Max = max(C,[],"all");
disp(Max)

% Max_Circle = conc_level;
% Total_Circle = conc_total;
% Max_Diamond = conc_level;
% Total_Diamond = conc_total;
% Max_Oval_Peristalsis = conc_level;
% Total_Oval_Peristalsis = conc_total;
% Max_Triangle = conc_level;
% Total_Triangle = conc_total;

% Max_Star = conc_level;
% Total_Star = conc_total;

% Max_Tear = conc_level;
% Total_Tear = conc_total;

