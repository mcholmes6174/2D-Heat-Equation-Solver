function [xgrid, ygrid, tgrid, u, E] = FTCS_IBVP_2DHEAT(...
     Lx, Ly, Nx, Ny, t0, dt, Nt, u0, gL, gR, gB, gT, s, ...
     solve_until, min_tol, max_tol, max_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the FTCS method to solve the 2D heat equation
% with source term:
%
%     u_t = u_xx + u_yy + s(u)   on   (0,L_x) x (0,L_y) for t > t0
%
% with initial and boundary conditions (Neumann @ x = 0, L_x and Dirichlet
% @ y = 0, L_y) given by
%
%               u(x,y,t0) = u_0(x,y)
% !dudx! -->  u_x(0,y,t)  = g_L(y,t)   and   u_x(L_x,y,t) = g_R(y,t)
%               u(x,0,t)  = g_B(x,t)   and     u(x,L_y,t) = g_T(x,t)
% INPUTS:
%   SCALARS
%   Lx -> specifies the x-interval (0,L_x)
%   Ly -> specifies the y-interval (0,L_y)
%   Nx -> specifies the # of subintervals to use in x-discretization
%   Ny -> specifies the # of subintervals to use in y-discretization
%   t0 -> specifies initial time
%   dt -> specifies size of time step to use
%   Nt -> specifies # of updates to perform (so final time tf = t0 + dt*Nt)
%
%   MORE SCALARS (ONLY NECESSARY IF solve_until == true)
%   min_tol  -> simulation will terminate if E_n < min_tol
%   max_tol  -> simulation will terminate if E_n > max_tol
%   max_time -> simulation will terminate if t_n > max_time
%
%   2D MATRIX
%   u0(:,:) -> contains the discretized initial condition u(x,y,t0)
%
%   FUNCTIONS
%   gL(y,t) ->  left  boundary condition (Neumann)
%   gR(y,t) ->  right boundary condition (Neumann)
%   gB(x,t) -> bottom boundary condition (Dirichlet)
%   gT(x,t) ->   top  boundary condition (Dirichlet)
%    s( u ) -> source term in PDE; a function of u
%
%   LOGICAL
%   solve_until -> if TRUE, then program enters a while loop, running until
%                   certain conditions are met. Otherwise, program runs a
%                   for loop until the final time specified by the user.
% OUTPUTS:
%   1D MATRICES
%   xgrid(:) -> discretized x-spatial interval
%   ygrid(:) -> discretized y-spatial interval
%   tgrid(:) -> discretized  temporal interval
%       E(:) -> E(n) will contain the infinity norm of u_t at time
%               t = t_{n-1} computed using a 1st-order approximation
%   2D MATRIX
%   u(:,:) -> the discretized solution at final the time tf = t0 + dt*Nt
%
%       !!!! CAUTION !!!! COLUMNS of u correspond to the x-values
%                      while ROWS of u correspond to the y-values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the numerical grid
% xgrid needs to accommodate the Neumann BCs at x = 0 and x = L_x, but 
% ygrid and tgrid are constructed as usual
dx = Lx / (Nx - 1); dy = Ly / Ny;
xgrid = -0.5*dx : dx : Lx + 0.5*dx; ygrid = 0 : dy : Ly;
if solve_until == false
    tgrid = t0 : dt : t0 + Nt*dt;
    E = zeros(length(tgrid),1);
end

%% Overwrite Nx and Ny so they accurately reflect the length of the grids
Nx = length(xgrid); Ny = length(ygrid);

%% Set up constants, solution matrix u and convergece array E, and check
%  for stability
r1 = dt/(dx)^2; r2 = dt/(dy)^2; tn = t0;
u = zeros(size(u0));
if r1 + r2 > 0.5
    disp('THE METHOD IS NOT STABLE')
    return
end

%% Run the FTCS method using either a for or while loop specified by the
%  value of solve_until. 
% * We must remember that 1st index (rows) == y and 2nd index (cols) == x *
if solve_until == false
    for n = 1:Nt % want to step forward in time Nt-many times
        % get time at which to compute u(x,y)
        tn = tn + dt;    

        % update interior of the domain
        u(2:Ny-1,2:Nx-1) = u0(2:Ny-1,2:Nx-1) + ...
    r1* ( u0(2:Ny-1,3:Nx) - 2*u0(2:Ny-1,2:Nx-1) + u0(2:Ny-1,1:Nx-2) ) + ...
    r2* ( u0(3:Ny,2:Nx-1) - 2*u0(2:Ny-1,2:Nx-1) + u0(1:Ny-2,2:Nx-1) ) + ...
    dt*s( u0(2:Ny-1,2:Nx-1) );

        % update boundary of the domain
        u(:,1)   = u(:,2)     - dx*gL(ygrid,tn)'; % left
        u(:,end) = u(:,end-1) + dx*gR(ygrid,tn)'; % right
        u(1,:)   = gB(xgrid,tn); % bottom
        u(end,:) = gT(xgrid,tn); % top

        % compute E(n)
        E(n) = max(abs(u-u0)./dt, [], 'all');

        % update u0
        u0 = u;

    end
    
elseif solve_until == true
    n = 0; % index to use for E
    tgrid(1) = t0;
    while tn <= max_time
        % get time at which to compute u(x,y)
        n = n + 1; tn = tn + dt; tgrid(n) = tn;

        % update interior of the domain
        u(2:Ny-1,2:Nx-1) = u0(2:Ny-1,2:Nx-1) + ...
    r1* ( u0(2:Ny-1,3:Nx) - 2*u0(2:Ny-1,2:Nx-1) + u0(2:Ny-1,1:Nx-2) ) + ...
    r2* ( u0(3:Ny,2:Nx-1) - 2*u0(2:Ny-1,2:Nx-1) + u0(1:Ny-2,2:Nx-1) ) + ...
    dt*s( u0(2:Ny-1,2:Nx-1) );

        % update boundary of the domain
        u(:,1)   = u(:,2)     - dx*gL(ygrid,tn)'; % left
        u(:,end) = u(:,end-1) + dx*gR(ygrid,tn)'; % right
        u(1,:)   = gB(xgrid,tn); % bottom
        u(end,:) = gT(xgrid,tn); % top

        % compute E(n)
        E(n) = max(abs(u-u0)./dt, [], 'all');

        % update u0
        u0 = u;
        
        % exit if requested
        if E(n) < min_tol
            disp('CONVERGENCE REACHED')
            return
        elseif E(n) > max_tol
            disp('SOLUTION IS DIVERGING')
            return
        elseif tn > max_time
            disp('TIME LIMIT REACHED')
        end

    end
end
end