function [PDE, F, domain, Method0, Method1, DiscMethod] = ClearOrder(Conditions, Function, Domain, Coefficients, Number_of_Points)


a = Domain(1); b = Domain(2);
v = Coefficients(1); betta = Coefficients(2); gamma = Coefficients(3);
n = Number_of_Points; % number of points
h = (b - a)/(n+1); % step
L = (1/(h^2))*toeplitz([2 -1 zeros(1, n-2)]);


x = a+h:h:b-h;
domain = x;
f = Function(x)';


%% Kinds of Conditions
if Conditions(1) == 'D'
    U0 = Conditions(2);
elseif Conditions(1) == 'N'
    UPrime0 = Conditions(2);
end

if Conditions(3) == 'D'
    U1 = Conditions(4);
elseif Conditions(3) == 'N'
    UPrime1 = Conditions(4);
end

%% Centered Model for first Differential    U(i+1) - U(i-1)

if Conditions(5) == 'C'
Centered_D = (1/(2*h))*toeplitz([0 -1 zeros(1,n-2)], [0; 1; zeros(n-2, 1)]);
end

%% Upwind Model for first Differential U(i) - U(i-1)

if Conditions(5) == 'U'
Upwind_D = (1/h)*toeplitz([1 -1 zeros(1, n-2)], [1;zeros(n-1, 1)]);
end

%% Boundry Condition for U(a) = U0 // Dirichlet BC  // Centered          1

if Conditions(5) == 'C' && Conditions(1) == 'D'
DirichletBC_centeredU0 = [v*U0/(h^2) + betta*U0/(2*h); zeros(n-1, 1)];
end

%% Boundry Condition for U(b) = U1 // Dirichlet BC  // Centered          2

if Conditions(5) == 'C' && Conditions(3) == 'D'
DirichletBC_centeredU1 = [zeros(n-1, 1); v*U1/(h^2) - betta*U1/(2*h) ];
end

%% Boundry Condition for U(a) = U0 // Dirichlet BC  // UPwind            3

if Conditions(5) == 'U' && Conditions(1) == 'D'
DirichletBC_upwindU0 = [v*U0/(h^2) + betta*U0/(h); zeros(n-1, 1)];
end

%% Boundry Condition for U(b)= U1 // Dirichlet BC  // UPwind             4

if Conditions(5) == 'U' && Conditions(3) == 'D'
DirichletBC_upwindU1 = [zeros(n-1, 1); v*U1/(h^2) - betta*U1/(h)];
end

%% Boundry Condition for U'(a) = UPrime0 // Neumann BC // Centered       5    U(i+1) - U(i) = h*UPrime0   for i = 0

if Conditions(5) == 'C' && Conditions(1) == 'N'
L_Numann0 = L - diag([(1/h^2) zeros(1, n-1)]);
Centered_D_Numann0 = Centered_D + diag([-1/(2*h) zeros(1,n-1)]);
NumannBC_Centered_difU0 = [UPrime0*v/h - betta*UPrime0/2 ; zeros(n-1, 1)];
end

%% Boundry Condition for U'(b) = UPrime1 // Neumann BC // Centered       6    U(i+1) - U(i) = h*UPrime1   for i = n

if Conditions(5) == 'C' && Conditions(3) == 'N'
L_Numann1 = L - diag([zeros(1, n-1) (1/h^2)]);
Centered_D_Numann1 = Centered_D + diag([zeros(1,n-1) 1/(2*h)]);
NumannBC_Centered_difU1 = [zeros(n-1, 1); UPrime1*v/h - betta*UPrime1/2 ];
end

%% Boundry Condition for U'(a) = UPrime0 // Neumann BC // Upwind         7   U(i+1) - U(i) = h*UPrime0   for i = 0

if Conditions(5) == 'U' && Conditions(1) == 'N'
L_Numann0 = L - diag([(1/h^2) zeros(1, n-1)]);
Upwind_D_Numman0 = Upwind_D - diag([1/h zeros(1, n-1)]);
NumannBC_Upwind_difU0 = [v*UPrime0/h - UPrime0 ; zeros(n-1, 1)];
end

%% Boundry Condition for U'(b) = UPrime1 // Neumann BC // Upwind         8    U(i+1) - U(i) = h*UPrime0   for i = n

if Conditions(5) == 'U' && Conditions(3) == 'N'
L_Numann1 = L - diag([zeros(1, n-1) (1/h^2)]);
Upwind_D_Numman1 = Upwind_D ; %% doesn't change
NumannBC_Upwind_difU1 = [zeros(n-1, 1); v*UPrime1/h];
end

%% just Copy Function to not change our main Function

F = f;

%% Classification for First Boundry Condition
if Conditions(1) == 'D'
    L_op = L;
    if Conditions(5) == 'U'
        F = F + DirichletBC_upwindU0;
        D_op = Upwind_D;
    elseif Conditions(5) == 'C'
        F = F + DirichletBC_centeredU0;
        D_op = Centered_D;
    end
elseif Conditions(1) == 'N'
    L_op = L_Numann0;
    if Conditions(5) == 'U'
        D_op = Upwind_D_Numman0;
        F = F + NumannBC_Upwind_difU0;
    elseif Conditions(5) == 'C'
        F = F + NumannBC_Centered_difU0;
        D_op = Centered_D_Numann0;
    end
end

%% Classification for Second Boundry Condition

if Conditions(3) == 'D'

    L_op = L;
    if Conditions(5) == 'U'
        F = F + DirichletBC_upwindU1;
        D_op = Upwind_D;
    elseif Conditions(5) == 'C'
        F = F + DirichletBC_centeredU1;
        D_op = Centered_D; 
    end
elseif Conditions(3) == 'N'
    L_op = L_Numann1;
    if Conditions(5) == 'U'
        F = F + NumannBC_Upwind_difU1;
        D_op = Upwind_D_Numman1;
    elseif Conditions(5) == 'C'
        F = F + NumannBC_Centered_difU1;
        D_op = Centered_D_Numann1;
    end
end

%% Create PDE Eq. and Some notations

if Conditions(1) == 'D'
    Method0 = 'Dirichlet';
else
    Method0 = 'Nemman';
end

if Conditions(3) == 'D'
    Method1 = 'Dirichlet';
else
    Method1 = 'Nemman';
end

if Conditions(5) == 'U'
    DiscMethod = 'Upwind';
elseif Conditions(5) == 'C'
    DiscMethod = 'Centered';
end

PDE = v*L_op + betta*D_op + gamma*eye(n);





