function J = Jacobian(n,P,C,M,Da,Pe,dt,dZ)

%  parameters
%n                            %number of nodes
%a       % m        % Characteristic pore/channel radius
%k       % m/Pa.s   % Hydraulic permeability  of the membrane
%T       % K        % Absolute temperature of the system
%Rg      % J/mol.K  % Universal gas constant 
%l       % m        % Total length of the phloem
%dz=l/n                        %spatial step size
%c0      % mol/L    % Inlet boundary concentration
%psil     % Pa              % xylem water potential
%mu       % Pa.s     % Dynamic viscosity of the fluid
% D       % m^2/s    % Molecular diffusion coefficient of the species
% tr is recation time scale



% location vector
i = 2:n-1;


% allocate J 
J = zeros(2*n,2*n);


%p equation terms
J(sub2ind(size(J),2.*i-1,2.*(i-1)-1))=-1/(M*dZ^2); % p-1 variable
J(sub2ind(size(J),2.*i-1,2.*i-1))=2/(M*dZ^2) + 1; % p variab;e
J(sub2ind(size(J),2.*i-1,2.*i))=-1; % c variable
J(sub2ind(size(J),2.*i-1,2.*(i+1)-1))=-1/(M*dZ^2); % p+1 variable


% c equation terms
J(sub2ind(size(J),2.*i,2.*(i-1)-1))=(1/dZ^2).* C(1:n-2); % p-1 variable
J(sub2ind(size(J),2.*i,2.*(i-1)))= 1/dZ^2.*( P(2:n-1)-P(1:n-2)) -1/(Pe*dZ^2) ; % c-1 variable
J(sub2ind(size(J),2.*i,2.*i-1))=  (1/dZ^2).*(C(2:n-1)+C(1:n-2)) ; % p variable
J(sub2ind(size(J),2.*i,2.*i))= 1/dt - 1/dZ^2.*(P(3:n) - P(2:n-1) )+2/(Pe*dZ^2) + Da; % c variable
J(sub2ind(size(J),2.*i,2.*(i+1)-1))= (-1/dZ^2).*C(2:n-1); % p+1 variable
J(sub2ind(size(J),2.*i,2.*(i+1)))=-1/(Pe*dZ^2); % c+1 variable

%BC
%at top
% p equation
J(1,1)=1; 
J(1,2)=-1;
% c equation
J(2,2)=1;

% at bottom
% p equation
J(2*n-1,2*n-1)=1;
J(2*n-1,2*n)= -1;
% c equation
J(2*n,2*n)=1;
end