%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Revised HMO with beta depending on bond length

function [ hmoSol, bndLen ] = hmo_rev_mob( coord, POW, R0, q, twist_bonds )


if( nargin == 1 )
    POW = 0.7;
    %R0 = 1.4; % Angstrom
    R0 = mean([1.39662275 1.39662275 1.39662275 1.39662275 1.39662200 1.39662200]);
    q = 0;
    twist_bonds = [];
end
if( nargin == 2 )
    %R0 = 1.4; % Angstrom
    R0 = mean([1.39662275 1.39662275 1.39662275 1.39662275 1.39662200 1.39662200]);
    q = 0;
    twist_bonds = [];
end
if( nargin == 3 )
    q = 0;
    twist_bonds = [];
end
if( nargin == 4 )
    twist_bonds = [];
end



lm = linkage( coord );
N = length( lm );
A = zeros( N, N );
bndLen = [];
for j = 1 : N
    for k = j+1 : N
        if lm(j,k) == 0
            continue
        end
        R = norm( coord(j,:) - coord(k,:) );
        bndLen = [ bndLen; R ];

        %A(j,k) = exp( -Alpha*(R-R0) );
        A(j,k) = (R/R0)^(-POW);

        A(k,j) = A(j,k);
    end
end

% Inverse phases:
for j = 1 : size( twist_bonds, 1 )
    A( twist_bonds(j,1), twist_bonds(j,2) ) = -A( twist_bonds(j,1), twist_bonds(j,2) );
    A( twist_bonds(j,2), twist_bonds(j,1) ) = -A( twist_bonds(j,2), twist_bonds(j,1) );
end


[ C, E ] = eig( A );

E = -diag(E); % in units of -beta
% sort E from negative to positive values
[ E, IX ] = sort( E );
% sort V with the order of E
C = C( :, IX );

% Applying Hund's rule and Aufbau principle to fill the levels with
% alpha and beta electrons:
% (1) Degeneracies of of levels:
j = 1;
for k = 1 : length(E)
    if k > 1 && abs( E1(j-1) - E(k)) < eps*20
        Deg(j-1) = Deg(j-1) + 1;
    else
        E1(j) = E(k);
        Deg(j) = 1;
        j = j + 1;
    end
end

% (2) Fill the electrons:
Ne = N - q; % total number of electrons
Occ_a = zeros( N, 1 ); % alpha MO occupancies:
Occ_b = zeros( N, 1 ); % beta MO occupancies:
counter = 0;
ia = 1;
ib = 1;
for k = 1 : length( Deg )
    deg_k = Deg(k);
    % Fill with alpha electrons:    
    for j = 1 : deg_k
        if counter == Ne % No eletron left
            break;
        end
        counter = counter + 1;
        Occ_a( ia ) = 1;
        ia = ia + 1;
    end
    % Fill with beta electrons:
    for j = 1 : deg_k
        if counter == Ne % No eletron left
            break;
        end
        counter = counter + 1;
        Occ_b( ib ) = 1;
        ib = ib + 1;
    end
end

Ne_a = sum( Occ_a ); % number of alpha electrons
Ne_b = sum( Occ_b ); % number of beta electrons

Occ = Occ_a + Occ_b;
Occ_spin = Occ_a - Occ_b;

Mult = sum( Occ_a - Occ_b ) + 1;

% Total energy
Etot = sum( Occ.*E );

% HOMO-LUMO gap
if( Ne_a ~= Ne_b ) % open-shell
    Gap = 0;
else
    Gap = E( Ne_a + 1 ) - E( Ne_a );
end

% Atomic charges:
Chg = ones(N,1); % initially each atom has one pi-electron
for iat = 1 : N
    Chg(iat) = Chg(iat) - sum( Occ'.*(C(iat,:).^2) );
end

% Atomic spin charges:
Spin = zeros(N,1);
for iat = 1 : N
    Spin(iat) = sum( Occ_spin'.*(C(iat,:).^2) );
end

% Bond orders
for i1 = 1 : N
    for i2 = 1 : N
        ifTwistBond = false;
        for k = 1 : size( twist_bonds, 1 )
            if isequal([i1 i2], twist_bonds(k,:)) || ...
                    isequal([i2 i1], twist_bonds(k,:))
                ifTwistBond = true;
                break
            end
        end
        if ifTwistBond
            BOs( i1, i2 ) = sum( -Occ'.*( C(i1,:).*C(i2,:) ) );
        else
            BOs( i1, i2 ) = sum( Occ'.*( C(i1,:).*C(i2,:) ) );
        end
    end
end

% Bond orders only for connected bonds:
nBond = 0;
for i1 = 1 : N
    for i2 = i1 + 1 : N
        if( lm(i1,i2) )
            nBond = nBond + 1;
            BO( nBond, 1 ) = i1;
            BO( nBond, 2 ) = i2;
            BO( nBond, 3 ) = BOs( i1, i2 );
        end
    end
end

% Free valence:
F = zeros(N,1); % initialize
for ib = 1 : size( BO, 1 )
    i1 = BO( ib, 1 );
    i2 = BO( ib, 2 );
    F(i1) = F(i1) + BO( ib, 3 );
    F(i2) = F(i2) + BO( ib, 3 );
end
F = (3+sqrt(3)) - (F+3); % there are 3 sigma-bonds already


% Density matrix:
Da = C(:,1:Ne_a)*C(:,1:Ne_a)';
Db = C(:,1:Ne_b)*C(:,1:Ne_b)';
D = Da + Db;

% Jahn-Teller distortion:
JT = false;
if abs( E(Ne_a) - E(Ne_a+1) ) < eps*20
    JT = true;
end
if abs( E(Ne_b) - E(Ne_b+1) ) < eps*20
    JT = true;
end

hmoSol.E = E;
hmoSol.C = C;
hmoSol.N = N;
hmoSol.Ne = Ne;
hmoSol.Ne_a = Ne_a;
hmoSol.Ne_b = Ne_b;
hmoSol.Occ_a = Occ_a;
hmoSol.Occ_b = Occ_b;
hmoSol.Occ = Occ;
hmoSol.Etot = Etot;
hmoSol.Gap = Gap;
hmoSol.Chg = Chg;
hmoSol.Spin = Spin;
hmoSol.Mult = Mult;
hmoSol.BOs = BOs;
hmoSol.BO = BO;
hmoSol.F = F;
hmoSol.D = D;
hmoSol.Da = Da;
hmoSol.Db = Db;
hmoSol.JT = JT;
hmoSol.xyz = coord;
hmoSol.twist_bonds = twist_bonds;

end
