%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Get extended bond types by considering NN bonds

function bndTypExt = extendedBondTypes( bonds_NN, CN )

Nbnd = length( bonds_NN );

bndTypExt = zeros( Nbnd, 1 );
for j = 1 : Nbnd
    bndNN = bonds_NN{j};
    Nbnd_NN = size( bndNN, 1 );
    at1 = bndNN(1,1);

    % I. C2-C2 bond:
    if Nbnd_NN == 2
        nb1 = bndNN(1,2);
        at2 = bndNN(2,1);
        nb2 = bndNN(2,2);
        %[ at1, at2, nb1, nb2 ]
        assert( CN(at1) == 2 && CN(at2) == 2 )

        if CN(nb1) == 2 && CN(nb2) == 2
            bndTypExt(j) = 1;
            % [ at1, at2, nb1, nb2 ]
            % pause
        elseif CN(nb1) == 2 && CN(nb2) == 3 || CN(nb1) == 3 && CN(nb2) == 2
            bndTypExt(j) = 2;
        elseif CN(nb1) == 3 && CN(nb2) == 3
            bndTypExt(j) = 3;
        else
            error( 'Impossible for C2-C2 bond!!!' )
        end
        continue
    end

    % II. C2-C3 bond:
    if Nbnd_NN == 3
        nb1 = bndNN(1,2);
        at2 = bndNN(2,1);
        nb2a = bndNN(2,2);
        nb2b = bndNN(3,2);
        %[ at1, at2, nb1, nb2a, nb2b ]
        assert( CN(at1) == 2 && CN(at2) == 3 )
        assert( CN(nb2a) == 3 )

        if CN(nb1) == 2 && CN(nb2b) == 2
            bndTypExt(j) = 4;
        elseif CN(nb1) == 2 && CN(nb2b) == 3
            bndTypExt(j) = 5;
        elseif CN(nb1) == 3 && CN(nb2b) == 2
            bndTypExt(j) = 6;
        elseif CN(nb1) == 3 && CN(nb2b) == 3
            bndTypExt(j) = 7;
        else
            error( 'Impossible for C2-C3 bond!!!' )
        end
        continue
    end

    % III. C3-C3 bond:
    if Nbnd_NN == 4
        nb1a = bndNN(1,2);
        nb1b = bndNN(2,2);
        at2 = bndNN(3,1);
        nb2a = bndNN(3,2);
        nb2b = bndNN(4,2);
        %[ at1, at2, nb1a, nb1b, nb2a, nb2b ]
        assert( CN(at1) == 3 && CN(at2) == 3 )

        if CN(nb1a) == 2 && CN(nb1b) == 2 && CN(nb2a) == 2 && CN(nb2b) == 2
            bndTypExt(j) = 8;
            % [ at1, at2, nb1a, nb1b, nb2a, nb2b ]
            % pause
        elseif CN(nb1a) + CN(nb1b) + CN(nb2a) + CN(nb2b) == 9
            bndTypExt(j) = 9;
        elseif ( CN(nb1a) == 2 && CN(nb1b) == 2 && CN(nb2a) == 3 && CN(nb2b) == 3 ) ...
                || ( CN(nb1a) == 3 && CN(nb1b) == 3 && CN(nb2a) == 2 && CN(nb2b) == 2 )
            bndTypExt(j) = 10;
        elseif ( CN(nb1a) == 3 && CN(nb1b) == 2 && CN(nb2a) == 3 && CN(nb2b) == 2 ) ...
                || ( CN(nb1a) == 2 && CN(nb1b) == 3 && CN(nb2a) == 2 && CN(nb2b) == 3 )
            bndTypExt(j) = 11;
        elseif ( CN(nb1a) == 2 && CN(nb1b) == 3 && CN(nb2a) == 3 && CN(nb2b) == 2 ) ...
                || ( CN(nb1a) == 3 && CN(nb1b) == 2 && CN(nb2a) == 2 && CN(nb2b) == 3 )
            bndTypExt(j) = 12;
        elseif CN(nb1a) + CN(nb1b) + CN(nb2a) + CN(nb2b) == 11
            bndTypExt(j) = 13;
        elseif CN(nb1a) == 3 && CN(nb1b) == 3 && CN(nb2a) == 3 && CN(nb2b) == 3
            bndTypExt(j) = 14;
        else
            error( 'Impossible for C3-C3 bond!!!' )
        end
        continue
    end
end

assert( all( bndTypExt > 0 ) )

end