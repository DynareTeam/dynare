function [ info, M, options, oo ] = gaussian_approximation( M, options, oo )
    % derived from CacheConditionalCovariancesAndAugmentedStateTransitionMatrices.m
    info = 0;
    [numeric_version] = return_dynare_version(dynare_version);
    if numeric_version >= 4.4 
        nstatic = M.nstatic;
        nspred = M.nspred; % note M.nspred = M.npred+M.nboth;
    else
        nstatic = oo.dr.nstatic;
        nspred = oo.dr.npred;
    end
    SelectState = ( nstatic + 1 ):( nstatic + nspred );
    nEndo = M.endo_nbr;
    
    % pre-calculations for order=1 terms
    A1 = sparse( nEndo, nEndo );
    A1( :, SelectState ) = oo.dr.ghx;
    A1( abs(A1)<eps ) = 0;

    B1 = spsparse( oo.dr.ghu );

    Sigma = spsparse( M.Sigma_e );

    % pre-calculations for order=2 terms
    nState = length( SelectState );
    nExo = M.exo_nbr;
    nExo2 = nExo * nExo;

    % Idx1 =  1:nEndo;
    % Idx2 = (nEndo+1):(2*nEndo);
    % Idx3 = (2*nEndo+1):LengthZ2;

    % A2( Idx1, Idx1 ) = A1;
    [ A2i, A2j, A2s ] = spfind( A1 );

    % A2( Idx2, Idx2 ) = A1;
    A2i = [ A2i; A2i + nEndo ];
    A2j = [ A2j; A2j + nEndo ];
    A2s = [ A2s; A2s ];

    % A2( Idx2, Idx3 ) = 0.5 * oo.dr.ghxx;
    beta22 = spsparse( oo.dr.ghxx );
    [ Tmpi, Tmpj, Tmps ] = find( beta22 );
    A2i = [ A2i; Tmpi + nEndo ];
    A2j = [ A2j; Tmpj + 2*nEndo ];
    A2s = [ A2s; 0.5 * Tmps ];

    % A2( Idx3, Idx3 ) = spkron( oo.dr.ghx( SelectState, : ), oo.dr.ghx( SelectState, : ) );
    A1S = A1( SelectState, SelectState );
    A1S2 = spkron( A1S, A1S );
    [ Tmpi, Tmpj, Tmps ] = find( A1S2 );
    A2i = [ A2i; Tmpi + 2*nEndo ];
    A2j = [ A2j; Tmpj + 2*nEndo ];
    A2s = [ A2s; Tmps ];

    nState2 = nState * nState;
    LengthZ2 = 2 * nEndo + nState2;

    A2 = sparse( A2i, A2j, A2s, LengthZ2, LengthZ2 );

    B1S = B1( SelectState, : );
    B1S2 = ( spkron( B1S, B1S ) );


    K_nState_nState = commutation_sparse( nState, nState );


    % Jdx1 = 1:nExo;
    % Jdx2 = (nExo+1):(nExo + nExo2);
    % Jdx3 = (nExo + nExo2 + 1):LengthXi;

    % B2( Idx1, Jdx1 ) = oo.dr.ghu;
    [ B2i, B2j, B2s ] = find( B1 );

    % B2( Idx2, Jdx2 ) = 0.5 * oo.dr.ghuu;
    [ Tmpi, Tmpj, Tmps ] = spfind( 0.5 * oo.dr.ghuu );
    B2i = [ B2i; Tmpi + nEndo ];
    B2j = [ B2j; Tmpj + nExo ];
    B2s = [ B2s; Tmps ];

    % B2( Idx2, Jdx3 ) = oo.dr.ghxu;
    [ Tmpi, Tmpj, Tmps ] = spfind( oo.dr.ghxu );
    B2i = [ B2i; Tmpi + nEndo ];
    B2j = [ B2j; Tmpj + nExo + nExo2 ];
    B2s = [ B2s; Tmps ];

    % B2( Idx3, Jdx2 ) = spkron( oo.dr.ghu( SelectState, : ), oo.dr.ghu( SelectState, : ) );
    [ Tmpi, Tmpj, Tmps ] = find( B1S2 );
    B2i = [ B2i; Tmpi + 2*nEndo ];
    B2j = [ B2j; Tmpj + nExo ];
    B2s = [ B2s; Tmps ];

    % B2( Idx3, Jdx3 ) = ( speye( nState * nState ) + commutation_sparse( nState, nState ) ) * spkron( oo.dr.ghx( SelectState, : ), oo.dr.ghu( SelectState, : ) );
    [ Tmpi, Tmpj, Tmps ] = find( ( speye( nState2 ) + K_nState_nState ) * spkron( A1S, B1S ) );
    B2i = [ B2i; Tmpi + 2*nEndo ];
    B2j = [ B2j; Tmpj + nExo + nExo2 ];
    B2s = [ B2s; Tmps ];

    LengthXi = nExo + nExo2 + nState * nExo;

    B2 = sparse( B2i, B2j, B2s, LengthZ2, LengthXi );  

    % BCovXiB{i}( Jdx1, Jdx1 ) = Sigma;
    [ Vi, Vj, Vs ] = find( Sigma );

    % BCovXiB{i}( Jdx2, Jdx2 ) = dynareOBC_.Variance_exe;
    [ Tmpi, Tmpj, Tmps ] = find( ( speye( nExo2 ) + commutation_sparse( nExo, nExo ) ) * spkron( Sigma, Sigma ) );
    Vi = [ Vi; Tmpi + nExo ];
    Vj = [ Vj; Tmpj + nExo ];
    Vs = [ Vs; Tmps ];

    Var_z1 = SparseLyapunovSymm( A1, B1*Sigma*B1' );
    [ Tmpi, Tmpj, Tmps ] = spkron( Var_z1( SelectState, SelectState ), Sigma );
    Ui = [ Vi; Tmpi + nExo + nExo2 ];
    Uj = [ Vj; Tmpj + nExo + nExo2 ];
    Us = [ Vs; Tmps ];

    UnconditionalVarXi = sparse( Ui, Uj, Us, LengthXi, LengthXi );

    % Var_z = SparseLyapunovSymm( A2, B2*UnconditionalVarXi*B2' );
    
    [ MOi, MOj, MOs ] = spfind( 0.5 * ( oo.dr.ghuu * vec( Sigma ) ) );
    [ Tmpi, Tmpj, Tmps ] = spfind( B1S2 * vec( Sigma ) );
    MOi = [ MOi + nEndo; Tmpi + 2*nEndo ];
    MOj = [ MOj; Tmpj ];
    MOs = [ MOs; Tmps ];
    c = sparse( MOi, MOj, MOs, LengthZ2, 1 );
    Mean_z = ( speye( LengthZ2 ) - A2 ) \ c;
    Mean_y = Mean_z( 1:nEndo ) + Mean_z( (nEndo+1):(2*nEndo) );
    Mean_y = oo.dr.ys + full( Mean_y( oo.dr.inv_order_var ) );
    
    options.original_order = options.order;
    options.order = 1;
    
    oo.original_dr = oo.dr;
    A2 = full( A2 ); % sub-optimal, but easier
    new_nEndo = nEndo + size( A2, 1 );
    oo.dr.ghx = [ A2( oo.dr.inv_order_var, : ) + A2( nEndo+oo.dr.inv_order_var, : ); A2 ];
    TempB2 = full( B2 ) * reduced_rank_cholesky( full( UnconditionalVarXi ) )';
    oo.dr.ghu = [ TempB2( oo.dr.inv_order_var, : ) + TempB2( nEndo+oo.dr.inv_order_var, : ); TempB2 ];
    oo.dr.ys = [ Mean_y; zeros( new_nEndo - nEndo, 1 ) ];
    oo.steady_state = oo.dr.ys;
    
    M.original_M = M;
    M.endo_nbr = new_nEndo;
    new_nExo = size( TempB2, 2 );
    M.Sigma_e = eye( new_nExo );
    M.H = [];
    M.sigma_e_is_diagonal = 1;
    M.nstatic = nEndo;
    M.nfwrd = 0;
    M.npred = LengthZ2;
    M.nboth = 0;
    M.nsfwrd = 0;
    M.nspred = LengthZ2;
    M.ndynamic = LengthZ2;
    oo.dr.eigval = [];
    oo.dr.order_var = [ oo.dr.order_var; ((nEndo+1):new_nEndo)' ];
    oo.dr.inv_order_var = [ oo.dr.inv_order_var; ((nEndo+1):new_nEndo)' ];
    oo.dr.kstate = [];
    oo.dr.restrict_var_list = [];
    oo.dr.restrict_columns = [];
    new_endo_names = repmat( ' ', new_nEndo, size( M.endo_names, 2 ) );
    new_endo_names( 1:nEndo, : ) = M.endo_names;
    M.endo_names = new_endo_names;
    new_exo_names = repmat( ' ', new_nExo, size( M.exo_names, 2 ) );
    new_exo_names( 1:M.exo_nbr, : ) = M.exo_names;
    M.exo_names = new_exo_names;
    M.lead_lag_incidence = zeros( 3, new_nEndo );
    M.lead_lag_incidence( 1, (nEndo+1):end ) = 1:(new_nEndo-nEndo);
    M.lead_lag_incidence( 2, : ) = (new_nEndo-nEndo+1):(new_nEndo-nEndo+new_nEndo);
    oo.dr = set_state_space( oo.dr, M, options );
end