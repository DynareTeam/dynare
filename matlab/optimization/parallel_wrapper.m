function RV = parallel_wrapper( objective_function, XV, varargin )
    n = size( XV, 2 );
    RV = zeros( 1, n );
          
    parfor i = 1 : n
        R = [];
        WarningState = warning( 'off', 'all' );
        try
            R = objective_function( XV( :, i ), varargin{:} ); %#ok<PFBNS>
        catch Error
            disp( Error );
            for j = 1 : length( Error.stack )
                disp( Error.stack( j ) );
            end
        end
        warning( WarningState );
        if isempty( R ) || ~isfinite( R( 1 ) ) || ( imag( R( 1 ) ) ~= 0 )
            R = NaN;
        end
        RV( i ) = real( R( 1 ) );
    end
end
