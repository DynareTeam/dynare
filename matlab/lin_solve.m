function x = lin_solve( A, b )
    if norm( b ) < sqrt( eps ) % then x = 0 is a solution
        x = 0;
        return
    end
    
    WarningState = warning( 'off', 'all' );
    
    try % wrap in a try-catch block to ensure that the WarningState is always reset
    
        x = A\b;
        x( ~isfinite( x ) ) = 0;
        relres = norm( b - A * x ) / norm( b );
        if relres > 1e-6
            fprintf( 'WARNING : Failed to find a solution to the linear system.\n' );
        end
        
    catch Error
        
        warning( WarningState );
        error( Error ); %rethrow
        
    end
    
    warning( WarningState );
    
end
