function [ x, flag, relres ] = lin_solve_robust( A, b )
    if norm( b ) < sqrt( eps ) % then x = 0 is a solution
        x = 0;
        flag = 0;
        relres = 0;
        return
    end
    
    WarningState = warning( 'off', 'all' );
    
    try % wrap in a try-catch block to ensure that the WarningState is always reset
    
        x = A\b;
        x( ~isfinite( x ) ) = 0;
        [ x, flag, relres ] = bicgstab( A, b, [], [], [], [], x ); % returns immediately if x is a solution
        if flag == 0
            warning( WarningState );
            return
        end
        
        disp( relres );

        fprintf( 'Initial bicgstab failed, trying alternative start point.\n' );
        old_x = x;
        old_relres = relres;
        [ x, flag, relres ] = bicgstab( A, b );
        if flag == 0
            warning( WarningState );
            return
        end

        fprintf( 'Alternative start point also failed with bicgstab, trying gmres.\n' );
        if old_relres < relres
            x = old_x;
        end
        [ x, flag, relres ] = gmres( A, b, [], [], [], [], [], x );
        if flag == 0
            warning( WarningState );
            return
        end

        fprintf( 'Initial gmres failed, trying alternative start point.\n' );
        old_x = x;
        old_relres = relres;
        [ x, flag, relres ] = gmres( A, b );
        if flag == 0
            warning( WarningState );
            return
        end

        fprintf( 'Alternative start point also failed with gmres, using the (SLOW) Moore-Penrose Pseudo-Inverse.\n' );
        if old_relres < relres
            x = old_x;
            relres = old_relres;
        end
        old_x = x;
        old_relres = relres;
        x = pinv( full( A ) ) * b;
        relres = norm( b - A * x ) / norm( b );
        if old_relres < relres
            x = old_x;
            relres = old_relres;
        end
        flag = relres > 1e-6;
        if flag ~= 0
            fprintf( 'WARNING : Failed to find a solution to the linear system\n' );
        end
        
    catch Error
        
        warning( WarningState );
        error( Error ); %rethrow
        
    end
    
    warning( WarningState );

end
