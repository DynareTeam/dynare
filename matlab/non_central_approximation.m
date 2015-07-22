function [ info, M, options, oo ] = non_central_approximation( M, options, oo )
    options.original_order = options.order;
    deflect_ = compute_deflected_linear_approximation( M, options, oo, options.non_central_approximation );
    options.order = 1;
    oo.dr.ys = deflect_.y;
    oo.steady_state = deflect_.y;
    oo.dr.ghx = deflect_.y_x;
    oo.dr.ghu = deflect_.y_u;
    
    info = 0;
    if any( isnan( deflect_.y ) )
        info = 12345;
    end
    if any( any( isnan( deflect_.y_x ) ) )
        info = 12345;
    end
    if any( any( isnan( deflect_.y_u ) ) )
        info = 12345;
    end
end