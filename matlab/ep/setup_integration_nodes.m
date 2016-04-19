function [nodes,weights,nnodes] = setup_integration_nodes(EpOptions,pfm)
    if EpOptions.stochastic.order
        % Compute weights and nodes for the stochastic version of the extended path.
        switch EpOptions.IntegrationAlgorithm
          case 'Tensor-Gaussian-Quadrature'
            % Get the nodes and weights from a univariate Gauss-Hermite quadrature.
            [nodes0,weights0] = gauss_hermite_weights_and_nodes(EpOptions.stochastic.quadrature.nodes);
            % Replicate the univariate nodes for each innovation and dates, and, if needed, correlate them. 
            nodes0 = repmat(nodes0,1,pfm.number_of_shocks*pfm.stochastic_order)*kron(eye(pfm.stochastic_order),pfm.Omega);
            % Put the nodes and weights in cells
            for i=1:pfm.number_of_shocks
                rr(i) = {nodes0(:,i)};
                ww(i) = {weights0};
            end
            % Build the tensorial grid
            nodes = cartesian_product_of_sets(rr{:});
            weights = prod(cartesian_product_of_sets(ww{:}),2);
            nnodes = length(weights);
          case 'Stroud-Cubature-3'
            [nodes,weights] = cubature_with_gaussian_weight(pfm.number_of_shocks*pfm.stochastic_order,3,'Stroud')
            nodes = kron(eye(pfm.stochastic_order),transpose(pfm.Omega))*nodes;
            weights = weights;
            nnodes = length(weights);
          case 'Stroud-Cubature-5'
            [nodes,weights] = cubature_with_gaussian_weight(pfm.number_of_shocks*pfm.stochastic_order,5,'Stroud')
            nodes = kron(eye(pfm.stochastic_order),transpose(pfm.Omega))*nodes;
            weights = weights;
            nnodes = length(weights);
          case 'Unscented'
            p = pfm.number_of_shocks;
            k = EpOptions.ut.k;
            C = sqrt(pfm.number_of_shocks + k)*pfm.Omega';
            nodes = [zeros(1,p); -C; C];
            weights = [k/(p+k); (1/(2*(p+k)))*ones(2*p,1)];
            nnodes = 2*p+1;
          otherwise
            error('Stochastic extended path:: Unknown integration algorithm!')
        end
    end
