function [res,A,info] = ep_problem_2(y,x,pm)

info = 0;
res = [];
A = [];

dynamic_model = pm.dynamic_model;
ny = pm.ny;
params = pm.params;
steady_state = pm.steady_state;
order = pm.order;
nodes = pm.nodes;
nnodes = pm.nnodes;
weights = pm.weights;
h_correction = pm.h_correction;
dimension = pm.dimension;
world_nbr = pm.world_nbr;
nnzA = pm.nnzA;
periods = pm.periods;
i_rows = pm.i_rows';
i_cols = pm.i_cols;
nyp = pm.nyp;
nyf = pm.nyf;
hybrid_order = pm.hybrid_order;
i_cols_1 = pm.i_cols_1;
i_cols_j = pm.i_cols_j;
icA = pm.icA;
i_cols_T = pm.i_cols_T;
eq_index = pm.eq_index;

i_cols_p = i_cols(1:nyp);
i_cols_s = i_cols(nyp+(1:ny));
i_cols_f = i_cols(nyp+ny+(1:nyf));
i_cols_A = i_cols;
i_cols_Ap0 = i_cols_p;
i_cols_As = i_cols_s;
i_cols_Af0 = i_cols_f - ny;
i_hc = i_cols_f - 2*ny;

nzA = cell(periods,world_nbr);
res = zeros(ny,periods,world_nbr);
Y = zeros(ny*(periods+2),world_nbr);
Y(1:ny,1) = pm.y0;
Y(end-ny+1:end,:) = repmat(steady_state,1,world_nbr);
Y(pm.i_upd_y) = y;
offset_r0 = 0;
for i = 1:order+1
    i_w_p = 1;
    for j = 1:(1+(nnodes-1)*(i-1))
        innovation = x;
        if i <= order && j == 1
            % first world, integrating future shocks
            if nargout > 1
                A1 = sparse([],[],[],i*(1+(nnodes-1)*(i-1))*ny,dimension,nnzA*world_nbr);
            end
            for k=1:nnodes
                if nargout > 1
                    if i == 2
                        i_cols_Ap = i_cols_Ap0;
                    elseif i > 2
                        i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes- ...
                                                          1)*(i-2)*(i-3)/2);
                    end
                    if k == 1
                        i_cols_Af = i_cols_Af0 + ny*(i-1+(nnodes-1)*i*(i-1)/ ...
                                                     2);
                    else
                        i_cols_Af = i_cols_Af0 + ny*(i-1+(nnodes-1)*i*(i-1)/ ...
                                                     2+(i-1)*(nnodes-1) ...
                                                     + k - 1);
                    end
                end
                if i > 1
                    innovation(i+1,:) = nodes(k,:);
                end
                if k == 1
                    k1 = 1;
                else
                    k1 = (nnodes-1)*(i-1)+k;
                end
                if hybrid_order == 2 && (k > 1 || i == order) 
                    z = [Y(i_cols_p,1);
                         Y(i_cols_s,1);
                         Y(i_cols_f,k1)+h_correction(i_hc)];
                else
                    z = [Y(i_cols_p,1);
                         Y(i_cols_s,1);
                         Y(i_cols_f,k1)];
                end
                if nargout > 1
                    [d1,jacobian] = dynamic_model(z,innovation,params,steady_state,i+1);
                    if i == 1
                        % in first period we don't keep track of
                        % predetermined variables
                        i_cols_A = [i_cols_As - ny; i_cols_Af];
                        A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(eq_index,i_cols_1);  
                    else
                        i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                        A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(eq_index,i_cols_j);  
                    end
                else
                    d1 = dynamic_model(z,innovation,params,steady_state,i+1);
                end
                res(:,i,1) = res(:,i,1)+weights(k)*d1(eq_index);
            end
            if nargout > 1
                [ir,ic,v] = find(A1);
                nzA{i,j} = [ir,ic,v]';
            end
        elseif j > 1 + (nnodes-1)*(i-2)
            % new world, using previous state of world 1 and hit
            % by shocks from nodes
            if nargout > 1
                i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes-1)*(i-2)*(i-3)/2);
                i_cols_Af = i_cols_Af0 + ny*(i+(nnodes-1)*i*(i-1)/2+j-2);
            end
            k = j - (nnodes-1)*(i-2);
            innovation(i+1,:) = nodes(k,:);
            z = [Y(i_cols_p,1);
                 Y(i_cols_s,j);
                 Y(i_cols_f,j)];
            if nargout > 1
                [d1,jacobian] = dynamic_model(z,innovation,params,steady_state,i+1);
                i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                [ir,ic,v] = find(jacobian(eq_index,i_cols_j));
                nzA{i,j} = [i_rows(ir),i_cols_A(ic), v]';
            else
                d1 = dynamic_model(z,innovation,params,steady_state,i+1);
            end
            res(:,i,j) = d1(eq_index);
            if nargout > 1
                i_cols_Af = i_cols_Af + ny;
            end
        else
            % existing worlds other than 1
            if nargout > 1
                i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes-1)*(i-2)*(i-3)/2+j-1);
                i_cols_Af = i_cols_Af0 + ny*(i+(nnodes-1)*i*(i-1)/2+j-2);
            end
            z = [Y(i_cols_p,j);
                 Y(i_cols_s,j);
                 Y(i_cols_f,j)];
            if nargout > 1
                [d1,jacobian] = dynamic_model(z,innovation,params,steady_state,i+1);
                i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                [ir,ic,v] = find(jacobian(eq_index,i_cols_j));
                nzA{i,j} = [i_rows(ir),i_cols_A(ic),v]';
                i_cols_Af = i_cols_Af + ny;
            else
                d1 = dynamic_model(z,innovation,params,steady_state,i+1);
            end
            res(:,i,j) = d1(eq_index);
        end
        i_rows = i_rows + ny;
        if mod(j,nnodes) == 0
            i_w_p = i_w_p + 1;
        end
        if nargout > 1 && i > 1
            i_cols_As = i_cols_As + ny;
        end
        offset_r0 = offset_r0 + ny;
    end
    i_cols_p = i_cols_p + ny;
    i_cols_s = i_cols_s + ny;
    i_cols_f = i_cols_f + ny;
end
for j=1:world_nbr
    i_rows_y = i_cols+(order+1)*ny;
    offset_c = ny*(order+(nnodes-1)*(order-1)*order/2+j-1);
    offset_r = offset_r0+(j-1)*ny;
    for i=order+2:periods
        if nargout > 1
            [d1,jacobian] = dynamic_model(Y(i_rows_y,j),x,params, ...
                                          steady_state,i+1);
            if i < periods
                [ir,ic,v] = find(jacobian(eq_index,i_cols_j));
            else
                [ir,ic,v] = find(jacobian(eq_index,i_cols_T));
            end
            nzA{i,j} = [offset_r+ir,offset_c+icA(ic), v]';
        else
            d1 = dynamic_model(Y(i_rows_y,j),x,params, ...
                                          steady_state,i+1);
        end
        res(:,i,j) = d1(eq_index);
        i_rows_y = i_rows_y + ny;
        offset_c = offset_c + world_nbr*ny;
        offset_r = offset_r + world_nbr*ny;
    end
end
if nargout > 1
    iA = [nzA{:}]';
    A = sparse(iA(:,1),iA(:,2),iA(:,3),dimension,dimension);
end
res = res(pm.i_upd_r);