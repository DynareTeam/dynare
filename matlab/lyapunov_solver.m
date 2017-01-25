function P=lyapunov_solver(T,R,Q,DynareOptions) % --*-- Unitary tests --*--
% function P=lyapunov_solver(T,R,Q,DynareOptions)
% Solves the Lyapunov equation P-T*P*T' = R*Q*R' arising in a state-space
% system, where P is the variance of the states 
%
% Inputs
%   T               [double]    n*n matrix.
%   R               [double]    n*m matrix.
%   Q               [double]    m*m matrix.
%   DynareOptions   [structure] Dynare options
% 
% Outputs
%   P               [double]    n*n matrix.
%
% Algorithms
%   Default, if none of the other algorithms is selected:
%       Reordered Schur decomposition (Bartels-Stewart algorithm)
%   DynareOptions.lyapunov_fp == 1
%       iteration-based fixed point algorithm
%   DynareOptions.lyapunov_db == 1
%       doubling algorithm
%   DynareOptions.lyapunov_srs == 1
%       Square-root solver for discrete-time Lyapunov equations (requires Matlab System Control toolbox 
%       or Octave control package)

% Copyright (C) 2016 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if DynareOptions.lyapunov_fp == 1
    P = lyapunov_symm(T,R*Q'*R',DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold, 3, DynareOptions.debug);
elseif DynareOptions.lyapunov_db == 1
    [P, errorflag] = disclyap_fast(T,R*Q*R',DynareOptions.lyapunov_doubling_tol);
    if errorflag %use Schur-based method
        P = lyapunov_symm(T,R*Q*R',DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold, [], DynareOptions.debug);
    end
elseif DynareOptions.lyapunov_srs == 1
    % works only with Matlab System Control toolbox or Octave control package,
    if isoctave
        if ~user_has_octave_forge_package('control')
            error('lyapunov=square_root_solver not available; you must install the control package from Octave Forge')
        end
    else
        if ~user_has_matlab_license('control_toolbox')
            error('lyapunov=square_root_solver not available; you must install the control system toolbox')
        end
    end
    chol_Q = R*chol(Q,'lower');
    R_P = dlyapchol(T,chol_Q);
    P = R_P' * R_P;    
else
    P = lyapunov_symm(T,R*Q*R',DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold, [], DynareOptions.debug);
end;

%@test:1
%$ t = NaN(10,1);
%$ options_.lyapunov_complex_threshold = 1e-15;
%$ options_.qz_zero_threshold = 1e-6;
%$ options_.qz_criterium=1-options_.qz_zero_threshold;
%$ options_.lyapunov_fixed_point_tol = 1e-10;
%$ options_.lyapunov_doubling_tol = 1e-16;
%$ options_.debug=0;
%$
%$ n_small=8;
%$ m_small=10;
%$ T_small=randn(n_small,n_small);
%$ T_small=0.99*T_small/max(abs(eigs(T_small)));
%$ tmp2=randn(m_small,m_small);
%$ Q_small=tmp2*tmp2';
%$ R_small=randn(n_small,m_small);
%$ 
%$ n_large=9;
%$ m_large=11;
%$ T_large=randn(n_large,n_large);
%$ T_large=0.99*T_large/max(abs(eigs(T_large)));
%$ tmp2=randn(m_large,m_large);
%$ Q_large=tmp2*tmp2';
%$ R_large=randn(n_large,m_large);
%$ 
%$ % DynareOptions.lyapunov_fp == 1
%$ options_.lyapunov_fp = 1;
%$ try
%$    Pstar1_small = lyapunov_solver(T_small,R_small,Q_small,options_);
%$    Pstar1_large = lyapunov_solver(T_large,R_large,Q_large,options_);
%$    t(1) = 1;
%$ catch
%$    t(1) = 0;
%$ end
%$
%$ % Dynareoptions.lyapunov_db == 1
%$ options_.lyapunov_fp = 0;
%$ options_.lyapunov_db = 1;
%$ try
%$    Pstar2_small = lyapunov_solver(T_small,R_small,Q_small,options_);
%$    Pstar2_large = lyapunov_solver(T_large,R_large,Q_large,options_);
%$    t(2) = 1;
%$ catch
%$    t(2) = 0;
%$ end
%$
%$ % Dynareoptions.lyapunov_srs == 1
%$ if (isoctave && user_has_octave_forge_package('control')) || (~isoctave && user_has_matlab_license('control_toolbox'))
%$     options_.lyapunov_db = 0;
%$     options_.lyapunov_srs = 1;
%$     try
%$        Pstar3_small = lyapunov_solver(T_small,R_small,Q_small,options_);
%$        Pstar3_large = lyapunov_solver(T_large,R_large,Q_large,options_);
%$        t(3) = 1;
%$     catch
%$        t(3) = 0;
%$     end
%$ else
%$     t(3) = 1;
%$ end
%$
%$ % Standard
%$     options_.lyapunov_srs = 0;
%$ try
%$    Pstar4_small = lyapunov_solver(T_small,R_small,Q_small,options_);
%$    Pstar4_large = lyapunov_solver(T_large,R_large,Q_large,options_);
%$    t(4) = 1;
%$ catch
%$    t(4) = 0;
%$ end
%$
%$ % Test the results.
%$
%$ if max(max(abs(Pstar1_small-Pstar2_small)))>1e-8
%$    t(5) = 0;
%$ else
%$    t(5) = 1;
%$ end
%$
%$ if (isoctave && user_has_octave_forge_package('control')) || (~isoctave && user_has_matlab_license('control_toolbox'))
%$    if max(max(abs(Pstar1_small-Pstar3_small)))>1e-8
%$       t(6) = 0;
%$    else
%$       t(6) = 1;
%$    end
%$ else
%$    t(6) = 1;
%$ end
%$
%$ if max(max(abs(Pstar1_small-Pstar4_small)))>1e-8
%$    t(7) = 0;
%$ else
%$    t(7) = 1;
%$ end
%$
%$ if max(max(abs(Pstar1_large-Pstar2_large)))>1e-8
%$    t(8) = 0;
%$ else
%$    t(8) = 1;
%$ end
%$
%$ if (isoctave && user_has_octave_forge_package('control')) || (~isoctave && user_has_matlab_license('control_toolbox'))
%$    if max(max(abs(Pstar1_large-Pstar3_large)))>1e-8
%$       t(9) = 0;
%$    else
%$       t(9) = 1;
%$    end
%$ else
%$    t(9) = 1;
%$ end
%$
%$ if max(max(abs(Pstar1_large-Pstar4_large)))>1e-8
%$    t(10) = 0;
%$ else
%$    t(10) = 1;
%$ end
%$
%$ T = all(t);
%@eof:1