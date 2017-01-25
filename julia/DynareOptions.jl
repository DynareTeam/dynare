module DynareOptions
##
 # Copyright (C) 2015 Dynare Team
 #
 # This file is part of Dynare.
 #
 # Dynare is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Dynare is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
##

export Options, dynare_options

type PFMSolver
    maxit::Int
    periods::Int
    tolx::Float64
    tolf::Float64
end

function pfmsolver_set_defaults()
    return PFMSolver(500,       # maxit (Maximum number of iterations in Newton algorithm)
                     400,       # periods (Number of periods to return to the steady state)
                     1e-6,      # tolx (Tolerance criterion on the paths for the endogenous variables)
                     1e-6       # tolf (Tolerance criterion on the stacked non linear equations)
                     )
end

type Options
    dynare_version::String
    linear::Bool
    pfmsolver::PFMSolver
end

function dynare_options()
    return Options("",                          # dynare_version
                   false,                       # linear
                   pfmsolver_set_defaults()     # pfmsolver
                   )
end

end
