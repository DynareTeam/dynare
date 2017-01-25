module DynareOutput
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

export Ouput, dynare_output

type Output
    dynare_version::String
    steady_state::Vector{Float64}
    exo_steady_state::Vector{Float64}
end

function dynare_output()
    return Output("",                # dynare_version
                  Array(Float64, 0), # steady_state
                  Array(Float64, 0)  # exo_steady_state
                 )
end

end
