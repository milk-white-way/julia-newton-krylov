module RHS

include("../form-boundary-conditions/form_bcs.jl")

using LinearAlgebra
using FormBCS

export rhs_newton

function rhs_newton(CSolution) 
  
  form_bcs(CSolution)


end

function rhs_calculation(CSolution)

  convection
  viscous_flux

end

function rhs_convection(CSolution, concoef)

  M = CSolution.m2
  N = CSolution.n2

  ucur_new_x = CSolution.ucur_x
  ucur_new_y = CSolution.ucur_y

  for runidx = 1:M-1
    for runidy = 1:N-1
      un = 0.5*ucur_new_x[runidx, runidy]
      um = un - abs(un)
      up = un + abs(un)
    end
  end

end
