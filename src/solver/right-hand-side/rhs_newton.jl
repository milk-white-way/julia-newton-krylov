module RHS

include("../form-boundary-conditions/form_bcs.jl")

using LinearAlgebra
using FormBCS

export rhs_newton

function rhs_newton(CSolution) 
  
  M = CSolution.m2
  N = CSolution.n2

  dt = CSolution.dt

  form_bcs(CSolution)
  rhs_calc(CSolution)

  local_rhs_x = CSolution.rhs_x
  local_rhs_y = CSolution.rhs_y

  local_ucur_x = CSolution.ucur_x
  local_ucur_y = CSolution.ucur_y

  local_rhs_x = local_rhs_x - ( 1.5/dt )*( uim_x - local_ucur_x ) + ( 0.5/dt )*(  )
  local_rhs_y = local_rhs_y - ( 1.5/dt )*( uim_y - local_ucur_y ) + ( 0.5/dt )*(  )

  for runidx = 1:M
    local_rhs_y[runidx,   1] = 0
    local_rhs_y[runidx,   N] = 0
    local_rhs_y[runidx, N-1] = 0
  end

  for runidy = 1:N
    local_rhs_x[  1, runidy] = 0
    local_rhs_x[  M, runidy] = 0
    local_rhs_x[M-1, runidy] = 0
  end

  global CSolution.rhs_x = local_rhs_x
  global CSolution.rhs_y = local_rhs_y

end

function rhs_calc(CSolution)

  M = CSolution.m2
  N = CSolution.n2

  rhs_convective_flux_calc(CSolution)
  rhs_viscous_flux_calc(CSolution)
  rhs_pressure_gradient_calc(CSolution)

  local_conv_flux_x = CSolution.conv_flux_x
  local_conv_flux_y = CSolution.conv_flux_y
  local_visc_flux_x = CSolution.visc_flux_x
  local_visc_flux_y = CSolution.visc_flux_y
  local_pres_grad_x = CSolution.pres_grad_x
  local_visc_flux_y = CSolution.pres_grad_y

  local_rhs_x = CSolution.rhs_x
  local_rhs_y = CSolution.rhs_y

  total_flux_x = - local_conv_flux_x + local_visc_flux_x - local_pres_grad_x
  total_flux_y = - local_conv_flux_y + local_visc_flux_y - local_pres_grad_y

  for runidx = 1:M-1
    for runidy = 1:N
      local_rhs_x[runidx, runidy] = 0.5*( total_flux_x[runidx, runidy] + total_flux_x[runidx+1, runidy] )
    end
  end

  for runidx = 1:M
    for runidy = 1:N-1
      local_rhs_y[runidx, runidy] = 0.5*( total_flux_y[runidx, runidy] + total_flux_y[runidx, runidy+1] )
    end
  end

  for runidx = 1:M
    local_rhs_x[runidx,   1] = 0
    local_rhs_x[runidx,   N] = 0

    local_rhs_y[runidx,   1] = 0
    local_rhs_y[runidx, N-1] = 0
    local_rhs_y[runidx,   N] = 0
  end

  for runidy = 1:N
    local_rhs_x[  M, runidy] = 0
    local_rhs_x[M-1, runidy] = 0
    local_rhs_x[  1, runidy] = 0

    local_rhs_y[  1, runidy] = 0
    local_rhs_y[  M, runidy] = 0
  end

  global CSolution.rhs_x = local_rhs_x
  global CSolution.rhs_y = local_rhs_y

end

function rhs_convective_flux_calc(CSolution, coef)

  M = CSolution.m2
  N = CSolution.n2

  dx = CSolution.dx
  dy = CSolution.dy

  local_ucat_x = CSolution.ucat_x
  local_ucat_y = CSolution.ucat_y

  local_ucur_x = CSolution.ucur_x
  local_ucur_y = CSolution.ucur_y

  local_conv_flux_x = CSolution.conv_flux_x
  local_conv_flux_y = CSolution.conv_flux_y

  for runidx = 1:M-1
    for runidy = 1:N
      un = 0.5*local_ucur_x[runidx, runidy]
      um = un - abs(un)
      up = un + abs(un)

      vn = 0.5*local_ucur_y[runidx, runidy]
      vm = vn - abs(vn)
      vp = vn + abs(vn)

      if ( runidx == 1 )

        fpx1[runidx,  runidy] = um*( coef*( - local_ucat_x[runidx+2, runidy] - 2*local_ucat_x[runidx+1, runidy] + 3*local_ucat_x[runidx, runidy] ) + local_ucat_x[runidx+1, runidy] ) + up*( coef*( - local_ucat_x[runidx, runidy] - 2*local_ucat_x[runidx, runidy] + 3*local_ucat_x[runidx+1, runidy] ) + local_ucat_x[runidx, runidy] )

        fpy1[runidx,  runidy] = vm*( coef*( - local_ucat_y[runidx+2, runidy] - 2*local_ucat_y[runidx+1, runidy] + 3*local_ucat_y[runidx, runidy] ) + local_ucat_y[runidx+1, runidy] ) + vp*( coef*( - local_ucat_y[runidx, runidy] - 2*local_ucat_y[runidx, runidy] + 3*local_ucat_y[runidx+1, runidy] ) + local_ucat_y[runidx, runidy] )

      elseif ( runidx == M-1 )

        fpx1[runidx,  runidy] = um*( coef*( - local_ucat_x[runidx+1, runidy] - 2*local_ucat_x[runidx+1, runidy] + 3*local_ucat_x[runidx, runidy] ) + local_ucat_x[runidx+1, runidy] ) + up*( coef*( - local_ucat_x[runidx-1, runidy] - 2*local_ucat_x[runidx, runidy] + 3*local_ucat_x[runidx+1, runidy] ) + local_ucat_x[runidx, runidy] )

        fpy1[runidx,  runidy] = vm*( coef*( - local_ucat_y[runidx+1, runidy] - 2*local_ucat_y[runidx+1, runidy] + 3*local_ucat_y[runidx, runidy] ) + local_ucat_y[runidx+1, runidy] ) + vp*( coef*( - local_ucat_y[runidx-1, runidy] - 2*local_ucat_y[runidx, runidy] + 3*local_ucat_y[runidx+1, runidy] ) + local_ucat_y[runidx, runidy] )

      else

        fpx1[runidx,  runidy] = um*( coef*( - local_ucat_x[runidx+2, runidy] - 2*local_ucat_x[runidx+1, runidy] + 3*local_ucat_x[runidx, runidy] ) + local_ucat_x[runidx+1, runidy] ) + up*( coef*( - local_ucat_x[runidx-1, runidy] - 2*local_ucat_x[runidx, runidy] + 3*local_ucat_x[runidx+1, runidy] ) + local_ucat_x[runidx, runidy] )

        fpy1[runidx,  runidy] = vm*( coef*( - local_ucat_y[runidx+2, runidy] - 2*local_ucat_y[runidx+1, runidy] + 3*local_ucat_y[runidx, runidy] ) + local_ucat_y[runidx+1, runidy] ) + vp*( coef*( - local_ucat_y[runidx-1, runidy] - 2*local_ucat_y[runidx, runidy] + 3*local_ucat_y[runidx+1, runidy] ) + local_ucat_y[runidx, runidy] )

      end

    end
  end

  for runidy = 1:N
    fpx1[M, runidy] = 0
    fpy1[M, runidy] = 0
  end

  for runidx = 1:M
    for runidy = 1:N-1
      un = 0.5*local_ucur_x[runidx, runidy]
      um = un - abs(un)
      up = un + abs(un)

      vn = 0.5*local_ucur_y[runidx, runidy]
      vm = vn - abs(vn)
      vp = vn + abs(vn)

      if ( runidy == 1 )

        fpx2[runidx,  runidy] = um*( coef*( - local_ucat_x[runidx, runidy+2] - 2*local_ucat_x[runidx, runidy+1] + 3*local_ucat_x[runidx, runidy] ) + local_ucat_x[runidx, runidy+1] ) + up*( coef*( - local_ucat_x[runidx, runidy] - 2*local_ucat_x[runidx, runidy] + 3*local_ucat_x[runidx, runidy+1] ) + local_ucat_x[runidx, runidy] )

        fpy2[runidx,  runidy] = vm*( coef*( - local_ucat_y[runidx, runidy+2] - 2*local_ucat_y[runidx, runidy+1] + 3*local_ucat_y[runidx, runidy] ) + local_ucat_y[runidx, runidy+1] ) + vp*( coef*( - local_ucat_y[runidx, runidy] - 2*local_ucat_y[runidx, runidy] + 3*local_ucat_y[runidx, runidy+1] ) + local_ucat_y[runidx, runidy] )

      elseif ( runidy == N-1 )

        fpx2[runidx,  runidy] = um*( coef*( - local_ucat_x[runidx, runidy+1] - 2*local_ucat_x[runidx, runidy+1] + 3*local_ucat_x[runidx, runidy] ) + local_ucat_x[runidx, runidy+1] ) + up*( coef*( - local_ucat_x[runidx, runidy-1] - 2*local_ucat_x[runidx, runidy] + 3*local_ucat_x[runidx, runidy+1] ) + local_ucat_x[runidx, runidy] )

        fpy2[runidx,  runidy] = vm*( coef*( - local_ucat_y[runidx, runidy+1] - 2*local_ucat_y[runidx, runidy+1] + 3*local_ucat_y[runidx, runidy] ) + local_ucat_y[runidx, runidy+1] ) + vp*( coef*( - local_ucat_y[runidx, runidy-1] - 2*local_ucat_y[runidx, runidy] + 3*local_ucat_y[runidx, runidy+1] ) + local_ucat_y[runidx, runidy] )

      else

        fpx2[runidx,  runidy] = um*( coef*( - local_ucat_x[runidx, runidy+2] - 2*local_ucat_x[runidx, runidy+1] + 3*local_ucat_x[runidx, runidy] ) + local_ucat_x[runidx, runidy+1] ) + up*( coef*( - local_ucat_x[runidx, runidy-1] - 2*local_ucat_x[runidx, runidy] + 3*local_ucat_x[runidx, runidy+1] ) + local_ucat_x[runidx, runidy] )

        fpy2[runidx,  runidy] = vm*( coef*( - local_ucat_y[runidx, runidy+2] - 2*local_ucat_y[runidx, runidy+1] + 3*local_ucat_y[runidx, runidy] ) + local_ucat_y[runidx, runidy+1] ) + vp*( coef*( - local_ucat_y[runidx, runidy-1] - 2*local_ucat_y[runidx, runidy] + 3*local_ucat_y[runidx, runidy+1] ) + local_ucat_y[runidx, runidy] )

      end

    end
  end

  for runidx = 1:M
    fpx2[runidx, N] = 0
    fpy2[runidx, N] = 0
  end

  for runidx = 2:M-1
    for runidy = 2:N-1
      local_conv_flux_x[runidx, runidy] = ( fpx1[runidx, runidy] - fpx1[runidx-1, runidy] )/dx + ( fpx2[runidx, runidy] - fpx2[runidx, runidy-1] )/dy
      local_conv_flux_y[runidx, runidy] = ( fpy1[runidx, runidy] - fpy1[runidx-1, runidy] )/dx + ( fpy2[runidx, runidy] - fpy2[runidx, runidy-1] )/dy
    end
  end

  for runidx = 1:M
    local_conv_flux_x[runidx, 1] = 0
    local_conv_flux_y[runidx, 1] = 0
    local_conv_flux_x[runidx, N] = 0
    local_conv_flux_y[runidx, N] = 0
  end

  for runidy = 1:N
    local_conv_flux_x[1, runidy] = 0
    local_conv_flux_y[1, runidy] = 0
    local_conv_flux_x[M, runidy] = 0
    local_conv_flux_y[M, runidy] = 0
  end

  global CSolution.conv_flux_x = local_conv_flux_x
  global CSolution.conv_flux_y = local_conv_flux_y

end

function rhs_viscous_flux_calc(CSolution)

  M = CSolution.m2
  N = CSolution.n2

  dx = CSolution.dx
  dy = CSolution.dy

  local_ucat_x = CSolution.ucat_x
  local_ucat_y = CSolution.ucat_y

  local_ucur_x = CSolution.ucur_x
  local_ucur_y = CSolution.ucur_y

  local_visc_flux_x = CSolution.visc_flux_x
  local_visc_flux_y = CSolution.visc_flux_y

  local_ren = CSolution.ren

  for runidx = 2:M-1
    for runidy = 2:N-1
      ue = local_ucat_x[runidx-1, runidy]
      up = local_ucat_x[  runidx, runidy]
      uw = local_ucat_x[runidx+1, runidy]

      un = local_ucat_x[runidx, runidy-1]
      us = local_ucat_x[runidx, runidy+1]

      ve = local_ucat_y[runidx-1, runidy]
      vp = local_ucat_y[  runidx, runidy]
      vw = local_ucat_y[runidx+1, runidy]

      vn = local_ucat_y[runidx, runidy-1]
      vs = local_ucat_y[runidx, runidy+1]

      local_visc_flux_x[runidx, runidy] = ( uw - 2*up + ue )./( dx*dx ) + ( us - 2*up + un )./( dy*dy )
      local_visc_flux_y[runidx, runidy] = ( vw - 2*vp + ve )./( dx*dx ) + ( vs - 2*vp + vn )./( dy*dy )
    end
  end

  for runidx = 1:M
    local_visc_flux_x[runidx, 1] = 0
    local_visc_flux_y[runidx, 1] = 0
    local_visc_flux_x[runidx, N] = 0
    local_visc_flux_y[runidx, N] = 0
  end
  
  for runidy = 1:N
    local_visc_flux_x[1, runidy] = 0
    local_visc_flux_y[1, runidy] = 0
    local_visc_flux_x[1, runidy] = 0
    local_visc_flux_y[1, runidy] = 0
  end

  global CSolution.visc_flux_x = local_visc_flux_x/local_ren
  global CSolution.visc_flux_y = local_visc_flux_y/local_ren

end

function rhs_pressure_gradient_calc(CSolution)

  M = CSolution.m2
  N = CSolution.n2

  local_pres_grad_x = CSolution.pres_grad_x
  local_pres_grad_y = CSolution.pres_grad_y

  local_pres = CSolution.pres

  for runidx = 1:M
    for runidy = 1:N
      if ( runidx == 1 ) 

        local_pres_grad_x[runidx, runidy] = 0

      elseif ( runidx == 2 )

        local_pres_grad_x[runidx, runidy] = ( local_pres[runidx+1, runidy] - local_pres[runidx, runidy] ) / ( dx )

      elseif ( runidx == M-1 )

        local_pres_grad_x[runidx, runidy] = ( local_pres[runidx, runidy] - local_pres[runidx-1, runidy] ) / ( dx )

      elseif ( runidx == M )

        local_pres_grad_x[runidx, runidy] = 0

      else

        local_pres_grad_x[runidx, runidy] = ( local_pres[runidx+1, runidy] - local_pres[runidx-1, runidy] ) / ( 2*dx )

      end


      if ( runidy == 1 ) 

        local_pres_grad_y[runidx, runidy] = 0

      elseif ( runidy == 2 )

        local_pres_grad_y[runidx, runidy] = ( local_pres[runidx, runidy+1] - local_pres[runidx, runidy] ) / ( dy )

      elseif ( runidy == N-1 )

        local_pres_grad_y[runidx, runidy] = ( local_pres[runidx, runidy] - local_pres[runidx, runidy-1] ) / ( dy )

      elseif ( runidy == N )

        local_pres_grad_y[runidx, runidy] = 0

      else

        local_pres_grad_y[runidx, runidy] = ( local_pres[runidx, runidy+1] - local_pres[runidx, runidy-1] ) / ( 2*dy )

      end
    end
  end

  global CSolution.pres_grad_x = local_pres_grad_x
  global CSolution.pres_grad_y = local_pres_grad_y

end

end
