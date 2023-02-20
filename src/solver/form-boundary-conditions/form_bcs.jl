module FormBCS

include("../converts/cur2cat.jl")
using .Converts

export form_bcs

function form_bcs(CSolution)

  dx = CSolution.dx
  dy = CSolution.dy

  ucur_new_x = CSolution.ucur_x
  ucur_new_y = CSolution.ucur_y

  M = CSolution.m2
  N = CSolution.n2

  for runidx = 1:M
    ucur_new_x[runidx,   1] = 0
    ucur_new_y[runidx, N-1] = 0

    jj = 1
    x = (runidx - 2)*dx
    y = (jj - 2 + 0.5)*dy
    ubcs_x_1[runidx] = 0
    ubcs_y_1[runidx] = 0

    jj = N-1
    x = (runidx - 2)*dx
    y = (jj - 2 + 0.5)*dy
    ubcs_x_2[runidx] = 0
    ubcs_y_2[runidx] = 0
  end

  for runidy = 1:N
    ucur_new_x[  1, runidy] = 0
    ucur_new_y[M-1, runidy] = 0
    
    ii = 1
    x = (ii - 2 + 0.5)*dx
    y = (runidy - 2)*dy
    ubcs_x_3[runidy] = 0
    ubcs_y_3[runidy] = 0

    ii = M-1
    x = (ii - 2 + 0.5)*dx
    y = (runidy - 2)*dy
    ubcs_x_4[runidy] = 0
    ubcs_y_4[runidy] = 0
  end 

  global CSolution.ucur_x = ucur_new_x
  global CSolution.ucur_y = ucur_new_y

  cur2cat(CSolution)

  ucat_new_x = CSolution.ucat_x
  ucat_new_y = CSolution.ucat_y

  for runidx = 1:M
    ucat_new_x[runidx, 1] = 2 * ubcs_x_1[runidx] - ucat_new_x[runidx, 2]
    ucat_new_y[runidx, 1] = 2 * ubcs_y_1[runidx] - ucat_new_y[runidx, 2]

    ucat_new_x[runidx, N] = 2 * ubcs_x_2[runidx] - ucat_new_x[runidx, N-]
    ucat_new_y[runidx, N] = 2 * ubcs_y_2[runidx] - ucat_new_y[runidx, N-]
  end

  for runidy = 1:N
    ucat_new_x[1, runidy] = 2 * ubcs_x_3[runidy] - ucat_new_x[2, runidy]
    ucat_new_y[1, runidy] = 2 * ubcs_y_3[runidy] - ucat_new_y[2, runidy]

    ucat_new_x[M, runidy] = 2 * ubcs_x_4[runidy] - ucat_new_x[M-1, runidy]
    ucat_new_y[M, runidy] = 2 * ubcs_y_4[runidy] - ucat_new_y[M-1, runidy]
  end

  global CSolution.ucat_x = ucat_new_x
  global CSolution.ucat_y = ucat_new_y

end
