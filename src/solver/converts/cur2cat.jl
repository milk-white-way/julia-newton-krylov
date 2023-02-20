module Converts

export cur2cat

function cur2cat(CSolution)

  M = length(ucur_x[:,1])
  N = length(ucur_x[1,:])

  ucur_x = CSolution.ucur_x
  ucur_y = CSolution.ucur_y

  for runidx = 2:M-1
    for runidy = 2:N-1
      ucat_x[runidx, runidy] = 0.5*( ucur_x[runidx, runidy] + ucur_x[runidx - 1, runidy] )
      ucat_y[runidx, runidy] = 0.5*( ucur_y[runidx, runidy] + ucur_y[runidx - 1, runidy] )
    end
  end

  for runidx = 1:M
    ucat_x[runidx, 1] = 0
    ucat_x[runidx, N] = 0

    ucat_y[runidx, 1] = 0
    ucat_y[runidx, N] = 0
  end

  for runidy = 1:N
    ucat_x[1, runidy] = 0
    ucat_x[M, runidy] = 0

    ucat_y[1, runidy] = 0
    ucat_y[M, runidy] = 0
  end

  global CSolution.ucat_x = ucat_x
  global CSolution.ucat_y = ycat_y

end

end
