module Initialization

export init

function init(CSolution, SYS2D_LENGTH_X, SYS2D_LENGTH_Y,
              SYS2D_GHOST_X, SYS2D_GHOST_Y,
              SYS2D_STEP_X, SYS2D_STEP_Y,
              STEP_TIME, ENDO_TIME,
              REN_NUM, DYN_VIS)


    Ucat_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Ucat_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Ucur_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Ucur_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Ubcs_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Ubcs_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Uini_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Uini_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
      Pres = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)

    Conv_Flux_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Conv_Flux_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Visc_Flux_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Visc_Flux_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Pres_Grad_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    Pres_Grad_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
  
    RHS_x = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)
    RHS_y = Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y)

    for runidx = 1:SYS2D_GHOST_X
        for runidy = 1:SYS2D_GHOST_Y
            Ucat_x[runidx, runidy] = 0
            Ucat_y[runidx, runidy] = 0
            Ucur_x[runidx, runidy] = 0
            Ucur_y[runidx, runidy] = 0
            Ubcs_x[runidx, runidy] = 0
            Ubcs_y[runidx, runidy] = 0
            Uini_x[runidx, runidy] = 0
            Uini_y[runidx, runidy] = 0
             Press[runidx, runidy] = 0

            Conv_Flux_x[runidx, runidy] = 0
            Conv_Flux_y[runidx, runidy] = 0
            Visc_Flux_x[runidx, runidy] = 0
            Visc_Flux_y[runidx, runidy] = 0
            Pres_Grad_x[runidx, runidy] = 0
            Pres_Grad_y[runidx, runidy] = 0
            
            RHS_x[runidx, runidy] = 0
            RHS_y[runidx, runidy] = 0
        end
    end

    global CSolution.length_x = SYS2D_LENGTH_X
    global CSolution.length_y = SYS2D_LENGTH_Y
    global CSolution.m2 = SYS2D_GHOST_X
    global CSolution.n2 = SYS2D_GHOST_Y
    global CSolution.dx = SYS2D_STEP_X
    global CSolution.dy = SYS2D_STEP_Y
    global CSolution.dt = STEP_TIME
    global CSolution.et = ENDO_TIME
    global CSolution.ren = REN_NUM
    global CSolution.vis = DYN_VIS
    global CSolution.ucat_x = Ucat_x
    global CSolution.ucat_y = Ucat_y
    global CSolution.ucur_x = Ucur_x
    global CSolution.ucur_y = Ucur_y
    global CSolution.ubcs_x = Ubcs_x
    global CSolution.ubcs_y = Ubcs_y
    global CSolution.uini_x = Uini_x
    global CSolution.uini_y = Uini_y
    global CSolution.pres   = Pres

    global CSolution.conv_flux_x = Conv_Flux_x
    global CSolution.conv_flux_y = Conv_Flux_y
    global CSolution.visc_flux_x = Visc_Flux_x
    global CSolution.visc_flux_y = Visc_Flux_y
    global CSolution.pres_grad_x = Pres_Grad_x
    global CSolution.pres_grad_y = Pres_Grad_y

    global CSolution.rhs_x = RHS_x
    global CSolution.rhs_y = RHS_y
end
end
