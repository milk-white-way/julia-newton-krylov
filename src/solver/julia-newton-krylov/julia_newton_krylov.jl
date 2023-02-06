## JULIAN CODE BY TAM THIEN NGUYEN
# ++ email: tam.thien.nguyen@ndus.edu
## The goal is to numerically solve the Navier-Stokes equation
# Assumption:
# ++ Fluid Mechanics:
# ====== Incompressible
# ====== <Item 2>
# ++ Fluid Dynamics:
# ====== 2 Dimensional
# ====== <Item 1>
# ++ Input
# ====== <Add input here>
# ++ Output
# ====== u_cart_x: x-component of the Cartesian velocity
# ====== u_cart_y: y-component of the Cartesian velocity
# ====== pressure: scalar pressure field
# ====== time: time

include("../initialization/init.jl")
## GETTING LIBRARIES
using LinearAlgebra
using Debugger
using .Initialization

## SCHEME CONSTANTS
const MAX::UInt8 = 100
const TOL::Float64 = 1.0e-5

## PROBLEM CONSTANTS
const SYS2D_LENGTH_X::Float64 = pi
const SYS2D_LENGTH_Y::Float64 = pi
const SYS2D_DIM_X::UInt8 = 61
const SYS2D_DIM_Y::UInt8 = 61
const SYS2D_GHOST_X::UInt8 = SYS2D_DIM_X + 2
const SYS2D_GHOST_Y::UInt8 = SYS2D_DIM_Y + 2
const SYS2D_STEP_X::Float64 = SYS2D_LENGTH_X / (SYS2D_DIM_X - 1)
const SYS2D_STEP_Y::Float64 = SYS2D_LENGTH_Y / (SYS2D_DIM_Y - 1)

const STEP_TIME::Float64 = 0.05
const ENDO_TIME::UInt8 = 10

const REN_NUM::UInt8 = 100
const DYN_VIS::UInt8 = 1

CONVECTIVE_COEFFICIENT::Float64 = 0.128 # 1/8

# Creating solution structure containing:
#   Constants: lengths, spatial dimensions, spatial steps, time step, end time, reynolds number, viscosity
#   Variables: cartesian velocity field, curvilinear velocity field, boundary conditions, initial velocity field (in cartesian), pressure field
#mutable struct Constant2D
mutable struct Solution2D
    length_x::Float64
    length_y::Float64
    m2::UInt8
    n2::UInt8
    dx::Float64
    dy::Float64
    dt::Float64
    et::UInt8
    ren::UInt8
    vis::Float64
#end

#mutable struct Variable2D
    ucat_x::Matrix{Float64}
    ucat_y::Matrix{Float64}
    ucur_x::Matrix{Float64}
    ucur_y::Matrix{Float64}
    ubcs_x::Matrix{Float64}
    ubcs_y::Matrix{Float64}
    uini_x::Matrix{Float64}
    uini_y::Matrix{Float64}
     press::Matrix{Float64}
#end

#mutable struct Solution2D
#    con::Constant2D
#    var::Variable2D
end

CSolution = Solution2D(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y),
                       Matrix{Float64}(undef, SYS2D_GHOST_X, SYS2D_GHOST_Y))


init(CSolution,
     SYS2D_LENGTH_X, SYS2D_LENGTH_Y,
     SYS2D_GHOST_X, SYS2D_GHOST_Y,
     SYS2D_STEP_X, SYS2D_STEP_Y,
     STEP_TIME, ENDO_TIME,
     REN_NUM, DYN_VIS)

#=
for steptime in 1:TOTALTIME
  # Calculate current time
    global time = steptime*TIMESTEP
  # Initialize previous solution
  Pstep = Solution(0, 0, u_cont_x, u_cont_y, 
                   u_cart_x, u_cart_y, 
                   u_bc_x, u_bc_y,
                   pressure, 
                   REN_NUM,
                   DYN_VIS,
                   XSTEP,
                   YSTEP,
                   TIMESTEP,
                   time,
                   X_AFTER_GHOST,
                   Y_AFTER_GHOST,
                   u_cont_x,
                   u_cont_y
                )
  # Forming Boundary Conditions
  # Calculate RHS vector
  ## Convective Flux
    for runix in 1:Pstep.m2
        for runiz in 1:Pstep.n2
            u_conv = 0.5 * Pstep.u_cont_x[runix, runiz]

            up = u_conv + abs(u_conv) #???
            um = u_conv - abs(u_conv) #???


        end
    end
  ## Viscous Flux

  ## Pressure Gradient

  # Solving the Momentum Equation using Jacobian-free Newton Krylov

end
=#

# Taylor-Green
# NEED TO UNDERSTAND ALGORITHM
