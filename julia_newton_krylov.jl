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
## GETTING LIBRARIES
using LinearAlgebra

## SCHEME CONSTANTS 
MAX = 100
TOL = 1.0e-5

## PROBLEM CONSTANTS
LENGTH_IN_X = pi
LENGTH_IN_Y = pi

DISCRETIZATION_IN_X = 61
DISCRETIZATION_IN_Y = 61

X_AFTER_GHOST = DISCRETIZATION_IN_X + 2
Y_AFTER_GHOST = DISCRETIZATION_IN_Y + 2

TOTALTIME = 10
TIMESTEP = 0.05
REN_NUM = 100
DYN_VIS = 1

CONVECTIVE_COEFFICIENT = 1/8

# dU_x=0 and du_Y=0
## INITIALIZATION STEP
# ++ Calculate bounds
# ++ Generate grid
# ++ Generate solution arrays on grid
# ++ Apply initial conditions (zero everywhere)
XSTEP = LENGTH_IN_X / (X_AFTER_GHOST - 3)
YSTEP = LENGTH_IN_Y / (Y_AFTER_GHOST - 3)

u_cart_x = zeros(m2, n2)
u_cart_y = zeros(m2, n2)
u_cont_x = zeros(m2, n2)
u_cont_y = zeros(m2, n2)
pressure = zeros(m2, n2)
u_bc_x = zeros(m2, n2)
u_bc_y = zeros(m2, n2)
time = zeros(TOTALTIME, 1)

mutable struct Solution
  du_x
  du_y
  u_cont_x
  u_cont_y
  u_cart_x
  u_cart_y
  u_bc_x
  u_bc_y
  pressure
  ren
  vis
  dx
  dy
  dt
  t
  m2
  n2
  u_init_x
  u_init_y
end

## ITERATIVE SCHEME
# ++ Steps here are:
# ====== Solving the Momentum equation => Get velocity field
# ====== Solving the Poisson equation => Get pressure fied
# ====== Update the solution
# ====== Apply Boundary Conditions
# ====== Repeat until reaching desired tolerance
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
## BCS
function forming_bcs(sol_field::Solution)
    # Fetching contravariant velocity field
    u_cont_new_x = sol_field.u_cont_x
    u_cont_new_y = sol_field.u_cont_y
    # Setup x contravariant component
    for runix in 1:sol_field.n2
        u_cont_new_x[1, runid] = 0
        u_cont_new_x[sol_field.m2-1, runix] = 0
    end
    # Setup y contravariant component
    for runix in 1:sol_field.n2
        u_cont_new_y[runix, 1] = 0
        u_cont_new_x[runix, sol_field.n2-1] = 0
    end
    # Zeroing the boundary (no-slip boundary condition)
    for runix in 1:sol_field.n2
        for runiz in 1:sol_field.m2
            sol_field.u_bc_x[runix, runiz] = 0
            sol_field.u_bc_y[runix, runiz] = 0
        end
    end
end
# Taylor-Green
# NEED TO UNDERSTAND ALGORITHM