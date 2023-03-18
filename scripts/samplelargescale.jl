using DrWatson
@quickactivate "IslandSizeRB"

using Dates
using Logging
using Oceananigans
using Printf

processor = CPU()
# processor = GPU()

@info "$(now()) - IslandSizeRB - Setting up grid"

xmin,xmax = -1500,1500; nx = 15000
ymin,ymax = -500,500; ny = 5000
zmin,zmax = 0,15; ny = 64

grid = RectilinearGrid(
    processor,
    size = (nx,ny),    x = (xmin,xmax), z = (zmin,zmax),
    topology=(Periodic, Flat, Bounded)
    # size = (nx,ny,nz), x = (xmin,xmax), z = (zmin,zmax), y = (ymin,ymax),
    # topology=(Periodic, Periodic, Bounded)
)

@info "$(now()) - IslandSizeRB - Setting up top and bottom buoyancy boundary conditions"

b_bc_top = ValueBoundaryCondition(0)
b_bc_bot = ValueBoundaryCondition(100)
b_bcs = FieldBoundaryConditions(top=b_bc_top, bottom=b_bc_bot)

const Ra = 1e8
const Pr = 1
const κ = sqrt(1 / (Ra * Pr))
const ν = κ * Pr

model = NonhydrostaticModel(;
    grid,
    timestepper = :RungeKutta3,
    advection = WENO(),
    tracers = :b,
    buoyancy = BuoyancyTracer(),
    boundary_conditions = (; b=b_bcs),
    closure = ScalarDiffusivity(; ν, κ)
)

@info "$(now()) - IslandSizeRB - Setting up initial buoyancy and vertical velocity conditions"

bᵢ(x,y,z) = 100 * (1-z)
wᵢ(x,y,z) = 1e-3 * (rand()-0.5)
set!(model, b=bᵢ, w=wᵢ)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)

start_time = time_ns()
progress(sim) = @info "$(now()) - Oceananigans.jl - Integration completed through $(@sprintf("%08.2f",sim.model.clock.time)) time units | Advective CFL: $(@sprintf("%0.3f",AdvectiveCFL(sim.Δt)(sim.model))) | Diffusive CFL: $(@sprintf("%0.3f",DiffusiveCFL(sim.Δt)(sim.model)))"

simulation = Simulation(model, Δt=1e-2, stop_time=100)

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b        # unpack buoyancy `Field`

simulation.output_writers[:fields] = NetCDFOutputWriter(
    model, (; u, v, w, b),
    indices = (2501:3500,:,:),
    filename=datadir("samplerayleighbenard.nc"),
    overwrite_existing=true,
    schedule=TimeInterval(0.1)
)

simulation.callbacks[:wizard]   = Callback(wizard,   IterationInterval(10))
simulation.callbacks[:progress] = Callback(progress, TimeInterval(1))

run!(simulation)