using DrWatson
@quickactivate "IslandSizeRB"

using Dates
using Logging
using Oceananigans
using Printf

processor = CPU()
# processor = GPU()

xmin,xmax = 0,6000; nx = 6000
ymin,ymax = -100,100; ny = 2000
zmin,zmax = 0,1; ny = 64

grid = RectilinearGrid(
    processor,
    size = (nx,ny),    x = (xmin,xmax), y = (ymin,ymax),
    topology=(Periodic, Flat, Bounded)
    # size = (nx,ny,nz), x = (xmin,xmax), y = (ymin,ymax), z = (zmin,zmax),
    # topology=(Periodic, Periodic, Bounded)
)

b_bc_top = ValueBoundaryCondition(0)
b_bc_bot = ValueBoundaryCondition(1)
b_bcs = FieldBoundaryConditions(top=b_bc_top, bottom=b_bc_bot)

const Ra = 1e7
const Pr = 1

# Ra = 1 / (κ² * Pr)
#   -> κ = sqrt(1 / (Ra * Pr))
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

wᵢ(x,y,z) = 1e-3 * (rand()-1); set!(model, w=wᵢ)

wizard = TimeStepWizard(cfl=0.7, max_change=1.1)

start_time = time_ns()
progress(sim) = @info "$(now()) - Oceananigans.jl - Integration completed through $(@sprintf("%07d",sim.model.clock.iteration)) steps | Model Time: $(@sprintf("%09.3f",sim.model.clock.time))"

simulation = Simulation(model, Δt=1e-2, stop_time=100)
simulation.callbacks[:wizard]   = Callback(wizard,   IterationInterval(10))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b        # unpack buoyancy `Field`

simulation.output_writers[:field_writer] = NetCDFOutputWriter(
    model, (; u, v, w, b),
    filename=datadir("more_fields.nc"),
    overwrite_existing=true,
    schedule=TimeInterval(0.1)
)

run!(simulation)