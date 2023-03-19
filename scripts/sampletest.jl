using DrWatson
@quickactivate "IslandSizeRB"

using Dates
using Logging
using Oceananigans
using Printf

processor = CPU()
#processor = GPU()

@info "$(now()) - IslandSizeRB - Setting up grid"

xmin,xmax = -30,30; nx = 300
ymin,ymax = -500,500; ny = 5000
zmin,zmax = 0,15; ny = 75

grid = RectilinearGrid(
    processor,
    size = (nx,ny),    x = (xmin,xmax), z = (zmin,zmax),
    topology=(Periodic, Flat, Bounded)
)

@info "$(now()) - IslandSizeRB - Setting up top buoyancy boundary conditions"

b_bc_bot = ValueBoundaryCondition(100)
b_bcs = FieldBoundaryConditions(top=b_bc_top)

@inline function b_forcing_func(x, y, z, t, b)
    if b > 0.05
        return -1.5 / 86.4
    else
        return -b / 432
    end
end
b_forcing = Forcing(b_forcing_func, field_dependencies=(:b))

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
    forcing=(b=b_forcing)
    closure = ScalarDiffusivity(; ν, κ)
)

@info "$(now()) - IslandSizeRB - Setting up initial buoyancy and vertical velocity conditions"

bᵢ(x,y,z) = 0.5 + 0.5 * (15-2*z) / 15
wᵢ(x,y,z) = 1e-3 * (rand()-0.5)
set!(model, b=bᵢ, w=wᵢ)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)

start_time = time_ns()
progress(sim) = @info "$(now()) - Oceananigans.jl - Integration completed through $(@sprintf("%08.2f",sim.model.clock.time)) time units | Advective CFL: $(@sprintf("%0.3f",AdvectiveCFL(sim.Δt)(sim.model))) | Diffusive CFL: $(@sprintf("%0.3f",DiffusiveCFL(sim.Δt)(sim.model)))"

simulation = Simulation(model, Δt=1e-2, stop_time=500)

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b        # unpack buoyancy `Field`

simulation.output_writers[:fields] = NetCDFOutputWriter(
    model, (; b),
    filename=datadir("sampletest.nc"),
    overwrite_existing=true,
    schedule=TimeInterval(1)
)

simulation.callbacks[:wizard]   = Callback(wizard,   IterationInterval(10))
simulation.callbacks[:progress] = Callback(progress, TimeInterval(0.01))

run!(simulation)