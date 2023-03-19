using DrWatson
@quickactivate "IslandSizeRB"

using Dates
using Logging
using NCDatasets
using Oceananigans
using Printf

processor = CPU()
#processor = GPU()

@info "$(now()) - IslandSizeRB - Setting up grid"

xmin,xmax = -2,2; nx = 200
zmin,zmax = 0,1; ny = 100

grid = RectilinearGrid(
    processor,
    size = (nx,ny),    x = (xmin,xmax), z = (zmin,zmax),
    topology=(Periodic, Flat, Bounded)
)

@info "$(now()) - IslandSizeRB - Setting up top and bottom buoyancy boundary conditions"

b_bc_top = ValueBoundaryCondition(0)
@inline b_bc_bot(x, y, t) = 100 + 10 * sin(t*10*pi);
b_bcs = FieldBoundaryConditions(top=b_bc_top, bottom=ValueBoundaryCondition(b_bc_bot))

const Ra = 1e8
const Pr = 1
const κ = sqrt(1 / (Ra * Pr))
const ν = κ * Pr

@info "$(now()) - IslandSizeRB - Loading up initial conditions from spinup"

ds = NCDataset(datadir("samplerayleighbenard.nc"))
const uᵢ = ds["u"][:,:,:,end]
const vᵢ = ds["v"][:,:,:,end]
const wᵢ = ds["w"][:,:,:,end]
const bᵢ = ds["b"][:,:,:,end]
close(ds)

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

set!(model, b=bᵢ, u=uᵢ, v=vᵢ, w=wᵢ)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)

start_time = time_ns()
progress(sim) = @info "$(now()) - Oceananigans.jl - Integration completed through $(@sprintf("%08.2f",sim.model.clock.time)) time units | Advective CFL: $(@sprintf("%0.3f",AdvectiveCFL(sim.Δt)(sim.model))) | Diffusive CFL: $(@sprintf("%0.3f",DiffusiveCFL(sim.Δt)(sim.model)))"

simulation = Simulation(model, Δt=1e-3, stop_time=2)

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b        # unpack buoyancy `Field`

simulation.output_writers[:fields] = NetCDFOutputWriter(
    model, (; b),
    filename=datadir("samplediurnal.nc"),
    overwrite_existing=true,
    schedule=TimeInterval(0.01)
)

simulation.callbacks[:wizard]   = Callback(wizard,   IterationInterval(10))
simulation.callbacks[:progress] = Callback(progress, TimeInterval(0.01))

run!(simulation)