using DrWatson
@quickactivate "IslandSizeRB"

using Dates
using Logging
using NCDatasets
using Oceananigans
using Printf

#processor = CPU()
processor = GPU()

@info "$(now()) - IslandSizeRB - Setting up grid"

xmin,xmax = -500,500; nx = 5000
ymin,ymax = -500,500; ny = 5000
zmin,zmax = 0,15; ny = 75

grid = RectilinearGrid(
    processor,
    size = (nx,ny),    x = (xmin,xmax), z = (zmin,zmax),
    topology=(Periodic, Flat, Bounded)
    # size = (nx,ny,nz), x = (xmin,xmax), z = (zmin,zmax), y = (ymin,ymax),
    # topology=(Periodic, Periodic, Bounded)
)

@info "$(now()) - IslandSizeRB - Setting up top and bottom buoyancy boundary conditions"

b_bc_top = ValueBoundaryCondition(0)
@inline function b_bc_bot(x, y, t)
    if abs(x) < 50
        return 100 + 100 * sin(t*2*pi / 86.4)
    else 
        return 100
    end
end
b_bcs = FieldBoundaryConditions(top=b_bc_top, bottom=ValueBoundaryCondition(b_bc_bot))

const Ra = 1e8
const Pr = 1
const κ = sqrt(1 / (Ra * Pr))
const ν = κ * Pr

ds = NCDataset(datadir("samplelargescale.nc"))
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
progress(sim) = @info "$(now()) - Oceananigans.jl - Integration completed through $(@sprintf("%07.2f",sim.model.clock.time)) T | Adv CFL: $(@sprintf("%0.2f",AdvectiveCFL(sim.Δt)(sim.model))) | Diff CFL: $(@sprintf("%0.2f",DiffusiveCFL(sim.Δt)(sim.model)))"

simulation = Simulation(model, Δt=1e-2, stop_time=3.6*24*10)

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b        # unpack buoyancy `Field`

simulation.output_writers[:fields] = NetCDFOutputWriter(
    model, (; u, v, w, b),
    indices = (2001:3000,:,:),
    filename=datadir("sampleisland.nc"),
    overwrite_existing=true,
    schedule=TimeInterval(0.9)
)

simulation.callbacks[:wizard]   = Callback(wizard,   IterationInterval(10))
simulation.callbacks[:progress] = Callback(progress, TimeInterval(1))

run!(simulation)