using DrWatson
@quickactivate "IslandSizeRB"

using Oceananigans

grid = RectilinearGrid(
    # size=(256,256,64), x=(0,4), y=(0,4), z=(0,1),
    size=(1024,64), x=(0,256), z=(0,1),
    topology=(Periodic, Flat, Bounded)
)

b_bc_top = ValueBoundaryCondition(0)
b_bc_bot = ValueBoundaryCondition(1)
b_bcs = FieldBoundaryConditions(top=b_bc_top, bottom=b_bc_bot)

Ra = 1e7
Pr = 1

# Ra = 1 / (κ² * Pr)
#   -> κ = sqrt(1 / (Ra * Pr))
κ = sqrt(1 / (Ra * Pr))
ν = κ * Pr

model = NonhydrostaticModel(; grid,
                            timestepper = :RungeKutta3,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; b=b_bcs),
                            closure = ScalarDiffusivity(; ν, κ))

wᵢ(x, y, z) = 1e-3 * (rand() - 1)
set!(model, w=wᵢ)

simulation = Simulation(model, Δt=1e-2, stop_time=100)

wizard = TimeStepWizard(cfl=0.7, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

progress(sim) = @info string("Iter: ", iteration(sim),
                             ", time: ", time(sim),
                             ", max w: ", maximum(abs, model.velocities.w))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b        # unpack buoyancy `Field`

simulation.output_writers[:field_writer] = NetCDFOutputWriter(
    model, (; u, v, w, b), filename=datadir("more_fields.nc"), overwrite_existing=true,
    schedule=TimeInterval(0.1)
)

run!(simulation)