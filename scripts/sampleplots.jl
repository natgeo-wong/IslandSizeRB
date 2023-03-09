using DrWatson
@quickactivate "IslandSizeRB"

using NCDatasets
using PyCall, LaTeXStrings
pplt = pyimport("proplot")

ds = NCDataset(datadir("more_fields.nc"))
x = ds["xC"][:]
z = ds["zC"][:]
t = ds["time"][:]
b = ds["b"][:,1,:,:]
close(ds)

for it = 1 : length(t)

    pplt.close(); fig,axs = pplt.subplots(aspect=8,axwidth=4)
    axs[1].pcolormesh(x,z,b[:,:,it]',levels=vcat(3:10)./10)
    fig.savefig(plotsdir("test-$it.png"),transparent=false,dpi=150)

end
