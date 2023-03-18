using DrWatson
@quickactivate "IslandSizeRB"

using NCDatasets
using PyCall, LaTeXStrings
pplt = pyimport("proplot")

ds = NCDataset(datadir("samplerayleighbenard.nc"))
x = ds["xC"][:]
z = ds["zC"][:]
t = ds["time"][:]
b = ds["b"][:,1,:,:]
close(ds)

for it = 1 : length(t)

    pplt.close()
	fig,axs = pplt.subplots([[1],[1],[2]],aspect=4,axwidth=6,hspace=1)
	
	c = axs[1].pcolormesh(x,z,b[:,:,it]',levels=vcat(1:4,4.5,5.5,6:9)*10,extend="both",cmap="RdBu_r")
	# axs[1].format(xlim=(-200,200))
	fig.colorbar(c)

	axs[2].plot(x,b[:,1,it])
	axs[2].format(ylim=(0,100))
    fig.savefig(plotsdir("test-$it.png"),transparent=false,dpi=150)

end
