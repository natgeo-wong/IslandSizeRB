using DrWatson
@quickactivate "IslandSizeRB"

using NCDatasets
using StatsBase
using PyCall, LaTeXStrings
pplt = pyimport("proplot")

ds = NCDataset(datadir("samplediurnal.nc"))
x = ds["xC"][:]
z = ds["zC"][:]
b = ds["b"][:,1,:,:]
close(ds)

for it = 401 : 2 : 1001

    pplt.close()
	fig,axs = pplt.subplots([[1],[1],[1],[1],[2]],aspect=4,axwidth=4,hspace=1)
	
	c = axs[1].pcolormesh(x,z,b[:,:,it]',levels=vcat(40:2:48,49,51,52:2:60),extend="both",cmap="RdBu_r")
	# axs[1].format(xlim=(-200,200))
	fig.colorbar(c)

	pn = axs[1].panel("r",width="3em")
	pn.plot(dropdims(mean(b[:,:,it],dims=1),dims=1),z)
	pn.format(xlim=(40,60),xlocator=50)

	axs[2].plot(x,b[:,1,it])
	axs[2].format(ylim=(0,100))
    fig.savefig(plotsdir("test-$it.png"),transparent=false,dpi=150)

end
