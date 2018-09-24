# @time using Plots
# plotly() # Choose the Plotly.jl backend for web interactivity
# plot(rand(5,5),linewidth=2,title="My Plot")


using PyPlot
const DATA_PATH = "/tmp/dump/data/"

x, y = np.loadtxt('example.txt', delimiter=',', unpack=True)
plt.plot(x,y, label='Loaded from file!')
