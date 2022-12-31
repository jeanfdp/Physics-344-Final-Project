using DelimitedFiles: readdlm, writedlm
using ProgressMeter
include("disjointSetGraphStruct.jl")


function findCorrLength!(lengths::Vector{Float64},n::Int,β::Float64,μ::Float64,tmpGraphs::Vector{disjointSetGraph},tmpNews::Vector{Vector{Float64}},edges::Vector{Vector{Int}},otherCords::Vector{Vector{Int}},cords)
    p=1/((1+exp(β*(μ)))^4 * (1+exp(β*(μ-1)))^4)
    Threads.@threads for i=eachindex(lengths)
        resetDisjointSetGraph!(tmpGraphs[Threads.threadid()],cords)
        otherCords[Threads.threadid()][1]=0
        otherCords[Threads.threadid()][2]=0
        for j=1:n^2
            if p>rand()
                for edge in edges
                    otherCords[Threads.threadid()][1]=cords[j][1]+edge[1]
                    otherCords[Threads.threadid()][2]=cords[j][2]+edge[2]
                    if otherCords[Threads.threadid()][1]>0 && otherCords[Threads.threadid()][1]<=n && otherCords[Threads.threadid()][2]>0 && otherCords[Threads.threadid()][2]<=n
                        mergeClusters!(tmpGraphs[Threads.threadid()],j,index(otherCords[Threads.threadid()],n),tmpNews[Threads.threadid()],cords)
                    end
                end
            end
        end
        lengths[i]=corrLength(tmpGraphs[Threads.threadid()],tmpNews[Threads.threadid()])
    end
end


println(Threads.nthreads())

const pointfile=string(ARGS[1])
const n=parse(Int,ARGS[2])
const k=parse(Int,ARGS[3])

const cords=[[mod(i-1,n)+1,div(i-1,n)+1] for i=1:n^2];#Define the coordinates of each site
tmpGraphs=[initialiseDisjointSetGraph(n,cords) for i=1:Threads.nthreads()];#Initialise the graphs in a threadsafe array to pr-allocate memory
tmpNews=[[0.,0.] for _=1:Threads.nthreads()];#Initialise array used to pre-allocate memory
edges=[[1,0],[0,1],[1,1],[1,-1]];#Define strain 1 and 2 edges respectively, being careful not to overcount
otherCords=[[0,0] for _=1:Threads.nthreads()];#Initialise array used to pre-allocate memory

const points=readdlm(pointfile,',',Float64,'\n');

lengths=[zeros(Float64,k) for _=1:div(length(points),2)];#Initialise array used to stor outputs

open("out0_"*string(n)*"_"*string(k)*"_"*pointfile,"w") do io#Create file if it doesn't exist, otherwise overrite existing dfile of same name.
    nothing
end
open("out0_"*string(n)*"_"*string(k)*"_"*pointfile,"a") do io #Open file for appending
    t=@time @showprogress 1 "Computing..." for i=eachindex(lengths)
        findCorrLength!(lengths[i],n,points[i,1],points[i,2],tmpGraphs,tmpNews,edges,otherCords,cords)
        writedlm(io,sum(lengths[i])/k,',')
    end
    println(t)
end



exit()


n=30
k=16

pointfile="line_points.csv"

@profview findCorrLength!(lengths[1],n,0.1,0.1,tmpGraphs,tmpNews,edges,currCords,otherCords,ps,cords)

@time @profview for i=eachindex(lengths)
    findCorrLength!(lengths[i],n,points[i,1],points[i,2],tmpGraphs,tmpNews,edges,currCords,otherCords,ps,cords)
end

@time @profview_allocs for i=eachindex(lengths)
    findCorrLength!(lengths[i],n,points[i,1],points[i,2],tmpGraphs,tmpNews,edges,currCords,otherCords,ps,cords)
end

