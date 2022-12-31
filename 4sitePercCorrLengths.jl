using DelimitedFiles: readdlm, writedlm
using ProgressMeter
include("disjointSetGraphStruct.jl")


function findCorrLength!(lengths::Vector{Float64},n::Int,β::Float64,μ::Float64,tmpGraphs::Vector{disjointSetGraph},
    tmpNews::Vector{Vector{Float64}},edges::Vector{Vector{Vector{Int}}},otherCords::Vector{Vector{Int}},ps::Vector{Float64},cords,vertexDegrees::Vector{Vector{Int}},closedEdges::Vector{Vector{Vector{Bool}}})
    ps[1]=1-1/(1+exp(β*(μ)))
    ps[2]=1-1/(1+exp(β*(μ-1)))
    Threads.@threads for i=eachindex(lengths)
        #Reset everything
        resetDisjointSetGraph!(tmpGraphs[Threads.threadid()],cords)
        for j=1:n^2
            vertexDegrees[Threads.threadid()][j]=0
        end
        for j=1:n^2
            for k=1:4
                closedEdges[Threads.threadid()][j][k]=false
            end
        end
        otherCords[Threads.threadid()][1]=0
        otherCords[Threads.threadid()][2]=0
        #Find the edges
        for j=1:n^2
            for i=eachindex(edges[1])
                if ps[1]>rand()
                    otherCords[Threads.threadid()][1]=cords[j][1]+edges[1][i][1]
                    otherCords[Threads.threadid()][2]=cords[j][2]+edges[1][i][2]
                    if otherCords[Threads.threadid()][1]>0 && otherCords[Threads.threadid()][1]<=n && otherCords[Threads.threadid()][2]>0 && otherCords[Threads.threadid()][2]<=n
                        vertexDegrees[Threads.threadid()][j]+=1
                        vertexDegrees[Threads.threadid()][index(otherCords[Threads.threadid()],n)]+=1
                        closedEdges[Threads.threadid()][j][i]=true
                    end
                end
            end

            for i=eachindex(edges[2])
                if ps[2]>rand()
                    otherCords[Threads.threadid()][1]=cords[j][1]+edges[2][i][1]
                    otherCords[Threads.threadid()][2]=cords[j][2]+edges[2][i][2]
                    if otherCords[Threads.threadid()][1]>0 && otherCords[Threads.threadid()][1]<=n && otherCords[Threads.threadid()][2]>0 && otherCords[Threads.threadid()][2]<=n
                        vertexDegrees[Threads.threadid()][j]+=1
                        vertexDegrees[Threads.threadid()][index(otherCords[Threads.threadid()],n)]+=1
                        closedEdges[Threads.threadid()][j][i+2]=true
                    end
                end
            end
        end
        #Find the clusters
        for j=1:n^2
            if vertexDegrees[Threads.threadid()][j]==4

                for i=eachindex(edges[1])
                    if closedEdges[Threads.threadid()][j][i]
                        otherCords[Threads.threadid()][1]=cords[j][1]+edges[1][i][1]
                        otherCords[Threads.threadid()][2]=cords[j][2]+edges[1][i][2]
                        if otherCords[Threads.threadid()][1]>0 && otherCords[Threads.threadid()][1]<=n && otherCords[Threads.threadid()][2]>0 && otherCords[Threads.threadid()][2]<=n
                            if vertexDegrees[Threads.threadid()][index(otherCords[Threads.threadid()],n)]==4
                                mergeClusters!(tmpGraphs[Threads.threadid()],j,index(otherCords[Threads.threadid()],n),tmpNews[Threads.threadid()],cords)
                            end
                        end
                    end
                end

                for i=eachindex(edges[2])
                    if closedEdges[Threads.threadid()][j][i+2]
                        otherCords[Threads.threadid()][1]=cords[j][1]+edges[2][i][1]
                        otherCords[Threads.threadid()][2]=cords[j][2]+edges[2][i][2]
                        if otherCords[Threads.threadid()][1]>0 && otherCords[Threads.threadid()][1]<=n && otherCords[Threads.threadid()][2]>0 && otherCords[Threads.threadid()][2]<=n
                            if vertexDegrees[Threads.threadid()][index(otherCords[Threads.threadid()],n)]==4
                                mergeClusters!(tmpGraphs[Threads.threadid()],j,index(otherCords[Threads.threadid()],n),tmpNews[Threads.threadid()],cords)
                            end
                        end
                    end
                end

            end
        end
        #Find the length
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
edges=[[[1,0],[0,1]],[[1,1],[1,-1]]];#Define strain 1 and 2 edges respectively, being careful not to overcount
otherCords=[[0,0] for _=1:Threads.nthreads()];#Initialise array used to pre-allocate memory
ps=[0.,0.];#Initialise array used to pre-allocate memory
vertexDegrees=[[0 for _=1:n^2] for i=1:Threads.nthreads()];#Initialise array used to pre-allocate memory
closedEdges=fill(fill(fill(false,4),n^2),Threads.nthreads());#Initialise array used to pre-allocate memory

const points=readdlm(pointfile,',',Float64,'\n');

lengths=[zeros(Float64,k) for _=1:div(length(points),2)];#Initialise array used to stor outputs

open("out4_"*string(n)*"_"*string(k)*"_"*pointfile,"w") do io#Create file if it doesn't exist, otherwise overrite existing dfile of same name.
    nothing
end
open("out4_"*string(n)*"_"*string(k)*"_"*pointfile,"a") do io #Open file for appending
    t=@time @showprogress 1 "Computing..." for i=eachindex(lengths)
        findCorrLength!(lengths[i],n,points[i,1],points[i,2],tmpGraphs,tmpNews,edges,otherCords,ps,cords,vertexDegrees,closedEdges)
        writedlm(io,sum(lengths[i])/k,',')
    end
    println(t)
end



exit()


n=30
k=16

pointfile="points.csv"

@profview findCorrLength!(lengths[1],n,points[2,1],points[2,2],tmpGraphs,tmpNews,edges,otherCords,ps,cords,vertexDegrees,closedEdges)

tmpGraphs[1]
corrLength(tmpGraphs[1],tmpNews[1])

@time @profview for i=eachindex(lengths)
    findCorrLength!(lengths[i],n,points[i,1],points[i,2],tmpGraphs,tmpNews,edges,currCords,otherCords,ps,cords)
end

@time @profview_allocs for i=eachindex(lengths)
    findCorrLength!(lengths[i],n,points[i,1],points[i,2],tmpGraphs,tmpNews,edges,currCords,otherCords,ps,cords)
end

