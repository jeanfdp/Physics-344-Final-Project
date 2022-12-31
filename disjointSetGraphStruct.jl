function euclideanDistSq(x, y)
    return (x[1]-y[1])^2+(x[2]-y[2])^2
end

function index(cord,n::Int)
    return (cord[2]-1)*n+cord[1]
end

struct disjointSetGraph
    length::Int
    parents::Vector{Int}
    centers::Vector{Vector{Float64}}
    touching::Vector{Vector{Bool}}
    masses::Vector{Float64}
    inertias::Vector{Float64}
end

function initialiseDisjointSetGraph(n::Int,cords)
    return disjointSetGraph(n,collect(1:n^2),[[cords[i][1],cords[i][2]] for i=1:n^2],[[cords[i][1]==1,cords[i][1]==n,cords[i][2]==1,cords[i][2]==n] for i=1:n^2],fill(1/n^2,n^2),zeros(Float64,n^2))
end

function resetDisjointSetGraph!(graph::disjointSetGraph,cords)
    for i=eachindex(graph.parents)
        graph.parents[i]=i
        graph.centers[i].=cords[i]
        graph.touching[i][1]=cords[i][1]==1
        graph.touching[i][2]=cords[i][1]==graph.length
        graph.touching[i][3]=cords[i][2]==1
        graph.touching[i][4]=cords[i][2]==graph.length        
        graph.masses[i]=1/graph.length^2
        graph.inertias[i]=0.
    end
end


function root(graph::disjointSetGraph,i)
    if graph.parents[i]==i
        return i
    else
        graph.parents[i]=root(graph,graph.parents[i])
        return graph.parents[i]
    end
end

function mergeClusters!(graph::disjointSetGraph,v1::Int,v2::Int,tmpNews,cords)
    v1=root(graph,v1)
    v2=root(graph,v2)
    if v1!=v2
        tmpNews.=graph.centers[v1].*graph.masses[v1]

        graph.masses[v1]+=graph.masses[v2]

        tmpNews.+=(graph.centers[v2].*graph.masses[v2])
        tmpNews./=graph.masses[v1]

        graph.parents[v2]=v1
        graph.inertias[v1]=graph.inertias[v1]+graph.inertias[v2]+(graph.masses[v1]-graph.masses[v2])*euclideanDistSq(graph.centers[v1],tmpNews)+graph.masses[v2]*euclideanDistSq(graph.centers[v2],tmpNews)
        graph.centers[v1].=tmpNews

        graph.touching[v1][1]=graph.touching[v1][1]||graph.touching[v2][1]
        graph.touching[v1][2]=graph.touching[v1][2]||graph.touching[v2][2]
        graph.touching[v1][3]=graph.touching[v1][3]||graph.touching[v2][3]
        graph.touching[v1][4]=graph.touching[v1][4]||graph.touching[v2][4]
    end
end
    
function corrLength(graph::disjointSetGraph,tmpNews)
    fill!(tmpNews,0.)
    for i in eachindex(graph.parents)
        if graph.parents[i]==i
            if !((graph.touching[i][1]&&graph.touching[i][2])||(graph.touching[i][3]&&graph.touching[i][4]))#Only include non-spanning components
                tmpNews[1]+=graph.masses[i]*graph.inertias[i]
                tmpNews[2]+=graph.masses[i]^2
            end
        end
    end
    if tmpNews[2]==0.
        return 0.
    else
        return sqrt(tmpNews[1]/tmpNews[2])
    end
end