const TypedSubArray =
    SubArray{T,1,MT,Tuple{Int64,ST},true} where {T,MT<:AbstractMatrix{T},ST}

# Wrapper around intersection functionality of a shape, for type stability
const Intersector = FunctionWrapper{
    Nothing,
    Tuple{
        VFS,
        VFS,
        VFS,
        VFS,
        VFS,
        VIS,
        VFS,
        TypedSubArray{Symbol,Matrix{Symbol},Base.Slice{Base.OneTo{Int64}}},
    },
} where {
    VFS<:TypedSubArray{
        F,
        MF,
        Base.Slice{Base.OneTo{Int64}},
    } where {F<:AbstractFloat,MF<:AbstractMatrix{F}},
    VIS<:TypedSubArray{
        Int,
        MI,
        Base.Slice{Base.OneTo{Int64}},
    } where {MI<:AbstractMatrix{Int}},
}

# Wrapper around texturing functionality of a shape, for type stability
const Texturer = FunctionWrapper{
    Nothing,
    Tuple{VFS,VIS,VFS,VFS,VFS,VFS},
} where {
    VFS<:TypedSubArray{
        F,
        MF,
        Base.Slice{Base.OneTo{Int64}},
    } where {F<:AbstractFloat,MF<:AbstractMatrix{F}},
    VIS<:TypedSubArray{
        Int,
        MI,
        Base.Slice{Base.OneTo{Int64}},
    } where {MI<:AbstractMatrix{Int}},
}

# Type stable function wrapper: R ↦ R
const ScalarFunc = FunctionWrapper{F,Tuple{F}} where {F<:AbstractFloat}

# Type stable function wrapper: R^n ↦ R^n (in-place)
const Transform = FunctionWrapper{
    Nothing,
    Tuple{TypedSubArray{F,MF,ST}},
} where {F<:AbstractFloat,MF<:AbstractMatrix{F},ST}

# Type stable function wrapper: R^n ↦ R
const ScalarField = FunctionWrapper{
    F,
    Tuple{TypedSubArray{F,MF,ST}},
} where {ST,F<:AbstractFloat,MF<:AbstractMatrix{F}}

# Type stable function wrapper: R^n ↦ R^n (in-place on first function argument)
const VectorField = FunctionWrapper{
    Nothing,
    Tuple{
        TypedSubArray{F,MF,Base.Slice{Base.OneTo{Int64}}},
        TypedSubArray{F,MF,UnitRange{Int64}},
    },
} where {F<:AbstractFloat,MF<:AbstractMatrix{F}}

mutable struct Dtimer{F<:AbstractFloat}
    t::F
end

function Base.convert(::Type{Dtimer{F}}, dtimer::Dtimer) where {F}
    return Dtimer(convert(F, dtimer.t))
end

function start!(dtimer::Dtimer)::Nothing
    dtimer.t = time_ns() / 1e9
    return nothing
end

function get_Δt!(dtimer::Dtimer{F})::F where {F}
    t = time_ns() / 1e9
    Δt = t - dtimer.t
    dtimer.t = t
    return Δt
end

function snake_case_name(T)::Symbol
    name = String(Base.typename(T).name)
    parts = eachmatch(r"[A-Z][a-z]*", name)
    name_snake_case = join([part.match for part ∈ parts], "_") |> lowercase
    return Symbol(name_snake_case)
end

function identity_matrix(F)::Matrix
    return F[1 0 0; 0 1 0; 0 0 1]
end

struct BoundingBox{F<:AbstractFloat}
    coordinates_min::Vector{F}
    coordinates_max::Vector{F}
end

"""
Get the center of a bounding box.
"""
function center(bounding_box::BoundingBox{F})::Vector{F} where {F}
    return (bounding_box.coordinates_min + bounding_box.coordinates_max) / 2
end

"""
Determine whether a ray intersects a bounding box.
"""
function bounding_box_intersect(
    ray_loc::AbstractVector{F},
    ray_dir::AbstractVector{F},
    bounding_box::BoundingBox{F},
)::Bool where {F}
    (; coordinates_min, coordinates_max) = bounding_box
    t_min = zero(F)
    t_max = F(Inf)
    for dim ∈ 1:3
        t_1 = (coordinates_min[dim] - ray_loc[dim]) / ray_dir[dim]
        t_2 = (coordinates_max[dim] - ray_loc[dim]) / ray_dir[dim]
        t_min_candidate = min(t_1, t_2)
        t_max_candidate = max(t_1, t_2)
        t_min = max(t_min, t_min_candidate)
        t_max = min(t_max, t_max_candidate)
        if t_min > t_max
            return false
        end
    end
    return true
end

"""
Node in the scene partition graph.
bounding_box: The bounding box corresponding to this node
children: The indices of nodes with bounding boxes contained in the bounding
	box of this node
shape_names: The names of the shapes whose bounding box intersects the bounding
	box of this node
depth: the depth in the graph of this node
dim_split: The dimension along which the bounding box of this node was split to
	form child nodes (if this node is a leaf-node then dim_split = 0)
"""
struct PartitionNode{F<:AbstractFloat}
    bounding_box::BoundingBox{F}
    children::Vector{Int}
    shape_names::Vector{Symbol}
    depth::Int
    dim_split::Int
end
