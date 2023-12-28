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
