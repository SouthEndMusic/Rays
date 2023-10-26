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