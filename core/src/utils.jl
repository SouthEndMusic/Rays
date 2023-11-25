mutable struct Dtimer{F <: AbstractFloat}
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
	name_snake_case = join([part.match for part in parts], "_") |> lowercase
	return Symbol(name_snake_case)
end

function identity_matrix(F)::Matrix
	return F[1 0 0; 0 1 0; 0 0 1]
end
