abstract type Parameters{F<:AbstractFloat} end

mutable struct Interactor{F<:AbstractFloat}
    window::Ptr{SDL_Window}
    renderer::Ptr{SDL_Renderer}
    scene::Scene{F}
    parameters::Parameters{F}
    get_render::Function
    affect_parameters_input!::Function
    affect_parameters_time!::Function
    dtimer::Dtimer
end

function Interactor(
    scene::Scene{F},
    parameters::Parameters{F},
    affect_parameters_input!::Function,
    affect_parameters_time!::Function,
    get_render::Function;
    window_name::String = "Interactive rendering",
)::Interactor{F} where {F}
    window = SDL.SDL_CreateWindow(
        window_name,
        SDL.SDL_WINDOWPOS_CENTERED,
        SDL.SDL_WINDOWPOS_CENTERED,
        scene.cameras[1].screen_res...,
        SDL.SDL_WINDOW_SHOWN,
    )
    SDL.SDL_SetWindowResizable(window, false)

    renderer = SDL.SDL_CreateRenderer(
        window,
        -1,
        SDL.SDL_RENDERER_ACCELERATED | SDL.SDL_RENDERER_PRESENTVSYNC,
    )

    dtimer = Dtimer(convert(F, 0))

    return Interactor(
        window,
        renderer,
        scene,
        parameters,
        get_render,
        affect_parameters_input!,
        affect_parameters_time!,
        dtimer,
    )
end

function Base.convert(::Type{SDL.LibSDL2.SDL_bool}, bool::Bool)
    return bool ? SDL.SDL_TRUE : SDL.SDL_FALSE
end

function set_render!(interactor::Interactor)::Nothing
    (; renderer, scene, parameters, get_render) = interactor
    render = convert.(UInt8, round.(255 .* get_render(scene, parameters)))
    render = permutedims(render, [1, 3, 2])
    depth = 24
    pitch = sizeof(eltype(render)) * 3 * size(render)[2]
    surface = SDL.SDL_CreateRGBSurfaceFrom(
        render,
        size(render)[2:3]...,
        depth,
        pitch,
        0xFF0000,
        0x00FF00,
        0x0000FF,
        0,
    )
    texture = SDL.SDL_CreateTextureFromSurface(renderer, surface)
    SDL.SDL_RenderCopy(renderer, texture, SDL.C_NULL, SDL.C_NULL)
    SDL.SDL_DestroyTexture(texture)
    SDL.SDL_RenderPresent(renderer)
    return nothing
end

function run!(interactor::Interactor)::Nothing
    (;
        dtimer,
        parameters,
        affect_parameters_input!,
        affect_parameters_time!,
        renderer,
        window,
    ) = interactor
    close = false
    try
        while !close
            changed = false
            Δt = get_Δt!(dtimer)
            changed |= affect_parameters_time!(parameters, Δt)
            event_ref = Ref{SDL.SDL_Event}()
            while Bool(SDL.SDL_PollEvent(event_ref))
                event = event_ref[]
                if event.type == SDL.SDL_QUIT
                    close = true
                    break
                else
                    changed |= affect_parameters_input!(parameters, event, Δt)
                end
            end
            if changed
                set_render!(interactor)
            end
        end
    finally
        SDL.SDL_DestroyRenderer(renderer)
        SDL.SDL_DestroyWindow(window)
        SDL.SDL_Quit()
    end
    return nothing
end
