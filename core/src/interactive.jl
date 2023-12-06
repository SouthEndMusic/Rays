abstract type Parameters{F<:AbstractFloat} end

const GetRender = FunctionWrapper{Array{F,3},Tuple{Scene{F},P}} where {F,P}
const AffectParametersInput = FunctionWrapper{Bool,Tuple{P,SDL_Event,F}} where {F,P}
const AffectParametersTime = FunctionWrapper{Bool,Tuple{P,F}} where {F,P}

"""
Object with all the data for interactive rendering.
"""
struct Interactor{F<:AbstractFloat,P<:Parameters{F}}
    window::Ptr{SDL_Window}
    renderer::Ptr{SDL_Renderer}
    render_UInt8::Array{UInt8,3}
    scene::Scene{F}
    parameters::P
    get_render::GetRender{F,P}
    affect_parameters_input!::AffectParametersInput{F,P}
    affect_parameters_time!::AffectParametersTime{F,P}
    dtimer::Dtimer
    camera_main::Camera{F}
end

"""
Construct an interactor.
Assumes the window resolution is equal to the resolution of the first camera.
"""
function Interactor(
    scene::Scene{F},
    parameters::P,
    affect_parameters_input!::Function,
    affect_parameters_time!::Function,
    get_render::Function;
    window_name::String = "Interactive rendering",
    name_main_camera::Union{Symbol,Nothing} = nothing,
)::Interactor{F,P} where {F,P<:Parameters{F}}
    if isnothing(name_main_camera)
        camera_main = first(values(scene.cameras))
    else
        camera_main = scene.cameras[name_main_camera]
    end
    window = SDL.SDL_CreateWindow(
        window_name,
        SDL.SDL_WINDOWPOS_CENTERED,
        SDL.SDL_WINDOWPOS_CENTERED,
        get_screen_res(camera_main)...,
        SDL.SDL_WINDOW_SHOWN,
    )
    SDL.SDL_SetWindowResizable(window, false)

    renderer = SDL.SDL_CreateRenderer(
        window,
        -1,
        SDL.SDL_RENDERER_ACCELERATED | SDL.SDL_RENDERER_PRESENTVSYNC,
    )

    render_UINT8 = zeros(UInt8, size(camera_main.canvas)...)

    dtimer = Dtimer(convert(F, 0))

    return Interactor(
        window,
        renderer,
        render_UINT8,
        scene,
        parameters,
        # Not sure why these conversions do not happen automatically
        GetRender{F,P}(get_render),
        AffectParametersInput{F,P}(affect_parameters_input!),
        AffectParametersTime{F,P}(affect_parameters_time!),
        dtimer,
        camera_main,
    )
end

"""
Conver between SDL and julia booleans.
"""
function Base.convert(::Type{SDL.LibSDL2.SDL_bool}, bool::Bool)
    return bool ? SDL.SDL_TRUE : SDL.SDL_FALSE
end

"""
Compute a render and put it to the SDL window.
"""
function set_render!(interactor::Interactor)::Nothing
    (; renderer, scene, parameters, get_render, render_UInt8, camera_main) = interactor
    canvas = camera_main.canvas
    copyto!(canvas, get_render(scene, parameters))
    render_UInt8 .= convert.(UInt8, round.(255 .* canvas))
    render = permutedims(render_UInt8, [1, 3, 2])
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

"""
The main loop of interactive rendering.
"""
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
            parameters_changed = false
            Δt = get_Δt!(dtimer)
            parameters_changed |= affect_parameters_time!(parameters, Δt)
            event_ref = Ref{SDL.SDL_Event}()
            while Bool(SDL.SDL_PollEvent(event_ref))
                event = event_ref[]
                if event.type == SDL.SDL_QUIT
                    close = true
                    break
                else
                    parameters_changed |= affect_parameters_input!(parameters, event, Δt)
                end
            end
            if parameters_changed
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
