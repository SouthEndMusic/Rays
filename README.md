# Rays
CPU based 3D graphics in Julia using raycasting.

![texture_examples](https://github.com/SouthEndMusic/Rays/assets/74617371/45b666bf-fd20-4a7f-8943-5fd3aa71a4df)

## Main features

### Shapes
- Sphere, cube, tetrahedron
- Implicit surface
- Surface of revolution
- Triangle mesh
- Fractal shapes, e.g. [Menger sponge](https://nl.wikipedia.org/wiki/Spons_van_Menger), [Sierpinski pyramid](https://en.wikipedia.org/wiki/Sierpi%C5%84ski_triangle#Analogues_in_higher_dimensions)

For examples see [here](https://github.com/SouthEndMusic/Rays/blob/master/examples/Getting_started_shapes.ipynb).

### Textures
- Uniform color
- Color by face
- Color by shape coordinates

For examples see [here](https://github.com/SouthEndMusic/Rays/blob/master/examples/Getting_started_textures.ipynb).

### Interactive rendering

Integraton of [SimpleDirectMediaLayer.jl](https://github.com/JuliaMultimedia/SimpleDirectMediaLayer.jl) provides a simple API to use `Rays.jl` essentially as a game engine. For an example script see [here](https://github.com/SouthEndMusic/Rays/blob/master/examples/Getting_started_interactive.ipynb).

## Contributing

This is a hobby project of mine, partly to learn how to write good Julia code. Therefore I do not accept contributions, only suggestions 🙂 However, please feel free to experiment with the code and show me your creations (kudos to whomever creates a clone of DOOM).
