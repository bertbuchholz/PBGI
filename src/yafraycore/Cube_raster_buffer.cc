#include <map>

#include <utilities/CubeRasterBuffer.h>

__BEGIN_YAFRAY

std::map<std::string, Cube_raster_buffer::Splat_type> create_splat_type_map()
{
    std::map<std::string, Cube_raster_buffer::Splat_type> m;
    m["Single_pixel"] = Cube_raster_buffer::Single_pixel;
    m["Disc_tracing"] = Cube_raster_buffer::Disc_tracing;
    m["AA_square"] = Cube_raster_buffer::AA_square;
    return m;
}

std::map<std::string, Cube_raster_buffer::Splat_type> Cube_raster_buffer::enum_splat_type_map = create_splat_type_map();



std::map<std::string, Cube_raster_buffer::Type> create_type_map()
{
    std::map<std::string, Cube_raster_buffer::Type> m;
    m["Simple"] = Cube_raster_buffer::Simple;
    m["Accumulating"] = Cube_raster_buffer::Accumulating;
    return m;
}

std::map<std::string, Cube_raster_buffer::Type> Cube_raster_buffer::enum_type_map = create_type_map();

__END_YAFRAY
