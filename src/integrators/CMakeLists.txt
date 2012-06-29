include_directories(${YAF_INCLUDE_DIRS})

add_definitions(-DBUILDING_YAFRAYPLUGIN)

add_library(directlight SHARED directlight.cc)
target_link_libraries(directlight yafaraycore)

add_library(photonmap SHARED photonintegr.cc)
target_link_libraries(photonmap yafaraycore)

add_library(pathtrace SHARED pathtracer.cc)
target_link_libraries(pathtrace yafaraycore)

add_library(bidirpath SHARED bidirpath.cc)
target_link_libraries(bidirpath yafaraycore)

add_library(EmissionIntegrator SHARED EmissionIntegrator.cc)
target_link_libraries(EmissionIntegrator yafaraycore)

add_library(SingleScatterIntegrator SHARED SingleScatterIntegrator.cc)
target_link_libraries(SingleScatterIntegrator yafaraycore)

add_library(EmptyVolumeIntegrator SHARED EmptyVolumeIntegrator.cc)
target_link_libraries(EmptyVolumeIntegrator yafaraycore)

add_library(DebugIntegrator SHARED DebugIntegrator.cc)
target_link_libraries(DebugIntegrator yafaraycore)

add_library(SkyIntegrator SHARED SkyIntegrator.cc)
target_link_libraries(SkyIntegrator yafaraycore)

add_library(sppm SHARED sppm.cc)
target_link_libraries(sppm yafaraycore)

add_library(pbgi_integrator SHARED Pbgi_integrator.cc Pbgi_integrator/Cube_raster_buffer.cc Pbgi_integrator/Scene_sampler.cpp Pbgi_integrator/Spherical_function.cc)
target_link_libraries(pbgi_integrator yafaraycore)

install (TARGETS directlight photonmap pathtrace bidirpath sppm pbgi_integrator
		EmissionIntegrator SingleScatterIntegrator EmptyVolumeIntegrator DebugIntegrator SkyIntegrator
		${YAF_TARGET_TYPE} DESTINATION ${YAF_PLUGIN_DIR})
