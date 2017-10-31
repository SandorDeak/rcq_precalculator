#include "terrain_generator.h"

#include <algorithm>
#include <random>

# define PI           3.14159265358979323846


terrain_generator::terrain_generator(const glm::uvec2& size) : m_terrain(size)
{
}


terrain_generator::~terrain_generator()
{
}

void terrain_generator::init_terrain()
{
	double sand_init = 0.2;
	double nutrition_init = 0.01;
	double water_init = sand_init*sand_water_capacity;

	auto center = m_terrain.center();
	auto size = m_terrain.size();

	for (uint32_t i = 0; i < size.x; ++i)
	{
		for (uint32_t j = 0; j < size.y; ++j)
		{
			glm::dvec2 pos = { static_cast<double>(i) - center.x, static_cast<double>(j) - center.y };
			double curvature = m_boundary_amplitude*sin(m_boundary_frequency*pos.y);
			pos.x += curvature;

			auto& c = m_terrain[i][j];
			c.rock_height = exp(-0.5*glm::dot(pos, m_tectonic_planes_covariance*pos))*m_tectonic_planes_normalisation_factor;
			c.sand = sand_init;
			c.nutrition = nutrition_init;
			c.water = water_init;
		}
	}
}

void terrain_generator::elevate_terrain()
{
	auto center = m_terrain.center();
	auto size = m_terrain.size();

	for (uint32_t i = 0; i < size.x; ++i)
	{
		for (uint32_t j = 0; j < size.y; ++j)
		{
			glm::dvec2 pos = { static_cast<double>(i) - center.x, static_cast<double>(j) - center.y };
			double curvature = m_boundary_amplitude*sin(m_boundary_frequency*pos.y);
			pos.x += curvature;

			auto& c = m_terrain[i][j];
			double elevation= exp(-0.5*glm::dot(pos, m_tectonic_planes_covariance*pos))*
				m_tectonic_planes_normalisation_factor*m_elevation_speed;
			c.rock_height += elevation;
		}
	}
}

void terrain_generator::coast()
{
	double coast_destroy_scale_sand = 0.6;
	double coast_destroy_scale_rock = 0.00001;
	double coast_destroy_nutrition = 0.7;

	glm::uvec2 size = m_terrain.size();
	for (uint32_t i = 0; i < size.x; ++i)
	{
		for (uint32_t j = 0; j < size.y; ++j)
		{
			auto& c= m_terrain[i][j];
			if (c.sand > 0.)
				c.sand *= coast_destroy_scale_sand;
			else
				c.rock_height -= coast_destroy_scale_rock;

			c.nutrition *= coast_destroy_nutrition;
			c.water = c.sand*sand_water_capacity;
		}
	}
}

void terrain_generator::nutrition()
{
	double nutrition_grow_factor = 1.01;

	glm::uvec2 size = m_terrain.size();
	for (uint32_t i = 0; i < size.x; ++i)
	{
		for (uint32_t j = 0; j < size.y; ++j)
		{
			m_terrain[i][j].nutrition *= nutrition_grow_factor;
		}
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<uint32_t> dis_x(0, size.x-1);
	std::uniform_int_distribution<uint32_t> dis_y(0, size.y-1);

	uint32_t random_events_count = 3;
	double random_nutrition = 0.05;
	for (uint32_t i = 0; i < random_events_count; ++i)
	{
		uint32_t x = dis_x(gen);
		uint32_t y = dis_y(gen);
		m_terrain[x][y].nutrition += random_nutrition;
	}
}

void terrain_generator::rain(double amount)
{
	glm::uvec2 size = m_terrain.size();

	for (uint32_t i = 0; i < size.x; ++i)
	{
		for (uint32_t j = 0; j < size.y; ++j)
		{
			auto& c = m_terrain[i][j];

			c.aggr.clear();
			c.aggr.water = amount;
			c.aggr.velocity = glm::dvec3(0., 0., -10.);
			c.aggr.exists = true;
			c.next_aggr.clear();
			c.next_aggr.exists=false;
		}
	}

	bool need_iteration = true;

	while (need_iteration)
	{

		m_terrain.calc_heights();
		m_terrain.calc_normals();

		m_terrain.collide_with_aggreagates();
		m_terrain.deposit_from_aggregates();

		m_terrain.calc_heights();
		m_terrain.update_velocities();
		m_terrain.distribute_aggregates();
		need_iteration = m_terrain.set_next_aggregate_to_aggregate();
	}
}

void terrain_generator::wind(const glm::dvec3 & dir, double velocity)
{

}