//#include "terrain_generator.h"
//
//#include <algorithm>
//
//# define PI           3.14159265358979323846
//
//
//terrain_generator::terrain_generator(const glm::uvec2& size) : m_terrain(size), m_rd(), m_rgen(m_rd()),
//	m_window(sf::VideoMode(size.x, size.y), "Debug info")
//{
//}
//
//
//terrain_generator::~terrain_generator()
//{
//}
//
//void terrain_generator::init_terrain()
//{
//	double sand_init = 0.2;
//	double nutrition_init = 0.01;
//	double sea_bottom = -100.;
//	double water_init = sand_init*sand_water_capacity;
//
//	auto center = m_terrain.center();
//	auto size = m_terrain.size();
//
//	m_boundary_amplitude = 10.;
//	m_boundary_frequency = 0.07;
//
//	double deviation_x = 80.;
//	double deviation_y = 100.;
//
//	m_elevation_speed = 1e4;
//
//
//	m_tectonic_planes_covariance = { 1./pow(deviation_x, 2.), 0., 0., 1./pow(deviation_y, 2.) };
//	m_tectonic_planes_normalisation_factor = 1. / (2.*PI*deviation_x*deviation_y);
//
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			glm::dvec2 pos = { static_cast<double>(i) - center.x, static_cast<double>(j) - center.y };
//			double curvature = m_boundary_amplitude*sin(m_boundary_frequency*pos.y);
//			pos.x += curvature;
//
//			auto& c = m_terrain[i][j];
//			double k = exp(-0.5*glm::dot(pos, m_tectonic_planes_covariance*pos));
//			c.rock_height = sea_bottom + exp(-0.5*glm::dot(pos, m_tectonic_planes_covariance*pos))*m_tectonic_planes_normalisation_factor;
//			c.sand = sand_init;
//			c.nutrition = nutrition_init;
//			c.water = water_init;
//		}
//	}
//}
//
//void terrain_generator::elevate_terrain()
//{
//	auto center = static_cast<glm::dvec2>(m_terrain.center());
//	auto size = m_terrain.size();
//
//	std::normal_distribution<> dis_ampl(0., m_boundary_amplitude*0.05);
//	std::normal_distribution<> dis_freq(0., m_boundary_frequency*0.05);
//	/*m_boundary_amplitude += dis_ampl(m_rgen);
//	m_boundary_frequency += dis_freq(m_rgen);*/
//
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			glm::dvec2 pos = { static_cast<double>(i) - center.x, static_cast<double>(j) - center.y };
//			double curvature = m_boundary_amplitude*sin(m_boundary_frequency*pos.y);
//			pos.x += curvature;
//
//			auto& c = m_terrain[i][j];
//			double elevation= exp(-0.5*glm::dot(pos, m_tectonic_planes_covariance*pos))*
//				m_tectonic_planes_normalisation_factor*m_elevation_speed;
//			c.rock_height += elevation;
//		}
//	}
//}
//
//void terrain_generator::coast()
//{
//	double coast_destroy_scale_sand = 0.6;
//	double coast_destroy_scale_rock = 0.00001;
//	double coast_destroy_nutrition = 0.7;
//
//	glm::uvec2 size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			auto& c= m_terrain[i][j];
//			if (c.sand > 0.)
//				c.sand *= coast_destroy_scale_sand;
//			else
//				c.rock_height -= coast_destroy_scale_rock;
//
//			c.nutrition *= coast_destroy_nutrition;
//			c.water = c.sand*sand_water_capacity;
//		}
//	}
//}
//
//void terrain_generator::nutrition()
//{
//	double nutrition_grow_factor = 1.01;
//
//	glm::uvec2 size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			m_terrain[i][j].nutrition *= nutrition_grow_factor;
//		}
//	}
//
//	std::uniform_int_distribution<uint32_t> dis_x(0, size.x-1);
//	std::uniform_int_distribution<uint32_t> dis_y(0, size.y-1);
//
//	uint32_t random_events_count = 3;
//	double random_nutrition = 0.05;
//	for (uint32_t i = 0; i < random_events_count; ++i)
//	{
//		uint32_t x = dis_x(m_rgen);
//		uint32_t y = dis_y(m_rgen);
//		m_terrain[x][y].nutrition += random_nutrition;
//	}
//}
//
//void terrain_generator::rain(double amount)
//{
//	glm::uvec2 size = m_terrain.size();
//
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			auto& c = m_terrain[i][j];
//
//			c.aggr.clear();
//			c.aggr.water = amount;
//			c.aggr.velocity = glm::dvec3(0., 0., -10.);
//			c.aggr.exists = true;
//			c.next_aggr.clear();
//			c.next_aggr.exists=false;
//		}
//	}
//
//	bool need_iteration = true;
//
//	while (need_iteration)
//	{
//
//		m_terrain.calc_heights();
//		m_terrain.calc_normals();
//
//		m_terrain.collide_with_aggreagates();
//		m_terrain.deposit_from_aggregates();
//
//		m_terrain.calc_heights();
//		m_terrain.update_velocities();
//		m_terrain.distribute_aggregates();
//		need_iteration = m_terrain.set_next_aggregate_to_aggregate();
//	}
//}
//
//void terrain_generator::wind(const glm::dvec3& wind, uint32_t iteration)
//{
//	constexpr double coeff0 = 10.;
//	constexpr double coeff1 = 0.001;
//
//	glm::uvec2 size = m_terrain.size();
//
//	m_terrain.calc_heights();
//	m_terrain.calc_normals();
//	std::uniform_real_distribution<> travelling_distance(0., glm::length(wind)*coeff0);
//	std::uniform_real_distribution<> bias(-glm::length(wind)*coeff1, glm::length(wind)*coeff1);
//
//	for (uint32_t iter = 0; iter < iteration; ++iter)
//	{
//		m_terrain.calc_heights();
//		m_terrain.calc_normals();
//
//		for (uint32_t i = 0; i < size.x; ++i)
//		{
//			for (uint32_t j = 0; j < size.y; ++j)
//			{
//				auto& c = m_terrain[i][j];
//				if (c.sand == 0.)
//				{
//					double eroded = std::max(glm::dot(c.normal, -wind), 0.0001)*wind_erosion_rock;
//					c.rock_height -= eroded;
//					glm::dvec3 travelling_dist = glm::normalize(wind)*travelling_distance(m_rgen) +
//						glm::dvec3(bias(m_rgen), bias(m_rgen), bias(m_rgen));
//					travelling_dist *= dt_for_wind;
//					glm::uvec2 indices = static_cast<glm::uvec2>(glm::floor(static_cast<glm::dvec2>(travelling_dist))) + glm::uvec2(i, j);
//					if (m_terrain.in_range(indices))
//						m_terrain[indices.x][indices.y].sand += eroded;
//				}
//				else
//				{
//					double travelling_sand = std::max(glm::dot(c.normal, -wind) / (c.nutrition + 1.),
//						0.0001 / (c.nutrition + 1.))*wind_sand_coeff;
//					travelling_sand *= c.sand_water_ratio();
//					if (c.sand < travelling_sand)
//						travelling_sand = c.sand;
//					c.sand -= travelling_sand;
//					glm::dvec3 travelling_dist = glm::normalize(wind)*travelling_distance(m_rgen) +
//						glm::dvec3(bias(m_rgen), bias(m_rgen), bias(m_rgen));
//					travelling_dist *= dt_for_wind;
//					glm::uvec2 indices = static_cast<glm::uvec2>(glm::floor(static_cast<glm::dvec2>(travelling_dist))) + glm::uvec2(i, j);
//					if (m_terrain.in_range(indices))
//						m_terrain[indices.x][indices.y].sand += travelling_sand;
//				}
//			}
//		}
//	}
//}
//
//void terrain_generator::sun(glm::dvec3 rays, double duration)
//{
//	constexpr double coeff0=0.001;
//	m_terrain.calc_heights();
//	m_terrain.calc_normals();
//
//	glm::uvec2 size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			auto& c = m_terrain[i][j];
//			if (c.height > 0.)
//			{
//				double vaporazition_factor = coeff0*std::max(glm::dot(-rays, c.normal), 0.)*duration;
//				vaporazition_factor = exp(-vaporazition_factor);
//				c.water *= vaporazition_factor;
//			}
//		}
//	}
//}
//
//void terrain_generator::generate()
//{
//	init_terrain();
//
//	m_terrain.calc_heights();
//	m_terrain.calc_normals();
//	show_normal();
//
//	for (uint32_t i = 0; i < 900; ++i)
//	{
//		show_height();
//		elevate_terrain();
//	}
//
//	for (uint32_t i = 0; i < 5; ++i)
//	{
//		//elevate_terrain();
//		//show_water();
//
//		std::normal_distribution<> dis_ray(0., 10.);
//		glm::dvec3 ray(dis_ray(m_rgen), dis_ray(m_rgen), dis_ray(m_rgen));
//		std::normal_distribution<> dis_duration(30., 10.);
//		double duration = dis_duration(m_rgen);
//		sun(ray, duration);
//		wind(ray, 10);
//		//show_sand();
//
//		//std::this_thread::sleep_for(std::chrono::seconds(3));
//	}
//
//	while (false)
//	{
//		std::normal_distribution<> dis_ray(6., 5.);
//		glm::dvec3 ray(dis_ray(m_rgen), dis_ray(m_rgen), dis_ray(m_rgen));
//		wind(ray, 10);
//		show_height();
//	}
//
//	write_to_file("output/t.terr");
//}
//
//void terrain_generator::show_height()
//{
//	sf::Event e;
//	while (m_window.pollEvent(e));
//
//	m_window.clear(sf::Color(100, 100, 100));
//
//	sf::RectangleShape r(sf::Vector2f(1, 1));
//
//	glm::uvec2 size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			int height = std::max(0, static_cast<int>(m_terrain[i][j].rock_height));
//			r.setPosition(j, i);
//			r.setFillColor(sf::Color(height, height, height, 255));
//			m_window.draw(r);
//		}
//	}
//	m_window.display();
//}
//
//void terrain_generator::show_water()
//{
//	sf::Event e;
//	while (m_window.pollEvent(e));
//
//	m_window.clear(sf::Color(100, 100, 100));
//
//	sf::RectangleShape r(sf::Vector2f(1, 1));
//
//	glm::uvec2 size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			int water = std::max(0, static_cast<int>(m_terrain[i][j].water * 2550.));
//			r.setPosition(j, i);
//			r.setFillColor(sf::Color(0, 0, water, 255));
//			m_window.draw(r);
//		}
//	}
//	m_window.display();
//}
//
//void terrain_generator::show_normal()
//{
//	sf::Event e;
//	while (m_window.pollEvent(e));
//
//	m_window.clear(sf::Color(100, 100, 100));
//
//	sf::RectangleShape r(sf::Vector2f(1, 1));
//
//	glm::uvec2 size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			glm::dvec3 n = m_terrain[i][j].normal;
//			n = 0.5*n + glm::dvec3(0.5);
//			n *= 255.;
//			r.setPosition(j, i);
//			r.setFillColor(sf::Color(n.x, n.y, n.z, 255));
//			m_window.draw(r);
//		}
//	}
//	m_window.display();
//}
//
//void terrain_generator::show_sand()
//{
//	sf::Event e;
//	while (m_window.pollEvent(e));
//
//	m_window.clear(sf::Color(100, 100, 100));
//
//	sf::RectangleShape r(sf::Vector2f(1, 1));
//
//	glm::uvec2 size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			double sand = m_terrain[i][j].sand*255.;
//			r.setPosition(j, i);
//			r.setFillColor(sf::Color(sand, sand, 0, 255));
//			m_window.draw(r);
//		}
//	}
//	m_window.display();
//}
//
//void terrain_generator::write_to_file(const std::string_view& filename)
//{
//	std::ofstream file(filename.data(), std::ios::binary);
//
//	auto size = m_terrain.size();
//	for (uint32_t i = 0; i < size.x; ++i)
//	{
//		for (uint32_t j = 0; j < size.y; ++j)
//		{
//			auto& c = m_terrain[i][j];
//			glm::vec4 val = { static_cast<float>(c.rock_height), static_cast<float>(c.water), 
//				static_cast<float>(c.sand), static_cast<float>(c.nutrition) };
//			file.write(reinterpret_cast<char*>(&val), sizeof(glm::vec4));
//		}
//	}
//	file.close();
//}