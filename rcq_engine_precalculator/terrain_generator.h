//#pragma once
//
//#include <vector>
//#include <glm.hpp>
//#include <SFML/Graphics.hpp>
//#include <thread>
//#include <random>
//#include <fstream>
//#include <string_view>
//
//constexpr double sand_density = 1.5;
//constexpr double rock_erosion_by_aggregate = 0.001;
//
//constexpr double sand_water_capacity = 0.25;
//constexpr double sand_water_wilting_point = 0.1;
//
//constexpr double terminal_velocity_for_sand = 0.5;
//
//constexpr double gravitational_acceleration = 9.81;
//
//constexpr double dt_for_rain = 1.;
//
//constexpr double dt_for_wind = 10.;
//
//constexpr double wind_erosion_rock=1e-3;
//constexpr double wind_sand_coeff = 0.1;
//
//struct aggregate
//{
//	double water;
//	double sand;
//	double nutrition;
//	glm::dvec3 velocity;
//	bool exists;
//
//	double mass() const
//	{
//		return water + sand*sand_density;
//	}
//
//	double sand_water_ratio()
//	{
//		return sand / (sand + water);
//	}
//
//	void clear()
//	{
//		water = 0.;
//		sand = 0.;
//		nutrition = 0.;
//		velocity = glm::dvec3(0.);
//	}
//
//};
//
//struct cell
//{
//	double rock_height;
//	double height;
//	double water;
//	double sand;
//	double nutrition;
//
//	//helper fileds
//	glm::dvec3 normal;
//	aggregate aggr;
//
//	aggregate next_aggr;
//
//	static glm::dvec2 iterate_with_mixing_equation(const glm::dvec2 start, double k1, double k2, double dt)
//	{
//		glm::dvec2 coeffs = glm::dmat2(1, k1, 1, -k2)*start / (k1 + k2);
//		return glm::dvec2(k2, k1)*coeffs.x + glm::dvec2(1., -1.)*coeffs.y*exp(-(k1 + k2)*dt);
//	}
//
//	double mass()
//	{
//		return water + sand*sand_density;
//	}
//
//	double sand_water_ratio()
//	{
//		return sand / (sand + water);
//	}
//
//	void collide_with_aggr()
//	{
//		double v_n_length = glm::dot(aggr.velocity, normal);
//		double v_n_length_square = v_n_length*v_n_length;
//		double pick_up_mass;
//		double v_t_length_square = glm::dot(aggr.velocity, aggr.velocity) - v_n_length_square;
//		aggr.velocity -= (normal*v_n_length);
//
//		if (v_t_length_square == 0.)
//		{
//			pick_up_mass = aggr.mass()*v_n_length_square;
//		}
//		else
//		{
//			pick_up_mass = aggr.mass()*v_n_length_square*(1. - 1. / (1. + v_t_length_square)) / v_t_length_square;
//		}
//
//		if (double m = mass(); m< pick_up_mass)
//		{
//			aggr.nutrition += nutrition;
//			nutrition = 0.;
//			aggr.sand += sand;
//			sand = 0.;
//			aggr.water += water;
//			water = 0.;
//			rock_height -= ((pick_up_mass - m)*rock_erosion_by_aggregate);
//		}
//		else
//		{
//			double plus_water = water*pick_up_mass / mass();
//			double plus_nutrition = nutrition*pick_up_mass / mass();
//			double plus_sand = sand*pick_up_mass / mass();
//
//			aggr.water += plus_water;
//			water -= plus_water;
//			aggr.sand += plus_sand;
//			sand -= plus_sand;
//			aggr.nutrition += plus_nutrition;
//			nutrition -= plus_nutrition;
//		}
//	}
//
//	void deposit_from_aggregate()
//	{
//		if (glm::length(aggr.velocity) < terminal_velocity_for_sand*aggr.sand_water_ratio())
//		{
//			//dissolve aggregate
//			nutrition += aggr.nutrition;
//			sand += aggr.sand;
//			water += aggr.water;
//			aggr.clear();
//			return;
//		}
//
//		if (aggr.sand_water_ratio() > sand_water_capacity && sand_water_ratio() < sand_water_capacity)
//		{
//			double solvable_sand = sand_water_capacity*aggr.water / (1. - sand_water_capacity);
//			double deposit_sand = aggr.sand - solvable_sand;
//
//			constexpr double nutrition_factor = 0.05;
//			double deposit_nutrition = nutrition_factor*aggr.nutrition*deposit_sand / (aggr.sand + aggr.water);
//
//			sand += deposit_sand;
//			aggr.sand = solvable_sand;
//
//			nutrition += deposit_nutrition;
//			aggr.nutrition -= deposit_nutrition;
//		}
//		else if (aggr.sand_water_ratio() < sand_water_capacity && sand_water_ratio() > sand_water_capacity)
//		{
//			double solvable_sand = sand_water_capacity*water / (1. - sand_water_capacity);
//			double deposit_sand = sand - solvable_sand;
//
//			constexpr double nutrition_factor = 0.05;
//			double deposit_nutrition = nutrition_factor*nutrition*deposit_sand / (sand + water);
//
//			aggr.sand += deposit_sand;
//			sand = solvable_sand;
//
//			aggr.nutrition += deposit_nutrition;
//			nutrition -= deposit_nutrition;
//		}
//		else if (aggr.sand_water_ratio() <= sand_water_capacity && sand_water_ratio() <= sand_water_capacity)
//		{
//			//double dt = 1. / std::max(glm::length(aggr.velocity), 0.000001);
//
//			//mixing nutritions
//			constexpr double nutrition_extraction_coeff = 0.001;
//
//			double k_nutrition_from_aggregate = nutrition_extraction_coeff*aggr.water / (aggr.water + aggr.sand);
//			double k_nutrition_from_cell = nutrition_extraction_coeff*water / (water + sand);
//
//			glm::dvec2 old_nutrition_contentrations = { aggr.nutrition / (aggr.sand + aggr.water), nutrition / (sand + water) };
//			glm::dvec2 new_nutrition_contentrations = iterate_with_mixing_equation(old_nutrition_contentrations,
//				k_nutrition_from_aggregate, k_nutrition_from_cell, dt_for_rain);
//
//			nutrition = new_nutrition_contentrations.y*(sand + water);
//			aggr.nutrition = new_nutrition_contentrations.x*(aggr.sand + aggr.water);
//		}
//	}
//
//	void add_to_next_aggregate(const aggregate& a, double ratio)
//	{
//		if (ratio == 0.)
//			return;
//
//		if (!next_aggr.exists)
//		{
//			next_aggr.exists = true;
//			next_aggr.nutrition = a.nutrition*ratio;
//			next_aggr.sand = a.sand*ratio;
//			next_aggr.water = a.water*ratio;
//			next_aggr.velocity = a.velocity;
//		}
//		else
//		{
//			next_aggr.velocity = (next_aggr.mass()*next_aggr.velocity + a.mass()*a.velocity*ratio)*(next_aggr.mass() + a.mass()*ratio);
//			next_aggr.nutrition += a.nutrition*ratio;
//			next_aggr.sand += a.sand*ratio;
//			next_aggr.water += a.water*ratio;
//		}
//	}
//};
//
//
//class terrain
//{
//public:
//	terrain(const glm::uvec2& size) : m_data(size.x)
//	{
//		for (auto& v : m_data)
//			v.resize(size.y); 
//		m_size = size;
//	}
//	glm::uvec2 size()
//	{
//		return m_size;
//	}
//	glm::uvec2 center()
//	{
//		return m_size / 2u;
//	}
//	std::vector<cell>& operator[](size_t i)
//	{
//		return m_data[i];
//	}
//
//	bool in_range(const glm::uvec2& index) const
//	{
//		return index.x >= 0 && index.y >= 0 && index.x < m_size.x && index.y < m_size.y;
//	}
//
//	void calc_heights()
//	{
//		for (uint32_t i = 0; i < m_data.size(); ++i)
//		{
//			for (uint32_t j = 0; j < m_data[i].size(); ++j)
//			{
//				auto& c = m_data[i][j];
//				c.height = c.rock_height + c.water + c.sand;
//			}
//		}
//	}
//
//	void calc_normals()
//	{
//		for (uint32_t i = 0; i < m_data.size(); ++i)
//		{
//			for (uint32_t j = 0; j < m_data[i].size(); ++j)
//			{
//				if (i > 0 && i < m_data.size() - 1 && j>0 && j < m_data[i].size() - 1)
//				{
//					glm::dvec3 normal = {
//						(m_data[i + 1][j].height - m_data[i - 1][j].height) / 4.,
//						(m_data[i][j + 1].height - m_data[i][j - 1].height) / 4.,
//						1.
//					};
//					m_data[i][j].normal = glm::normalize(normal);
//				}
//				else
//					m_data[i][j].normal = { 0., 0., 1. };
//			}
//		}
//	}
//
//	void update_velocities()
//	{
//		for (uint32_t i = 0; i < m_data.size(); ++i)
//		{
//			for (uint32_t j = 0; j < m_data[i].size(); ++j)
//			{
//				glm::dvec2 grad = { (m_data[i + 1][j].height - m_data[i - 1][j].height) / 4.,
//					(m_data[i][j + 1].height - m_data[i][j - 1].height) / 4.
//				};
//
//				glm::dvec3 acceleration(-grad, -glm::dot(grad, grad));
//
//				m_data[i][j].aggr.velocity += (acceleration*dt_for_rain);
//			}
//		}
//	}
//
//	void collide_with_aggreagates()
//	{
//		for (uint32_t i = 0; i < m_data.size(); ++i)
//		{
//			for (uint32_t j = 0; j < m_data[i].size(); ++j)
//			{
//				m_data[i][j].collide_with_aggr();
//			}
//		}
//	}
//
//	void deposit_from_aggregates()
//	{
//		for (uint32_t i = 0; i < m_data.size(); ++i)
//		{
//			for (uint32_t j = 0; j < m_data[i].size(); ++j)
//			{
//				m_data[i][j].deposit_from_aggregate();
//			}
//		}
//	}
//
//	void distribute_aggregates()
//	{
//		for (uint32_t i = 0; i < m_data.size(); ++i)
//		{
//			for (uint32_t j = 0; j < m_data[i].size(); ++j)
//			{
//				if (!m_data[i][j].aggr.exists)
//					continue;
//
//				glm::dvec2 pos = { static_cast<double>(i), static_cast<double>(j) };
//				pos += static_cast<glm::dvec2>(m_data[i][j].aggr.velocity)*dt_for_rain;
//
//				glm::dvec2 center = glm::ceil(pos);
//				glm::uvec2 index = static_cast<glm::uvec2>(center);
//				
//				if (in_range({index.x-1, index.y-1}))
//					m_data[index.x-1][index.y-1].add_to_next_aggregate(m_data[i][j].aggr, (center.x-pos.x)*(center.y-pos.y));
//				if (in_range({ index.x, index.y - 1 }))
//					m_data[index.x][index.y - 1].add_to_next_aggregate(m_data[i][j].aggr, (pos.x + 1. - center.x)*(center.y - pos.y));
//				if (in_range({ index.x - 1, index.y }))
//					m_data[index.x - 1][index.y].add_to_next_aggregate(m_data[i][j].aggr, (center.x - pos.x)*(pos.y + 1. - center.y));
//				if (in_range({ index.x, index.y }))
//					m_data[index.x][index.y].add_to_next_aggregate(m_data[i][j].aggr, (pos.x + 1. - center.x)*(pos.y + 1. - center.y));
//
//			}
//		}
//	}
//
//	bool set_next_aggregate_to_aggregate()
//	{
//		bool aggregate_exists = false;
//
//		for (uint32_t i = 0; i < m_data.size(); ++i)
//		{
//			for (uint32_t j = 0; j < m_data[i].size(); ++j)
//			{
//				if (m_data[i][j].next_aggr.exists)
//				{
//					aggregate_exists = true;
//					m_data[i][j].aggr = m_data[i][j].next_aggr;
//					m_data[i][j].next_aggr.exists = false;
//				}
//			}
//		}
//		return aggregate_exists;
//	}
//
//private:
//	std::vector<std::vector<cell>> m_data;
//	glm::uvec2 m_size;
//};
//
//class terrain_generator
//{
//public:
//	terrain_generator(const glm::uvec2& size);
//	~terrain_generator();
//
//	void generate();
//
//private:
//	enum
//	{
//		WIND,
//		RAIN,
//		NUTRITION,
//		COAST,
//		SUN
//	};
//
//	void init_terrain();
//	void elevate_terrain();
//	void wind(const glm::dvec3& wind, uint32_t iteration);
//	void rain(double amount);
//	void nutrition();
//	void coast();
//	void sun(glm::dvec3 rays, double duration);
//
//	void show_height();
//	void show_water();
//	void show_normal();
//	void show_sand();
//
//	void write_to_file(const std::string_view& filename);
//
//
//	terrain m_terrain;
//	//tectonic planes params
//	glm::dmat2 m_tectonic_planes_covariance;
//	double m_tectonic_planes_normalisation_factor;
//	double m_boundary_amplitude;
//	double m_boundary_frequency;
//	double m_elevation_speed;
//
//	//erosion params
//	double m_sand_terminal_velocity;
//	double m_high_tide;
//
//	//random
//	std::random_device m_rd;
//	std::mt19937 m_rgen;
//
//	//SFML
//	sf::RenderWindow m_window;
//};
//
