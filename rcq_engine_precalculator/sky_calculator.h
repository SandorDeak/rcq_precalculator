#pragma once

#include <vector>
#include <string_view>
#include <glm.hpp>
#include <SFML/Graphics.hpp>
#include <optional>

class image3d
{
public:
	void resize(size_t size)
	{
		m_data.resize(size);
		for (auto& p : m_data)
		{
			p.resize(size);
			for (auto& q : p)
				q.resize(size);
		}
		m_size = size;
	}
	uint32_t size() const { return m_size; }
	std::vector<std::vector<glm::dvec3>>& operator[](size_t i)
	{
		return m_data[i];
	}
	const std::vector<std::vector<glm::dvec3>>& operator[](size_t i) const
	{
		return m_data[i];
	}

private:
	std::vector<std::vector<std::vector<glm::dvec3>>> m_data;
	uint32_t m_size;
};

class image3d_float
{
public:
	void resize(size_t size)
	{
		m_data.resize(size);
		for (auto& p : m_data)
		{
			p.resize(size);
			for (auto& q : p)
				q.resize(size);
		}
		m_size = size;
	}
	uint32_t size() const { return m_size; }
	std::vector<std::vector<float>>& operator[](size_t i)
	{
		return m_data[i];
	}
	const std::vector<std::vector<float>>& operator[](size_t i) const
	{
		return m_data[i];
	}
private:
	std::vector<std::vector<std::vector<float>>> m_data;
	uint32_t m_size;
};

class image2d
{
public:
	void resize(size_t size)
	{
		m_data.resize(size);
		for (auto& p : m_data)
			p.resize(size);
		m_size = size;
	}
	uint32_t size() const { return m_size; }
	std::vector<glm::dvec3>& operator[](size_t i)
	{
		return m_data[i];
	}
	const std::vector<glm::dvec3>& operator[](size_t i) const
	{
		return m_data[i];
	}

private:
	std::vector<std::vector<glm::dvec3>> m_data;
	uint32_t m_size;
};


class sky_calculator
{
public:
	sky_calculator(size_t size, size_t transmittance_size);
	~sky_calculator();

	void set_params(double Rayleigh_particle_density, double Mie_particle_density, double assymetry_factor);
	void set_params_directly(const glm::dvec3& Rayleigh_scattering_coefficient, double Mie_scattering_coefficient, 
		double assymetry_factor, const glm::dvec3& ozone_absorbtion_coefficient);
	void load_Lebedev_params(const std::string_view& filename);
	void compute(uint32_t scattering_count, const std::string_view& filename);
	void load_from_file(const std::string_view& filename, size_t size);

private:
	static const double Rayleigh_scale_height;
	static const double Mie_scale_height;
	static const double index_of_refraction;
	static const double Earth_radius;
	static const double Atmosphere_radius;
	static const double Atmosphere_thickness;

	static double Rayleigh_phase(double theta);
	static double Mie_phase(double theta, double g);
	static double Rayleigh_density(double height);
	static double Mie_density(double height);
	static std::optional<std::pair<double, double>> atmosphere_intersection(const glm::dvec3& P0, const glm::dvec3& v);
	static std::optional<std::pair<double, double>> Earth_intersection(const glm::dvec3& P0, const glm::dvec3& v);
	std::pair<double, double> calc_integration_interval(const glm::dvec3& P0, const glm::dvec3& v);
	std::pair<double, double> calc_integration_interval_for_light(const glm::dvec3& P0, const glm::dvec3& v);
	void calc_transmittance();
	void write_to_file(const std::string_view& filename);
	void write_transmittance_to_file(const std::string_view& filename);
	
	glm::dvec3 Rayleigh_scattering_intensity(double cos_theta, double height);
	double Mie_scattering_intensity(double cos_theta, double height);

	//glm::dvec3 Rayleigh_transmittance(const glm::dvec3& P0, const glm::dvec3& v);//deprecated
	//double Mie_transmittance(const glm::dvec3& P0, const glm::dvec3& v);//deprecated

	glm::dvec3 transmittance(const glm::dvec3& P0, const glm::dvec3& P1);

	std::pair<glm::dvec3, float> Rayleigh_single_scattering(const glm::dvec3& params);
	glm::dvec3 Mie_single_scattering(const glm::dvec3& params);

	glm::dvec3 Rayleigh_gathered_scattered(const glm::dvec3& params);
	glm::dvec3 Mie_gathered_scattered(const glm::dvec3& params);

	glm::dvec3 Rayleigh_multiple_scattering(const glm::dvec3& params);
	glm::dvec3 Mie_multiple_scattering(const glm::dvec3& params);

	glm::dvec3 indices_to_params(const glm::dvec3& indices);
	glm::dvec3 params_to_indices(const glm::dvec3& params);
	glm::dvec2 tex_coords_to_params(const glm::dvec2& indices);

	glm::dvec3 sample(const image3d& tex, const glm::dvec3& params);

	void calc_single_Rayleigh_scattering();
	void calc_single_Mie_scattering();
	void iterate();
	void show(const image3d& pic, size_t i);
	void show_transmittance();

	image3d m_intensity;
	image3d_float m_ray_length;
	image3d m_accumulated_intensity;
	image3d m_gathered_scattered_light_Rayleigh;
	image3d m_gathered_scattered_light_Mie;
	image2d m_transmittance;
	double m_size;

	std::vector<std::pair<glm::dvec3, double>> m_Lebedev_coords; //x_i y_i z_i and 4*pi*w_i

	//params
	double m_Rayleigh_particle_density; //at sea level, mol/m^3
	double m_Mie_particle_density; //at sea level, mol/m^3
	double m_assymetry_factor; //for Rayleigh scattering
	
	glm::dvec3 m_Rayleigh_scattering_coefficient;
	double m_Mie_scattering_coefficient;
	glm::dvec3 m_Rayleigh_extinction_coefficient; //=Rayleigh scattering coeff
	double m_Mie_extinction_coefficient; //=1.1*Mie scattering coeff
	glm::dvec3 m_ozone_extinction_coefficient;

	//SFML
	sf::RenderWindow m_window;
	sf::RenderWindow m_window_transmittance;
};

