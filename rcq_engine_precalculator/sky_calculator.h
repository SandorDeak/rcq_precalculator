#pragma once

#include <vector>
#include <glm.hpp>

class image3d
{
public:
	void resize(size_t size)
	{
		data.resize(size);
		for (auto& p : data)
		{
			p.resize(size);
			for (auto& q : p)
				q.resize(size);
		}
	}
private:
	std::vector<std::vector<std::vector<glm::dvec3>>> data;
};


class sky_calculator
{
public:
	sky_calculator(size_t size);
	~sky_calculator();

private:
	static const double Rayleigh_scale_height;
	static const double Mie_scale_height;
	static const double index_of_refraction;

	static double Rayleigh_phase(double theta);
	static double Mie_phase(double theta, double g);
	static double Rayleigh_density(double height);
	static double Mie_density(double height);

	glm::dvec3 Rayleigh_scattering_intensity(double cos_theta, double height);
	double Mie_scattering_intensity(double cos_theta, double height);

	glm::dvec3 Rayleigh_transmittance(const glm::dvec3& start, const glm::dvec3& end);
	double Mie_transmittance(const glm::dvec3& start, const glm::dvec3& end);

	glm::dvec3 Rayleigh_single_scattering(double altitude, double cos_view_zenit, double cos_sun_zeit);
	double Mie_single_scattering(double altitude, double cos_view_zenit, double cos_sun_zenit);


	image3d m_image;
	size_t m_size;

	//params
	double m_Rayleigh_particle_density; //at sea level, mol/m^3
	double m_Mie_particle_density; //at sea level, mol/m^3
	double m_assymetry_factor; //for Rayleigh scattering
	
	glm::dvec3 m_Rayleigh_scattering_coefficient;
	double m_Mie_scattering_coefficient;
};

