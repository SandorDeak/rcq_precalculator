#include "sky_calculator.h"

const double sky_calculator::Rayleigh_scale_height=8000.;
const double sky_calculator::Mie_scale_height=1200.;
const double sky_calculator::index_of_refraction=1.0003;


sky_calculator::sky_calculator(size_t size)
{
	m_image.resize(size);
}


sky_calculator::~sky_calculator()
{

}

double sky_calculator::Rayleigh_phase(double cos_theta)
{
	return (3.*(1. + cos_theta*cos_theta)) / 4.;
}

double sky_calculator::Mie_phase(double cos_theta, double g)
{
	double g_square = g*g;
	double nom = 3. * (1. - g_square)*(1. + cos_theta*cos_theta);
	double denom = 2. * (2. + g_square)*pow(1. + g_square - 2. * g*cos_theta, 1.5);
	return nom / denom;
}

double sky_calculator::Rayleigh_density(double height)
{
	return exp(-height / Rayleigh_scale_height);
}
double sky_calculator::Mie_density(double height)
{
	return exp(-height / Mie_scale_height);
}

glm::dvec3 sky_calculator::Rayleigh_scattering_intensity(double cos_theta, double height)
{
	return (Rayleigh_density(height)*Rayleigh_phase(cos_theta))*m_Rayleigh_scattering_coefficient;
}

double sky_calculator::Mie_scattering_intensity(double cos_theta, double height)
{
	return Mie_density(height)*Mie_phase(cos_theta, m_assymetry_factor)*m_Mie_scattering_coefficient;
}

glm::dvec3 sky_calculator::Rayleigh_transmittance(const glm::dvec3& start, const glm::dvec3& end)
{
	glm::dvec3 v = end - start;
	double length = glm::length(v);
	
	double a1 = -(end.y - start.y) / Rayleigh_scale_height;
	double a0 = -start.y / Rayleigh_scale_height;

	return exp(-abs(length*exp(a0)*(exp(a1) - 1) / a1)*m_Rayleigh_scattering_coefficient);
}

double sky_calculator::Mie_transmittance(const glm::dvec3& start, const glm::dvec3& end)
{
	glm::dvec3 v = end - start;
	double length = glm::length(v);

	double a1 = -(end.y - start.y) / Mie_scale_height;
	double a0 = -start.y / Mie_scale_height;

	return exp(-abs(length*exp(a0)*(exp(a1) - 1) / a1)*m_Mie_scattering_coefficient);
}

glm::dvec3 sky_calculator::Rayleigh_single_scattering(double altitude, double cos_view_zenit, double cos_sun_zeit)
{

}
double sky_calculator::Mie_single_scattering(double altitude, double cos_view_zenit, double cos_sun_zenit)
{

}
