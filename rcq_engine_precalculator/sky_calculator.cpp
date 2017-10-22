#include "sky_calculator.h"

# define M_PI           3.14159265358979323846

const double sky_calculator::Rayleigh_scale_height=8000.;
const double sky_calculator::Mie_scale_height=1200.;
const double sky_calculator::index_of_refraction=1.0003;
const double sky_calculator::Earth_radius=6371000;
const double sky_calculator::Atmosphere_thickness=80000;
const double sky_calculator::Atmosphere_radius=Earth_radius+sky_calculator::Atmosphere_thickness;


sky_calculator::sky_calculator(size_t size)
{
	m_image.resize(size);
	m_intensity.resize(size);
	m_accumulated_intensities.resize(size);
	m_gathered_scattered_light.resize(size);
	m_size = static_cast<double>(size);
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
	double a = start.y - Earth_radius;
	double b = end.y - Earth_radius;
	return ((exp(a / Rayleigh_scale_height) - exp(b / Rayleigh_scale_height))*Rayleigh_scale_height / (b - a))*
		m_Rayleigh_scattering_coefficient;
}

double sky_calculator::Mie_transmittance(const glm::dvec3& start, const glm::dvec3& end)
{
	double a = start.y - Earth_radius;
	double b = end.y - Earth_radius;
	return ((exp(a / Mie_scale_height) - exp(b / Mie_scale_height))*Mie_scale_height / (b - a))*
		m_Mie_scattering_coefficient;
}

glm::dvec3 sky_calculator::Rayleigh_single_scattering(const glm::dvec3& params)
{
	static const double sample_count = 100.;

	double height = params.x;
	double  cos_view_zenit = params.y;
	double cos_sun_zenit = params.z;

	glm::dvec3 P0(0., height + Earth_radius, 0.);
	glm::dvec3 v(sqrt(1.-cos_view_zenit*cos_view_zenit), cos_view_zenit, 0.);
	glm::dvec3 l(sqrt(1. - cos_sun_zenit*cos_sun_zenit), cos_sun_zenit, 0.);

	double T = atmosphere_intersection(P0, v);

	double dt = T / (sample_count - 1.);
	glm::dvec3 step = dt*v;

	glm::dvec3 P = P0;
	double S = atmosphere_intersection(P, l);
	glm::dvec3 prev_val = exp(-Rayleigh_transmittance(P, P + S*l));

	glm::vec3 res(0.);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		glm::dvec3 f1(-(P.y - Earth_radius) / Rayleigh_scale_height);
		glm::dvec3 f2 = -Rayleigh_transmittance(P0, P);
		S = atmosphere_intersection(P, l);
		glm::dvec3 f3 = -Rayleigh_transmittance(P, P + S*l);
		glm::dvec3 current_value = exp(f1 + f2 + f3);

		glm::dvec3 area = dt*(current_value + prev_val) / 2.;
		res += area;
		prev_val = current_value;
	}
	res *= (m_Rayleigh_scattering_coefficient / (4.*M_PI));

	return res;
}

double sky_calculator::Mie_single_scattering(const glm::dvec3& params)
{
	static const double sample_count = 100.;

	double height = params.x;
	double  cos_view_zenit = params.y;
	double cos_sun_zenit = params.z;


	glm::dvec3 P0(0., height + Earth_radius, 0.);
	glm::dvec3 v(sqrt(1.-cos_view_zenit*cos_view_zenit), cos_view_zenit, 0.);
	glm::dvec3 l(sqrt(1.-cos_sun_zenit*cos_sun_zenit), cos_sun_zenit, 0.);

	double T = atmosphere_intersection(P0, v);

	double dt = T / (sample_count - 1.);
	glm::dvec3 step = dt*v;

	glm::dvec3 P = P0;
	double S = atmosphere_intersection(P, l);
	double prev_val = exp(-Mie_transmittance(P, P + S*l));

	double res(0.);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		double f1=-(P.y - Earth_radius) / Mie_scale_height;
		double f2 = -Mie_transmittance(P0, P);
		S = atmosphere_intersection(P, l);
		double f3 = -Mie_transmittance(P, P + S*l);
		double current_value = exp(f1 + f2 + f3);

		double area = dt*(current_value + prev_val) / 2.;
		res += area;
		prev_val = current_value;
	}
	res *= (m_Mie_scattering_coefficient / (4.*M_PI));

	return res;
}

double sky_calculator::atmosphere_intersection(const glm::dvec3 & P0, const glm::dvec3 & v)
{
	double v_dot_P0 = glm::dot(P0, v);
	return -v_dot_P0 + sqrt(v_dot_P0*v_dot_P0 - glm::dot(P0, P0) + Atmosphere_radius*Atmosphere_radius);
}

glm::dvec3 sky_calculator::map_indices(const glm::dvec3& indices)
{
	glm::dvec3 tex_coords =indices/(m_size-1.);

	glm::dvec3 ret;
	ret.x = Atmosphere_thickness*tex_coords.x*tex_coords.x;

	double cos_horizon = -sqrt(ret.x*(2 * Earth_radius + ret.x)) / (Earth_radius + ret.x);

	if (tex_coords.y > 0.5)
	{
		ret.y = pow(2.*tex_coords.y - 1., 5.)*(1.f - cos_horizon) + cos_horizon;
	}
	else
	{
		ret.y= -pow(2.*tex_coords.y, 5.)*(1.f + cos_horizon) + cos_horizon;
	}
	ret.y=glm::clamp(ret.y, -1., 1.);

	ret.z = tan((2.*tex_coords.z - 0.74)*1.1) / tan(1.26*1.1);
	ret.z = glm::clamp(ret.z, -1., 1.);

	return ret;
}


void sky_calculator::do_first_iteration()
{
	//fill intensitiy
	uint32_t size = static_cast<uint32_t>(m_size);
	for (uint32_t i = 0; i < size; ++i)
	{
		for (uint32_t j = 0; j < size; ++j)
		{
			for (uint32_t k = 0; k < size; ++k)
			{
				glm::dvec3 indices = { static_cast<double>(k), static_cast<double>(j), static_cast<double>(i) };
				glm::dvec3 params = map_indices(indices);
				glm::dvec3 Rayleigh = Rayleigh_single_scattering(params);
				glm::dvec3 Mie(Mie_single_scattering(params));
				m_intensity[i][j][k] = Mie + Rayleigh;
			}
		}

	}

	//fill gathered scattered light
}
