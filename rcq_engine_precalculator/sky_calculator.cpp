#include "sky_calculator.h"
#include <fstream>
#include <algorithm>
#include <string.h>
#include <gtx\transform.hpp>

# define PI           3.14159265358979323846
# define RED			680.e-9
# define GREEN			550.e-9
# define BLUE			440.e-9

const double sky_calculator::Rayleigh_scale_height=8000.;
const double sky_calculator::Mie_scale_height=1200.;
const double sky_calculator::index_of_refraction=1.0003;
const double sky_calculator::Earth_radius=6371000;
const double sky_calculator::Atmosphere_thickness=80000;
const double sky_calculator::Atmosphere_radius=Earth_radius+sky_calculator::Atmosphere_thickness;


sky_calculator::sky_calculator(size_t size) : 
	m_window(sf::VideoMode(size,size), "Debug info")
{
	m_image.create(size, size, sf::Color());
	m_image_sprite.setTexture(m_texture, true);
	m_intensity.resize(size);
	m_accumulated_intensities.resize(size);
	m_gathered_scattered_light.resize(size);
	m_size = static_cast<double>(size);
}


sky_calculator::~sky_calculator()
{

}

void sky_calculator::set_params(double Rayleigh_particle_density, double Mie_particle_density, double assymetry_factor)
{
	m_assymetry_factor = assymetry_factor;

	double const1 = 2.*PI*PI*pow(index_of_refraction*index_of_refraction - 1., 2.) / 3.;
	double Rayleigh_particle_polarisability = const1 / pow(Rayleigh_particle_density, 2.);
	double Mie_particle_polarisability = const1 / pow(Mie_particle_density, 2.);

	m_Rayleigh_scattering_coefficient = {
		4.*PI*Rayleigh_particle_polarisability*Rayleigh_particle_density / pow(RED, 4.),
		4.*PI*Rayleigh_particle_polarisability*Rayleigh_particle_density / pow(GREEN, 4.),
		4.*PI*Rayleigh_particle_polarisability*Rayleigh_particle_density / pow(BLUE, 4.)
	};

	m_Mie_scattering_coefficient = 4.*PI*Mie_particle_polarisability*Mie_particle_density;
}

void sky_calculator::load_Lebedev_params(const std::string_view& filename)
{
	std::ifstream file(filename.data());
	char line[512];
	while (file.getline(line, 512))
	{
		char* begin=line;
		double theta = strtod(begin, &begin)*PI/180.;
		double phi = strtod(begin, &begin)*PI/180.;
		double w = strtod(begin, &begin)*PI/180.;

		glm::dvec3 p = {
			cos(theta)*sin(phi),
			sin(theta)*sin(phi),
			cos(phi)
		};
		m_Lebedev_coords.emplace_back(p, w);
	}
	file.close();
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

glm::dvec3 sky_calculator::Rayleigh_transmittance(const glm::dvec3& P0, const glm::dvec3& P1)
{
	constexpr double sample_count = 100.;
	constexpr double dt = 1. / (sample_count - 1.);
	glm::dvec3 v = P1 - P0;
	glm::dvec3 P = P0;
	glm::dvec3 step = dt*v;
	double sum = 0.;
	double prev_val = exp((Earth_radius-glm::length(P))/Rayleigh_scale_height);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		double curr_val = exp((Earth_radius - glm::length(P)) / Rayleigh_scale_height);
		double area = dt*(prev_val + curr_val) / 2.;
		sum += area;
		prev_val = curr_val;
	}

	return sum*glm::length(v)*m_Rayleigh_scattering_coefficient;
}

double sky_calculator::Mie_transmittance(const glm::dvec3& P0, const glm::dvec3& P1)
{
	constexpr double sample_count = 100.;
	constexpr double dt = 1. / (sample_count - 1.);
	glm::dvec3 v = P1 - P0;
	glm::dvec3 P = P0;
	glm::dvec3 step = dt*v;
	double sum = 0.;
	double prev_val = exp((Earth_radius - glm::length(P)) / Mie_scale_height);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		double curr_val = exp((Earth_radius - glm::length(P)) / Mie_scale_height);
		double area = dt*(prev_val + curr_val) / 2.;
		sum += area;
		prev_val = curr_val;
	}

	return sum*glm::length(v)*m_Mie_scattering_coefficient;
}

std::pair<double, double> sky_calculator::calc_integration_interval(const glm::dvec3 & P0, const glm::dvec3 & v)
{
	auto T_Atmo = atmosphere_intersection(P0, v);
	double T_start, T_end;

	if (!T_Atmo)
		return { 0., 0. };
	if (T_Atmo.value().second <= 0.f)
		return { 0., 0. };

	T_start = std::max(0., T_Atmo.value().first);
	T_end = T_Atmo.value().second;

	auto T_Earth = Earth_intersection(P0, v);
	if (T_Earth)
	{
		if (T_Earth.value().first <= 0. && T_Earth.value().second > 0.) //inside Earth
			return { 0., 0. };
		if (T_Earth.value().first > 0.)
			T_end = std::min(T_end, T_Earth.value().first);
	}
	return { T_start, T_end };
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

	auto[T_start, T_end] = calc_integration_interval(P0, v);
	if (T_start == T_end)
		return glm::dvec3(0.);


	double dt = (T_end-T_start) / (sample_count - 1.);
	glm::dvec3 step = dt*v;
	P0 = P0 + T_start*step;
	glm::dvec3 P = P0;

	glm::dvec3 prev_val;	
	if (auto[S_start, S_end] = calc_integration_interval(P, l); S_start == S_end)
		prev_val = glm::dvec3(0.f);
	else
		prev_val = exp(-Rayleigh_transmittance(P+S_start*l, P + S_end*l));

	glm::vec3 res(0.);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		glm::dvec3 f1(-(glm::length(P) - Earth_radius) / Rayleigh_scale_height);
		glm::dvec3 f2 = -Rayleigh_transmittance(P0, P);
		glm::dvec3 f3;
		if (auto[S_start, S_end] = calc_integration_interval(P, l); S_start == S_end)
			f3 = glm::dvec3(0.f);
		else
			f3 = -Rayleigh_transmittance(P + S_start*l, P + S_end*l);
		glm::dvec3 current_value = exp(f1 + f2 + f3);

		glm::dvec3 area = dt*(current_value + prev_val) / 2.;
		res += area;
		prev_val = current_value;
	}
	res *= (m_Rayleigh_scattering_coefficient / (4.*PI));

	return res;
}

double sky_calculator::Mie_single_scattering(const glm::dvec3& params)
{
	static const double sample_count = 100.;

	double height = params.x;
	double  cos_view_zenit = params.y;
	double cos_sun_zenit = params.z;

	glm::dvec3 P0(0., height + Earth_radius, 0.);
	glm::dvec3 v(sqrt(1. - cos_view_zenit*cos_view_zenit), cos_view_zenit, 0.);
	glm::dvec3 l(sqrt(1. - cos_sun_zenit*cos_sun_zenit), cos_sun_zenit, 0.);

	auto[T_start, T_end] = calc_integration_interval(P0, v);
	if (T_start == T_end)
		return 0.;


	double dt = (T_end - T_start) / (sample_count - 1.);
	glm::dvec3 step = dt*v;
	P0 = P0 + T_start*step;
	glm::dvec3 P = P0;

	double prev_val;
	if (auto[S_start, S_end] = calc_integration_interval(P, l); S_start == S_end)
		prev_val = 0.;
	else
		prev_val = exp(-Mie_transmittance(P + S_start*l, P + S_end*l));

	double res(0.);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		double f1=-(glm::length(P) - Earth_radius) / Mie_scale_height;
		double f2 = -Mie_transmittance(P0, P);
		double f3;
		if (auto[S_start, S_end] = calc_integration_interval(P, l); S_start == S_end)
			f3 = 0.;
		else
			f3 = -Mie_transmittance(P + S_start*l, P + S_end*l);
		double current_value = exp(f1 + f2 + f3);

		double area = dt*(current_value + prev_val) / 2.;
		res += area;
		prev_val = current_value;
	}
	res *= (m_Mie_scattering_coefficient / (4.*PI));

	return res;
}

std::optional<std::pair<double, double>> sky_calculator::atmosphere_intersection(const glm::dvec3 & P0, const glm::dvec3 & v)
{
	std::optional<std::pair<double, double>> ret;

	double v_dot_P0 = glm::dot(P0, v);
	if (double d = v_dot_P0*v_dot_P0 - glm::dot(P0, P0) + Atmosphere_radius*Atmosphere_radius; d > 0.)
	{
		double sqrt_d = sqrt(d);
		double t1 = -v_dot_P0 - sqrt_d;
		double t2 = -v_dot_P0 + sqrt_d;
		ret.emplace(t1, t2);
	}
	return ret;
}

std::optional<std::pair<double, double>> sky_calculator::Earth_intersection(const glm::dvec3 & P0, const glm::dvec3 & v)
{
	std::optional<std::pair<double, double>> ret;

	double v_dot_P0 = glm::dot(P0, v);
	if (double d = v_dot_P0*v_dot_P0 - glm::dot(P0, P0) + Earth_radius*Earth_radius; d > 0.)
	{
		double sqrt_d = sqrt(d);
		double t1 = -v_dot_P0 - sqrt_d;
		double t2 = -v_dot_P0 + sqrt_d;
		ret.emplace(t1, t2);
	}
	return ret;
}

glm::dvec3 sky_calculator::indices_to_params(const glm::dvec3& indices)
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

glm::dvec3 sky_calculator::params_to_indices(const glm::dvec3& params)
{
	glm::dvec3 indices;
	indices.x = sqrt(params.x / Atmosphere_thickness);
	
	double cos_horizon= -sqrt(params.x*(2 * Earth_radius + params.x)) / (Earth_radius + params.x);
	if (params.y > cos_horizon)
	{
		indices.y = 0.5*pow((params.y - cos_horizon) / (1. - cos_horizon), 0.2) + 0.5;
	}
	else
	{
		indices.y= 0.5*pow((cos_horizon - params.y) / (1. + cos_horizon), 0.2);
	}

	indices.z = 0.5*(atan(std::max(params.z, -0.1975)*tan(1.26*1.1)) / 1.1 + 0.74);

	return indices*(m_size-1.);
}


glm::dvec3 sky_calculator::sample_intensity(const glm::dvec3& params)
{
	glm::dvec3 indices = params_to_indices(params);
	uint32_t i = static_cast<uint32_t>(floor(indices.z));
	uint32_t j = static_cast<uint32_t>(floor(indices.y));
	uint32_t k = static_cast<uint32_t>(floor(indices.x));

	double ai = indices.z - ceil(indices.z);
	double aj = indices.y - ceil(indices.y);
	double ak = indices.x - ceil(indices.x);
	ai = 1. - ai;
	aj = 1. - aj;
	ak = 1. - ak;

	uint32_t iplus = std::min(m_intensity.size()-1, i + 1);
	uint32_t jplus = std::min(m_intensity.size()-1, j + 1);
	uint32_t kplus = std::min(m_intensity.size()-1, k + 1);

	return
		ai*aj*ak*m_intensity[i][j][k] +
		ai*aj*(1. - ak)*m_intensity[i][j][kplus] +
		ai*(1. - aj)*ak*m_intensity[i][jplus][k] +
		ai*(1. - aj)*(1. - ak)*m_intensity[i][jplus][kplus] +
		(1. - ai)*aj*ak*m_intensity[i][j][k] +
		(1. - ai)*aj*(1. - ak)*m_intensity[iplus][j][kplus] +
		(1. - ai)*(1. - aj)*ak*m_intensity[iplus][jplus][k] +
		(1. - ai)*(1. - aj)*(1. - ak)*m_intensity[iplus][jplus][kplus];
}

glm::dvec3 sky_calculator::sample_gathered_scattered(const glm::dvec3& params)
{
	glm::dvec3 indices = params_to_indices(params);
	uint32_t i = static_cast<uint32_t>(floor(indices.z));
	uint32_t j = static_cast<uint32_t>(floor(indices.y));
	uint32_t k = static_cast<uint32_t>(floor(indices.x));

	double ai = indices.z - ceil(indices.z);
	double aj = indices.y - ceil(indices.y);
	double ak = indices.x - ceil(indices.x);
	ai = 1. - ai;
	aj = 1. - aj;
	ak = 1. - ak;

	uint32_t iplus = std::min(m_gathered_scattered_light.size()-1, i + 1);
	uint32_t jplus = std::min(m_gathered_scattered_light.size()-1, j + 1);
	uint32_t kplus = std::min(m_gathered_scattered_light.size()-1, k + 1);

	return
		ai*aj*ak*m_gathered_scattered_light[i][j][k] +
		ai*aj*(1. - ak)*m_gathered_scattered_light[i][j][kplus] +
		ai*(1. - aj)*ak*m_gathered_scattered_light[i][jplus][k] +
		ai*(1. - aj)*(1. - ak)*m_gathered_scattered_light[i][jplus][kplus] +
		(1. - ai)*aj*ak*m_gathered_scattered_light[i][j][k] +
		(1. - ai)*aj*(1. - ak)*m_gathered_scattered_light[iplus][j][kplus] +
		(1. - ai)*(1. - aj)*ak*m_gathered_scattered_light[iplus][jplus][k] +
		(1. - ai)*(1. - aj)*(1. - ak)*m_gathered_scattered_light[iplus][jplus][kplus];
}



glm::dvec3 sky_calculator::Rayleigh_gathered_scattered(const glm::dvec3& params)
{
	double height = params.x;
	double  cos_view_zenit = params.y;
	double cos_sun_zenit = params.z;

	glm::dvec3 P(0., height + Earth_radius, 0.);
	glm::dvec3 v(sqrt(1. - cos_view_zenit*cos_view_zenit), cos_view_zenit, 0.);
	glm::dvec3 l(sqrt(1. - cos_sun_zenit*cos_sun_zenit), cos_sun_zenit, 0.);

	glm::dvec3 sum(0.);
	for (uint32_t i = 0; i < m_Lebedev_coords.size(); ++i)
	{
		glm::dvec3 omega = m_Lebedev_coords[i].first;
		double w = m_Lebedev_coords[i].second;

		sum += (w*Rayleigh_phase(glm::dot(v, omega))*
			Rayleigh_phase(glm::dot(l, omega))*sample_intensity({params.x, omega.y, params.z}));
	}

	return sum;
}

glm::dvec3 sky_calculator::Mie_gathered_scattered(const glm::dvec3& params)
{
	double height = params.x;
	double  cos_view_zenit = params.y;
	double cos_sun_zenit = params.z;

	glm::dvec3 P(0., height + Earth_radius, 0.);
	glm::dvec3 v(sqrt(1. - cos_view_zenit*cos_view_zenit), cos_view_zenit, 0.);
	glm::dvec3 l(sqrt(1. - cos_sun_zenit*cos_sun_zenit), cos_sun_zenit, 0.);

	glm::dvec3 sum(0.);
	for (uint32_t i = 0; i < m_Lebedev_coords.size(); ++i)
	{
		glm::dvec3 omega = m_Lebedev_coords[i].first;
		double w = m_Lebedev_coords[i].second;

		sum += (w*Mie_phase(glm::dot(v, omega), m_assymetry_factor)*
			Mie_phase(glm::dot(l, omega), m_assymetry_factor)*sample_intensity({ params.x, omega.y, params.z }));
	}
	return sum;
}

glm::dvec3 sky_calculator::Rayleigh_multiple_scattering(const glm::dvec3& params)
{
	/*static const double sample_count = 100.;

	double height = params.x;
	double  cos_view_zenit = params.y;
	double cos_sun_zenit = params.z;


	glm::dvec3 P0(0., height + Earth_radius, 0.);
	glm::dvec3 v(sqrt(1. - cos_view_zenit*cos_view_zenit), cos_view_zenit, 0.);
	glm::dvec3 l(sqrt(1. - cos_sun_zenit*cos_sun_zenit), cos_sun_zenit, 0.);

	double T = atmosphere_intersection(P0, v);

	double dt = T / (sample_count - 1.);
	glm::dvec3 step = dt*v;

	glm::dvec3 P = P0;
	glm::dvec3 prev_val = sample_gathered_scattered(params);

	glm::dvec3 res(0.);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		double rot_angle=asin(glm::cross(P, P0).z/(glm::length(P)*glm::length(P0)));
		glm::dmat3 rot = static_cast<glm::dmat3>(glm::rotate(rot_angle, glm::dvec3(0., 0, 1.)));
		glm::dvec3 light = rot*l;
		glm::dvec3 view = rot*v;

		double view_zenit = view.y;
		double sun_zenit = light.y;
		double new_height = P.y - Earth_radius;
		glm::dvec3 gathered_scattered = sample_gathered_scattered({ new_height, view_zenit, sun_zenit });
		double f1 = -P.y / Rayleigh_scale_height;
		glm::dvec3 f2 = -Rayleigh_transmittance(P0, P);
		glm::dvec3 new_value = gathered_scattered*exp(f1 + f2);
		glm::dvec3 area = dt*(prev_val + new_value) / 2.;
		res += area;
		prev_val = new_value;
	}

	return res*m_Rayleigh_scattering_coefficient/4.*PI;*/
	return glm::dvec3();
}
glm::dvec3 sky_calculator::Mie_multiple_scattering(const glm::dvec3& params)
{
	/*static const double sample_count = 100.;

	double height = params.x;
	double  cos_view_zenit = params.y;
	double cos_sun_zenit = params.z;


	glm::dvec3 P0(0., height + Earth_radius, 0.);
	glm::dvec3 v(sqrt(1. - cos_view_zenit*cos_view_zenit), cos_view_zenit, 0.);
	glm::dvec3 l(sqrt(1. - cos_sun_zenit*cos_sun_zenit), cos_sun_zenit, 0.);

	double T = atmosphere_intersection(P0, v);

	double dt = T / (sample_count - 1.);
	glm::dvec3 step = dt*v;

	glm::dvec3 P = P0;
	glm::dvec3 prev_val = sample_gathered_scattered(params);

	glm::dvec3 res(0.);

	for (uint32_t i = 1; i < sample_count; ++i)
	{
		P += step;
		double rot_angle = asin(glm::cross(P, P0).z / (glm::length(P)*glm::length(P0)));
		glm::dmat3 rot = static_cast<glm::dmat3>(glm::rotate(rot_angle, glm::dvec3(0., 0, 1.)));
		glm::dvec3 light = rot*l;
		glm::dvec3 view = rot*v;

		double view_zenit = view.y;
double sun_zenit = light.y;
double new_height = P.y - Earth_radius;
glm::dvec3 gathered_scattered = sample_gathered_scattered({ new_height, view_zenit, sun_zenit });
double f1 = -P.y / Mie_scale_height;
glm::dvec3 f2(-Mie_transmittance(P0, P));
glm::dvec3 new_value = gathered_scattered*exp(f1 + f2);
glm::dvec3 area = dt*(prev_val + new_value) / 2.;
res += area;
prev_val = new_value;
	}

	return res*m_Mie_scattering_coefficient / 4.*PI;*/
	return glm::dvec3();
}

void sky_calculator::calc_single_scattering()
{
	uint32_t size = static_cast<uint32_t>(m_size);
	for (uint32_t i = 0; i < size; ++i)
	{
		for (uint32_t j = 0; j < size; ++j)
		{
			for (uint32_t k = 0; k < size; ++k)
			{
				glm::dvec3 indices = { static_cast<double>(k), static_cast<double>(j), static_cast<double>(i) };
				glm::dvec3 params = indices_to_params(indices);
				glm::dvec3 Rayleigh = Rayleigh_single_scattering(params);
				glm::dvec3 Mie(Mie_single_scattering(params));
				m_intensity[i][j][k] = Mie + Rayleigh;
				m_accumulated_intensities[i][j][k] = Mie + Rayleigh;
			}
		}

		show(m_intensity, i);

	}
}

void sky_calculator::iterate()
{
	//compute gathered scattered light
	uint32_t size = static_cast<uint32_t>(m_size);
	for (uint32_t i = 0; i < size; ++i)
	{
		for (uint32_t j = 0; j < size; ++j)
		{
			for (uint32_t k = 0; k < size; ++k)
			{
				glm::dvec3 indices = { static_cast<double>(k), static_cast<double>(j), static_cast<double>(i) };
				glm::dvec3 params = indices_to_params(indices);
				glm::dvec3 Rayleigh = Rayleigh_gathered_scattered(params);
				glm::dvec3 Mie(Mie_gathered_scattered(params));
				m_gathered_scattered_light[i][j][k] = Rayleigh + Mie;
			}
		}
	}

	// calculate new inensities
	for (uint32_t i = 0; i < size; ++i)
	{
		for (uint32_t j = 0; j < size; ++j)
		{
			for (uint32_t k = 0; k < size; ++k)
			{
				glm::dvec3 indices = { static_cast<double>(k), static_cast<double>(j), static_cast<double>(i) };
				glm::dvec3 params = indices_to_params(indices);
				glm::dvec3 Rayleigh = Rayleigh_multiple_scattering(params);
				glm::dvec3 Mie = Mie_multiple_scattering(params);
				m_intensity[i][j][k] = Rayleigh + Mie;
			}
		}
	}

	//accumulate intensity
	for (uint32_t i = 0; i < size; ++i)
	{
		for (uint32_t j = 0; j < size; ++j)
		{
			for (uint32_t k = 0; k < size; ++k)
			{
				m_accumulated_intensities[i][j][k] += m_intensity[i][j][k];
			}
		}
	}
}

void sky_calculator::compute(uint32_t scattering_count)
{
	calc_single_scattering();

	show(m_intensity, 0);

	for (uint32_t i = 0; i < scattering_count - 1; ++i)
	{
		iterate();
	}
}

void sky_calculator::write_to_file(const std::string_view& filename)
{
	std::ofstream file(filename.data(), std::ios::binary);
	uint32_t size = static_cast<uint32_t>(m_size);

	for (uint32_t i = 0; i < size; ++i)
	{
		for (uint32_t j = 0; j < size; ++j)
		{
			for (uint32_t k = 0; k < size; ++k)
			{
				auto& p = m_accumulated_intensities[i][j][k];
				glm::vec4 q = { static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z), 1.f };
				file.write(reinterpret_cast<char*>(&q), sizeof(glm::vec4));
			}
		}
	}
	std::flush(file);
	file.close();
}

void sky_calculator::show(const image3d& pic, size_t index)
{
	sf::Event e;
	while (m_window.pollEvent(e));

	m_window.clear(sf::Color(100,100,100));
	
	sf::RectangleShape r(sf::Vector2f(1,1));

	for (uint32_t i = 0; i < pic[index].size(); ++i)
	{
		for (uint32_t j = 0; j < pic[index][i].size(); ++j)
		{
			const auto& p = pic[index][i][j];
			r.setPosition(j, i);
			r.setFillColor(sf::Color(255.*p.x, 255.*p.y, 255. * p.z, 255));
			m_window.draw(r);
			auto& vmi = m_image.getPixel(j, i);
			int g = 4 + 9;
		}
	}
	m_window.display();
}

void sky_calculator::load_from_file(const std::string_view & filename, size_t size1)
{
	std::ifstream file(filename.data(), std::ios::ate | std::ios::binary);

	if (!file.is_open())
		throw std::runtime_error("failed to open file: " + std::string(filename));

	size_t size = (size_t)file.tellg();
	std::vector<glm::vec4> data(size/sizeof(glm::vec4));
	file.seekg(0);
	file.read(reinterpret_cast<char*>(data.data()), size);
	file.close();
	uint32_t count = 0;

	m_window.clear(sf::Color(100, 100, 100));

	sf::RectangleShape r(sf::Vector2f(1, 1));
	for (uint32_t i = 0; i < size1; ++i)
	{
		for (uint32_t j = 0; j < size1; ++j)
		{;
			r.setPosition(j, i);
			auto& p = data[count];
			r.setFillColor(sf::Color(255.*p.x, 255.*p.y, 255. * p.z, 255));
			m_window.draw(r);
			auto& vmi = m_image.getPixel(j, i);
			++count;
		}
	}
	m_window.display();
}