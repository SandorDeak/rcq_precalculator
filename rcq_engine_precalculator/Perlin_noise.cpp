#include "Perlin_noise.h"

#include <random>

# define PI           3.14159265358979323846

Perlin_noise::Perlin_noise()
{
	reset();
}


Perlin_noise::~Perlin_noise()
{
}

void Perlin_noise::reset()
{

	std::random_device rd;
	std::mt19937 rgen(rd());
	std::uniform_real_distribution<> angle(0., 2.*PI);

	for (uint32_t i = 0; i < GRID_SIZE; ++i)
	{
		for (uint32_t j = 0; j < GRID_SIZE; ++j)
		{
			double theta = angle(rgen);
			m_grads[i][j] = { cos(theta), sin(theta) };
		}
	}
}

double Perlin_noise::fade(double t)
{
	return t*t*t*(10. + t*(6.*t - 15.));
}

double Perlin_noise::noise(double x, double y)
{
	double x0 = floor(x);
	double x1 = ceil(x);
	double y0 = floor(y);
	double y1 = ceil(y);

	x -= x0;
	y -= y0;

	double a = fade(1. - x)*glm::dot(m_grads[static_cast<size_t>(x0)][static_cast<size_t>(y0)], { x, y }) +
		fade(x)*glm::dot(m_grads[static_cast<size_t>(x1)][static_cast<size_t>(y0)], { x-1., y });
	double b= fade(1. - x)*glm::dot(m_grads[static_cast<size_t>(x0)][static_cast<size_t>(y1)], { x, y-1. }) +
		fade(x)*glm::dot(m_grads[static_cast<size_t>(x1)][static_cast<size_t>(y1)], { x-1., y-1. });

	return fade(1. - y)*a + fade(y)*b;
}

