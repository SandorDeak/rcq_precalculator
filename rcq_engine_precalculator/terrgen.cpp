#include "terrgen.h"



terrgen::terrgen(const glm::uvec2& size): m_terrain(size), m_window(sf::VideoMode(size.x, size.y), "Debug info")
{
}


terrgen::~terrgen()
{
}

void terrgen::fill_with_noise(const glm::dvec2& scale)
{
	auto size = m_terrain.size();
	glm::dvec2 mult = scale / static_cast<glm::dvec2>(size);

	for (size_t i = 0; i < size.x; ++i)
	{
		for (size_t j = 0; j < size.y; ++j)
		{
			m_terrain[i][j].height = m_perlin.noise(static_cast<double>(i)*mult.x, static_cast<double>(j)*mult.y);
		}
	}
}

void terrgen::show_height()
{
	sf::RectangleShape r(sf::Vector2f(1, 1));

	glm::uvec2 size = m_terrain.size();
	for (uint32_t i = 0; i < size.x; ++i)
	{
		for (uint32_t j = 0; j < size.y; ++j)
		{		
			double val = (m_terrain[i][j].height+1.)*127.;
			r.setPosition(i, j);
			r.setFillColor(sf::Color(val, val, val, 255));
			m_window.draw(r);
		}
	}
	m_window.display();
}