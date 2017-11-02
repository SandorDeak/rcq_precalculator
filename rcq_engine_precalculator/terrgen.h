#pragma once

#include <glm.hpp>
#include <vector>
#include <SFML\Graphics.hpp>
#include "Perlin_noise.h"

struct cell
{
	double height;
};

class terrain
{
public:
	terrain(const glm::uvec2& size) : m_data(size.x)
	{
		for (auto& v : m_data)
			v.resize(size.y);
		m_size = size;
	}
	glm::uvec2 size()
	{
		return m_size;
	}
	glm::uvec2 center()
	{
		return m_size / 2u;
	}
	std::vector<cell>& operator[](size_t i)
	{
		return m_data[i];
	}
private:
	std::vector<std::vector<cell>> m_data;
	glm::uvec2 m_size;
};


class terrgen
{
public:
	terrgen(const glm::uvec2& size);
	~terrgen();

	void fill_with_noise(const glm::dvec2& scale );
	void show_height();


private:
	terrain m_terrain;

	Perlin_noise m_perlin;

	sf::RenderWindow m_window;
};

