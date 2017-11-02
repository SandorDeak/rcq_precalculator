#pragma once

#include <glm.hpp>
#include <array>

constexpr size_t GRID_SIZE = 256;

class Perlin_noise
{
public:
	Perlin_noise();
	~Perlin_noise();

	void reset();
	double noise(double x, double y);

private:
	double fade(double t);

	std::array<std::array<glm::vec2, GRID_SIZE>, GRID_SIZE> m_grads;
};

