#include "sky_calculator.h"
#include "terrain_generator.h"
#include "terrgen.h"
#include <random>
#include <iostream>

#define TERRAIN

#ifdef SKY

int main()
{
	sky_calculator sky(32, 512);
	//sky.load_from_file("output/try.sky", 32);

	//sky.set_params(2.545e25, 5e25, 0.75);
	sky.set_params_directly({ 5.8e-6, 1.35e-5, 3.31e-5 }, 2.e-6, 0.75, glm::dvec3(3.426, 8.298, 0.356)*0.06e-5);
	//sky.set_params_directly({ 5.8e-6, 1.35e-5, 3.31e-5 }, 2.e-4, 0.3);
	sky.load_Lebedev_params("Lebedev_data/lebedev110_17.txt");
	sky.compute(0, "output/try");
  	return 0;
}

#endif //SKY

#ifdef  TERRAIN

int main()
{
	/*terrain_generator terrain({ 256, 256 });
	terrain.generate();*/

	terrgen t({ 512, 512 });
	t.fill_with_noise({ 31.9, 31.9 });
	t.show_height();


	return 0;
}

#endif //  TERRAIN
