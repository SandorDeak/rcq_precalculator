#include "sky_calculator.h"

int main()
{
	sky_calculator sky(256);
	sky.set_params(2.545e25, 1e24, 0.75);
	sky.load_Lebedev_params("Lebedev_data/lebedev110_17.txt");
	sky.compute(4);
	sky.write_to_file("output/try.sky");

	return 0;
}