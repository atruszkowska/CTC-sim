#include "../../include/flow/flowfield.h"
#include "../../include/utils.h"
#include "../../include/test_utils.h"

bool no_field();

int main()
{
	test_pass(no_field(), "Zero velocity field");
}

bool no_field()
{
	const size_t Nx = 10, Ny = 20;
	FlowField flow(Nx, Ny);
	return true;
}
