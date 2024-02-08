#include <iostream>
#include <sstream>
#include <sys/stat.h>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../src/actionData.h"
#include "../src/speGrid.h"
/*
 * Created by Andrew S Maxwell 12/12/2023
 *
 * NSDI code using a sin square pulse envelope
 */

int main(int argc, char **argv) {
    speGridParameters gridParam; //Trivial code so tests, work//TODO add catch2 with main to CMake
	// --- Testing suite --- //
	int result = 0;
	//uncomment below to run all unit tests
	result = Catch::Session().run( argc, argv );
	// --- Testing End --- //

	
	return result;
}

TEST_CASE( "1: All test cases reside in other .cpp files (empty)", "[multi-file:1]" ) {
}
