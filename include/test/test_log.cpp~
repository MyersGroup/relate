#include "catch.hpp"

#include "fast_log.hpp"

TEST_CASE( "Testing fast log" ){
  float tolerance = 0.007;

  //check accuracy over range of 1e-3 to 1e3
  float test_num = 1e-5;
  for(int i = 1; i < 10; i++){
    test_num *= 10;
    REQUIRE( std::fabs( fast_log(test_num) - std::log(test_num) ) < tolerance );
  }
}


