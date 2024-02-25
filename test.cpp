/***********************************************************************
 * Source File:
 *    Test : Test runner
 * Author:
 *    Br. Helfrich
 * Summary:
 *    The test runner for all the unit tests
 ************************************************************************/

#include "test.h"
#include "testPosition.h"
#include "testGround.h"

/*****************************************************************
 * TEST RUNNER
 * Runs all the unit tests
 ****************************************************************/
void testRunner()
{
   TestPosition().run();
   TestGround().run();
}
