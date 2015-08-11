/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#ifndef VALUE_H
#define VALUE_H

#include "Kernel.h"

//forward declarations
class Value;

/**
 * Value returns the variable value times the test function.   To be used
 * with UserForcingFunction for performing an L2 projection.
 */

template<>
InputParameters validParams<Value>();

class Value : public Kernel
{
public:

  Value(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

private:
};

#endif //VALUE_H
