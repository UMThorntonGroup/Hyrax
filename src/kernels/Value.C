/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#include "Value.h"

template<>
InputParameters validParams<Value>()
{
  InputParameters params = validParams<Kernel>();

  return params;
}

Value::Value(const InputParameters & parameters)
    : Kernel(parameters)
{
}

Real
Value::computeQpResidual()
{
  return _u[_qp]*_test[_i][_qp];
}
