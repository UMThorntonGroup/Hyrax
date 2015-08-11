/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 October 2012
*
*************************************************************************/

#include "ACBulkSimmons.h"

template<>
InputParameters validParams<ACBulkSimmons>()
{
  InputParameters params = validParams<ACBulkCoupled>();
  return params;
}

ACBulkSimmons::ACBulkSimmons(const InputParameters & parameters)
  :ACBulkCoupled(parameters)
{
}

Real
ACBulkSimmons::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
  case Residual:
    return _a2[_qp]*(_c2[_qp] - _coupled_CH_var[_qp])*_u[_qp] - _a3[_qp]*_u[_qp]*_u[_qp]*_u[_qp]
           + _a4[_qp]*_u[_qp]*_u[_qp]*_u[_qp]*_u[_qp]*_u[_qp] ;

  case Jacobian:
    return _phi[_j][_qp]*( _a2[_qp]*(_c2[_qp] - _coupled_CH_var[_qp]) - 3.0*_a3[_qp]*_u[_qp]*_u[_qp] +
           5.0*_a4[_qp]*_u[_qp]*_u[_qp]*_u[_qp]*_u[_qp] ) ;
  }

  mooseError("Invalid type passed in");
}
