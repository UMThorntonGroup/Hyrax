/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  19 December 2012
*
*  ErrorFractionMaxHMarker is a clone of ErrorFractionMarker with the
*  added bonus of setting a maximum refinement level.
*
*************************************************************************/

#include "ErrorFractionMaxHMarker.h"

template<>
InputParameters validParams<ErrorFractionMaxHMarker>()
{
  InputParameters params = validParams<ErrorFractionMarker>();

  params.addRequiredParam<unsigned int>("max_h_level", "Maximum h-adapt level for the mesh");

  return params;
}


ErrorFractionMaxHMarker::ErrorFractionMaxHMarker(const InputParameters & parameters) :
    ErrorFractionMarker(parameters),
    _max_h_level(getParam<unsigned int>("max_h_level"))
{
}

Marker::MarkerValue
ErrorFractionMaxHMarker::computeElementMarker()
{
  Real error = _error_vector[_current_elem->id()];

  if((error > _refine_cutoff) && (_current_elem->level() < _max_h_level))
      return REFINE;
    else if(error < _coarsen_cutoff)
      return COARSEN;

  return DO_NOTHING;
}

