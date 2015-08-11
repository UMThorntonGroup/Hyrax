/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 July 2013
*
*************************************************************************/

#ifndef OPVARIANTKERNELACTION_H
#define OPVARIANTKERNELACTION_H

#include "InputParameters.h"
#include "Action.h"

//forward declaration
class OPVariantKernelAction;

template<>
InputParameters validParams<OPVariantKernelAction>();


class OPVariantKernelAction : public Action
{
public:

  OPVariantKernelAction(InputParameters params);

  virtual void act();

protected:

private:

  unsigned int _num_OPs;
  std::string _OP_name_base;
  std::string _kappa_name_OP;
  std::string _coupled_CH_var;

  bool _use_elasticity;
};

#endif //OPVARIANTKERNELACTION_H
