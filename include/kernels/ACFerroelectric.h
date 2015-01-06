/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 August 2013
*
*************************************************************************/

#ifndef ACFERROELECTRIC_H
#define ACFERROELECTRIC_H

#include "ACBulk.h"

//forward declarations
class ACFerroelectric;

template<>
InputParameters validParams<ACFerroelectric>();

class ACFerroelectric : public ACBulk
{
public:
  ACFerroelectric(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

  MaterialProperty<Real> & _a1;
  MaterialProperty<Real> & _a11;
  MaterialProperty<Real> & _a12;

  unsigned int _n_OP_vars;
  unsigned int _OP_number;

  std::vector<VariableValue *> _coupled_OP_vars;

private:

};


#endif //ACFERROELECTRIC_H
