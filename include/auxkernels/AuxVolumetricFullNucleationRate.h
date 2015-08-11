/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
* 13 January 2015
*
*************************************************************************/

#ifndef AUXVOLUMETRICFULLNUCLEATIONRATE_H
#define AUXVOLUMETRICFULLNUCLEATIONRATE_H

#include "AuxFullNucleationRate.h"

class AuxVolumetricFullNucleationRate;

template<>
InputParameters validParams<AuxVolumetricFullNucleationRate>();

class AuxVolumetricFullNucleationRate : public AuxFullNucleationRate
{
public:
  AuxVolumetricFullNucleationRate(const InputParameters & parameters);

protected:

//  virtual Real computeValue();

//  virtual void computeCriticalRadius();
//  virtual void computeCriticalEnergy();

//  virtual void computeZeldovichFactor();
//  virtual void computeCriticalFrequency();

  virtual void computeNumAtoms();

  /*
  const unsigned int _mesh_dimension;

  VariableValue & _coupled_energy;          // J/m^3 free energy change of transformation

  Real _Z;                                  // Zeldovich non-equilibrium factor
  Real _N;                                  // Number of atoms in the computational volume
  Real _beta_star;                          // frequency factor (critical -> supercritical nucleus)
  Real _linear_density;                     // linear atomic density of material


  Real _r_star;                             // Critical radius (2 or 3D)
  Real _G_star;                             // critical activation energy (J)


  Real _gamma;
  Real _Kb;
  //Real _temperature;
  Real _scale_factor;

  unsigned int _n_OP_variables;
  std::vector<VariableValue *> _OP;
  VariableValue & _T;
  VariableValue & _X;

  MaterialProperty<Real> & _D;

  Real _jump_distance;

  MaterialProperty<Real> & _Omega;
  Real _OP_threshold;
  Real _length_scale;
  */
private:
  Real _rate_volume;

};

#endif //AUXVOLUMETRICFULLNUCLEATIONRATE_H
