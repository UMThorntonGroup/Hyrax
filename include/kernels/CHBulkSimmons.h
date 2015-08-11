/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 October 2012
*
*************************************************************************/

#ifndef CHBULKSIMMONS_H
#define CHBULKSIMMONS_H

#include "CHBulkCoupled.h"

/** CHBulkSimmons handles the conserved order parameter(probably concentration),
 * evolved using the Cahn-Hilliard equation.  It couples to an order
 * parameter from the Alan-Cahn equation.  It uses the PFMobilityLandau materials class.
 * This is designed for the free energy functional in Simmons 2001 concurrent nucleation
 * and growth paper.
 */

//Forward Declarations
class CHBulkSimmons;

template<>
InputParameters validParams<CHBulkSimmons>();

class CHBulkSimmons : public CHBulkCoupled
{
public:

  CHBulkSimmons(const InputParameters & parameters);

protected:

  /**
   * computeGradDFDCons()
   * @return returns the GRADIENT of the partial(bulk free energy)/partial(c).  Don't screw that up.
   */

  virtual RealGradient computeGradDFDCons(PFFunctionType type);

private:

};

#endif //CHBULKSIMMONS_H
