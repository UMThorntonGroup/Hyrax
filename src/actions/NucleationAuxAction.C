/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 July 2013
*
*************************************************************************/

#include "NucleationAuxAction.h"
#include "Factory.h"
#include "FEProblem.h"

template<>
InputParameters validParams<NucleationAuxAction>()
{
  InputParameters params = validParams<Action>();

  //for the action
  params.addRequiredParam<unsigned int>("number_OPs", "specifies the number of order parameters");
  params.addRequiredParam<std::string>("OP_name_base", "name base for the OP variables");
  params.addRequiredParam<std::string>("bulk_energy_name_base", "name base for bulk energy change auxvar (i.e., Auxchemelastic");
  params.addRequiredParam<std::string>("nucleation_rate_name_base", "name base for nucleation rate auxvar");
  params.addRequiredParam<std::string>("nucleation_probability_name_base", "name base for nucleation probability auxvar");
  params.addRequiredParam<std::string>("deltaGstar_name_base", "name base for activation energy auxvar");

  //for bulk free energy - AuxChem or AuxChemElastic
  //params.addParam<bool>("use_auxchem", false, "false to use AuxChemElastic, true to use AuxChem");
  params.addRequiredParam<std::string>("bulk_energy_name", "name of which bulk energy auxkernel to use");
  params.addRequiredParam<std::string>("coupled_conserved_var", "coupled conserved field variable");
  params.addRequiredParam<Real>("precip_conserved", "value of the equilibrium 2nd phase conserved field variable");
  params.addRequiredParam<Real>("precip_nonconserved", "value of the equilibrium 2nd phase nonconserved field variable");

  //for NucleationRate
  params.addRequiredParam<Real>("gamma", "Surface energy");
  params.addParam<Real>("Kb", 1.3806503e-23, "Boltzmann's constant, make sure units all match");
  params.addParam<Real>("temperature", 473, "Temperature");
  params.addRequiredParam<Real>("scale_factor","factor to scale energy/dimensions by");
  params.addRequiredParam<Real>("Z", "Non-equilibrium Zeldovitch factor");
  params.addRequiredParam<Real>("Beta_star", "1/characteristic nucleation time");
  params.addRequiredParam<Real>("linear_density", "linear atomic density of matrix");

  //for NucleationProbability
  params.addParam<Real>("OP_threshold", 0.01, "Threshold value for determining existence within 2nd phase");

  return params;
}


NucleationAuxAction::NucleationAuxAction(InputParameters params) :
    Action(params),
    _num_OPs(getParam<unsigned int>("number_OPs")),
    _OP_name_base(getParam<std::string>("OP_name_base")),
    _bulk_energy_name_base(getParam<std::string>("bulk_energy_name_base")),
    _nucleation_rate_name_base(getParam<std::string>("nucleation_rate_name_base")),
    _nucleation_probability_name_base(getParam<std::string>("nucleation_probability_name_base")),
    _deltaGstar_name_base(getParam<std::string>("deltaGstar_name_base")),

    //_use_auxchem(getParam<bool>("use_auxchem")),
    _bulk_energy_name(getParam<std::string>("bulk_energy_name")),
    _coupled_conserved_var(getParam<std::string>("coupled_conserved_var")),
    _precip_conserved(getParam<Real>("precip_conserved")),
    _precip_nonconserved(getParam<Real>("precip_nonconserved")),

    _Z(getParam<Real>("Z")),
    _beta_star(getParam<Real>("Beta_star")),
    _linear_density(getParam<Real>("linear_density")),
    _gamma(getParam<Real>("gamma")),
    _Kb(getParam<Real>("Kb")),
    _temperature(getParam<Real>("temperature")),
    _scale_factor(getParam<Real>("scale_factor")),

    _OP_threshold(getParam<Real>("OP_threshold"))
{
}

void
NucleationAuxAction::act()
{
  //create vector of OP variable names
  std::vector<VariableName> OP_vector;
  OP_vector.resize(_num_OPs);

  for(unsigned int i = 1; i <= _num_OPs; ++i)
  {
    std::stringstream out;
    out << i;

    //create OP variable names
    std::string OP_var_name = _OP_name_base;
    OP_var_name.append(out.str());

    OP_vector[i-1] = OP_var_name;
  }

  //do all the auxkernel actions
  for (unsigned int i = 1; i <= _num_OPs; ++i)
  {
    //create auxvariable names
    std::string rate_name = _nucleation_rate_name_base;
    std::string probability_name = _nucleation_probability_name_base;
    std::string G_star_name = _deltaGstar_name_base;
    std::string bulk_energy_name = _bulk_energy_name_base;

    rate_name.append(OP_vector[i-1]);
    bulk_energy_name.append(OP_vector[i-1]);
    probability_name.append(OP_vector[i-1]);
    G_star_name.append(OP_vector[i-1]);

    std::string auxkernel_name;

    //Bulk Energy input parameters
    //get the parameters for AuxChem, CURRENTLY can use for AuxChem AND AuxChemElastic
    InputParameters action_params = _factory.getValidParams("AuxChem");

    action_params.set<AuxVariableName>("variable") = bulk_energy_name;
    action_params.set<std::vector<VariableName> >("coupled_nonconserved_var")
      = std::vector<VariableName> (1, OP_vector[i-1]);
    action_params.set<int>("nonconserved_var_number") = i;

    action_params.set<std::vector<VariableName> >("coupled_conserved_var")
      = std::vector<VariableName> (1, _coupled_conserved_var);
    action_params.set<Real>("precip_conserved") = _precip_conserved;
    action_params.set<Real>("precip_nonconserved") = _precip_nonconserved;
    action_params.set<int>("nonconserved_var_number") = i;

    //figure out if we need AuxChem or AuxChemElastic
    if(_bulk_energy_name == "AuxChem")
    {
      auxkernel_name = "AuxChem_";
      auxkernel_name.append(OP_vector[i-1]);
      _problem->addAuxKernel("AuxChem", auxkernel_name, action_params);
    }
    else if(_bulk_energy_name == "AuxGuoEnergy")
    {
      auxkernel_name = "AuxGuoEnergy_";
      auxkernel_name.append(OP_vector[i-1]);
      _problem->addAuxKernel("AuxGuoEnergy", auxkernel_name, action_params);
    }
    else if(_bulk_energy_name == "AuxCalphadEnergy")
    {
      auxkernel_name = "AuxCalphadEnergy_";
      auxkernel_name.append(OP_vector[i-1]);
      _problem->addAuxKernel("AuxCalphadEnergy", auxkernel_name, action_params);
    }
    else
      mooseError("Please enter AuxKernel bulk energy name for NucleationAuxAction");


    //get the parameters for nucleation rate
    action_params = _factory.getValidParams("AuxNucleationRate");
    action_params.set<AuxVariableName>("variable") = rate_name;
    action_params.set<Real>("gamma") = _gamma;
    action_params.set<Real>("Kb") = _Kb;
    action_params.set<Real>("temperature") = _temperature;
    action_params.set<Real>("scale_factor") = _scale_factor;
    action_params.set<Real>("Z") = _Z;
    action_params.set<Real>("Beta_star") = _beta_star;
    action_params.set<Real>("linear_density") = _linear_density;
    action_params.set<std::vector<VariableName> >("coupled_aux_var")
      = std::vector<VariableName> (1, bulk_energy_name);

    auxkernel_name = "NucleationRate_";
    auxkernel_name.append(OP_vector[i-1]);

    _problem->addAuxKernel("AuxNucleationRate", auxkernel_name, action_params);

   //get the parameters for nucleation probability
    action_params = _factory.getValidParams("AuxNucleationProbability");
    action_params.set<AuxVariableName>("variable") = probability_name;
    action_params.set<std::vector<VariableName> >("coupled_aux_var")
      = std::vector<VariableName> (1, rate_name);
    action_params.set<int>("n_OP_vars") = _num_OPs;
    action_params.set<std::vector<VariableName> >("coupled_variables") = OP_vector;
    action_params.set<Real>("OP_threshold") = _OP_threshold;

    auxkernel_name = "NucleationProbability_";
    auxkernel_name.append(OP_vector[i-1]);

    _problem->addAuxKernel("AuxNucleationProbability", auxkernel_name, action_params);

    //get the parameters for delta G*
    action_params = _factory.getValidParams("AuxDeltaGStar");
    action_params.set<AuxVariableName>("variable") = G_star_name;
    action_params.set<std::vector<VariableName> >("coupled_aux_var")
      = std::vector<VariableName> (1, bulk_energy_name);
    action_params.set<Real>("gamma") = _gamma;

    auxkernel_name = "DeltaGStar_";
    auxkernel_name.append(OP_vector[i-1]);

    _problem->addAuxKernel("AuxDeltaGStar", auxkernel_name, action_params);
  }
}
