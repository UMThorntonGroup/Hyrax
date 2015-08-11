/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 July 2013
*
*************************************************************************/

#include "OPVariantKernelAction.h"
#include "Factory.h"
#include "FEProblem.h"

template<>
InputParameters validParams<OPVariantKernelAction>()
{
  InputParameters params = validParams<Action>();

  params.addRequiredParam<unsigned int>("number_OPs", "specifies the number of order parameters");
  params.addRequiredParam<std::string>("OP_name_base", "name base of the OP variable");
  params.addParam<std::string>("kappa_name_OP", "kappa_n", "the OP kappa used with ACInterface");
  params.addRequiredParam<std::string>("coupled_CH_var", "coupled conserved variable for ACBulkPolyCoupled");
  params.addParam<bool>("use_elasticity", true, "true if using ACTransformElasticDF, false if leaving it out");

  return params;
}

OPVariantKernelAction::OPVariantKernelAction(InputParameters params) :
    Action(params),
    _num_OPs(getParam<unsigned int>("number_OPs")),
    _OP_name_base(getParam<std::string>("OP_name_base")),
    _kappa_name_OP(getParam<std::string>("kappa_name_OP")),
    _coupled_CH_var(getParam<std::string>("coupled_CH_var")),
    _use_elasticity(getParam<bool>("use_elasticity"))
{
}

void
OPVariantKernelAction::act()
{
  //create vector of OP names
  std::vector<VariableName> OP_vector;
  OP_vector.resize(_num_OPs);

  for(unsigned int i = 1; i <= _num_OPs; ++i)
  {
    //create variable names
    std::string variable_name = _OP_name_base;
    std::stringstream out;
    out << i;
    variable_name.append(out.str());

    OP_vector[i-1] = variable_name;
  }

  // do all the kernel actions
  for(unsigned int i = 1; i <= _num_OPs; ++i)
  {
    std::string kernel_name;

    //Get the parameters for the ACBulkPolyCoupled kernel and add kernels to problem
    InputParameters action_params = _factory.getValidParams("ACBulkPolyCoupled");
    action_params.set<NonlinearVariableName>("variable") = OP_vector[i-1];
    action_params.set<std::vector<VariableName> >("OP_var_names") = OP_vector;
    action_params.set<int>("n_OP_vars") = _num_OPs;
    action_params.set<int>("OP_number") = i;

    action_params.set<std::vector<VariableName> >("coupled_CH_var")
      = std::vector<VariableName> (1, _coupled_CH_var);

    kernel_name = "ACBulkPolyCoupled_";
    kernel_name.append(OP_vector[i-1]);
    _problem->addKernel("ACBulkPolyCoupled", kernel_name, action_params);

    if(_use_elasticity)
    {
      //Get the parameters for the ACTransformElasticDF kernel and add kernels to problem
      action_params = _factory.getValidParams("ACTransformElasticDF");
      action_params.set<NonlinearVariableName>("variable") = OP_vector[i-1];
      action_params.set<std::vector<VariableName> >("OP_var_names") = OP_vector;
      action_params.set<int>("n_OP_vars") = _num_OPs;
      action_params.set<int>("OP_number") = i;

      kernel_name = "ACTransformElasticDF_";
      kernel_name.append(OP_vector[i-1]);
      _problem->addKernel("ACTransformElasticDF", kernel_name, action_params);
    }

    //Get the parameters for the ACInterface kernel and add kernels to problem
    action_params = _factory.getValidParams("ACInterface");
    action_params.set<NonlinearVariableName>("variable") = OP_vector[i-1];
    action_params.set<std::string>("kappa_name") = _kappa_name_OP;

    kernel_name = "ACInterface_";
    kernel_name.append(OP_vector[i-1]);
    _problem->addKernel("ACInterface", kernel_name, action_params);

    //Get the parameters for the TimeDerivative kernel and add kernels to problem
    action_params = _factory.getValidParams("TimeDerivative");
    action_params.set<NonlinearVariableName>("variable") = OP_vector[i-1];

    kernel_name = "TimeDeriv_";
    kernel_name.append(OP_vector[i-1]);
    _problem->addKernel("TimeDerivative", kernel_name, action_params);


    //anything else you need to add?
  }
}
