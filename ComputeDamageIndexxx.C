// https://github.com/idaholab/blackbear/blob/devel/src/materials/MazarsDamage.C

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                       BlackBear                              */
/*                                                              */
/*           (c) 2017 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ComputeDamageIndexxx.h"




registerMooseObject("first_appApp", ComputeDamageIndexxx);
registerMooseObject("first_appApp", ADComputeDamageIndexxx);

template <bool is_ad>
InputParameters
ComputeDamageIndexxxTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute Damage Evolution");

  params.addRequiredParam<MaterialPropertyName>(
      "Yin", "Name of material property defining the Threshold parameter controlling the initiation of failure");
  params.addRequiredParam<MaterialPropertyName>(
      "P1", "Name of material property that control how quickly failure occurs after initiating");  
  params.addRequiredParam<MaterialPropertyName>(
      "P2", "Name of material property that control the general shape of the stress response after failure");
  params.addRequiredParam<MaterialPropertyName>(
      "Mu", "Name of material property defining the Damage consistency");
  
 
 params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addParam<std::string>("stress_name",
                               "stress",
                               "Optional parameter that allows the user to use "
                               "different stresses on the same material system. "
                               "For example, when we have a degraded_stress and an intact_stress, "
                               "we want to compute the degraded strain energy density and "
                               "the intact strain energy density.");
 
  return params;
}

template <bool is_ad>
ComputeDamageIndexxxTempl<is_ad>::ComputeDamageIndexxxTempl(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),

  // Damage law parameters
  _Yin(getMaterialProperty<Real>("Yin")),
  _P1(getMaterialProperty<Real>("P1")),
  _P2(getMaterialProperty<Real>("P2")),
  _Mu(getMaterialProperty<Real>("Mu")),
  
 

  // strain energy density 
 
  _damage_indexx(declareProperty<Real>("damage_indexx")),
  _damage_indexx_old(getMaterialPropertyOld<Real>("damage_indexx")),
  
  
        
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress_name(getParam<std::string>("stress_name")),
   
   /// _stress(getGenericMaterialProperty<RankTwoTensor, is_ad>(_base_name + _stress_name)),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + _stress_name)),
    _mechanical_strain(
        getGenericMaterialProperty<RankTwoTensor, is_ad>(_base_name + "mechanical_strain")),
    _strain_increment(
        getGenericOptionalMaterialProperty<RankTwoTensor, is_ad>(_base_name + "strain_increment")),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name))  ,
    _stresss(declareProperty<RankTwoTensor>(_base_name + "stresss")),
  
   _dt(_fe_problem.dt())
  
   
 
{
}



template <bool is_ad>
void
ComputeDamageIndexxxTempl<is_ad>::initQpStatefulProperties()
{
  
   _damage_indexx[_qp] = 0.0;
  
}



template <bool is_ad>
void
ComputeDamageIndexxxTempl<is_ad>::computeQpProperties()
{
  unsigned int i = 2;
  unsigned int j = 2;
 // stress = C * e
 std::cout << "_eleastic_tensor" << std::endl;
 _elasticity_tensor[_qp].print();
 
  // undamages stress
  
   _stresss[_qp] = MetaPhysicL::raw_value(_elasticity_tensor[_qp]).contractionIj(i,j,(MetaPhysicL::raw_value(_mechanical_strain[_qp]))); 
   
   std::cout << "_stresss" << std::endl;
   _stresss[_qp].print();
   
   std::cout << "_mechanical_strain" << std::endl;
  _mechanical_strain[_qp].print();
   
   Real strain_energy_density = MetaPhysicL::raw_value(_mechanical_strain[_qp])
           .doubleContraction(MetaPhysicL::raw_value(_stresss[_qp])) /
       2.0;
  
 
 
  // access damage at previous timestep
  Real d_old = _damage_indexx_old[_qp];
 
 
  // compute val_1
  Real val_1 = - (std::pow((((strain_energy_density / (1 - d_old))- _Yin[_qp])/(_P1[_qp] * _Yin[_qp])), _P2[_qp]));

  // compute G
  Real G_value = (1 - std::exp(val_1));
  
 // compute g 
  Real g_value = (G_value - _damage_indexx_old[_qp]);
  
  Real dt = _dt;
  
 
 
  //Calculation of damage evolution equation 
 if (g_value >= 0)
            
     _damage_indexx[_qp] = ((d_old + (_Mu[_qp]*dt*G_value))/ (1 + (_Mu[_qp]*dt)));
     
 else
    _damage_indexx[_qp] = d_old;
 
 

}

template class ComputeDamageIndexxxTempl<false>;
template class ComputeDamageIndexxxTempl<true>;
