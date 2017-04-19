//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  IPouplana $
//   Last modified by:    $Co-Author:             JMCarbonell $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MODIFIED_EXPONENTIAL_DAMAGE_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_MODIFIED_EXPONENTIAL_DAMAGE_HARDENING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_laws/hardening_law.hpp"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ModifiedExponentialDamageHardeningLaw
    : public HardeningLaw
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModifiedExponentialDamageHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION( ModifiedExponentialDamageHardeningLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModifiedExponentialDamageHardeningLaw();


    /// Copy constructor.
    ModifiedExponentialDamageHardeningLaw(ModifiedExponentialDamageHardeningLaw const& rOther);

    /// Assignment operator.
    ModifiedExponentialDamageHardeningLaw& operator=(ModifiedExponentialDamageHardeningLaw const& rOther);

    /// Clone.
    virtual HardeningLaw::Pointer Clone() const override;
    
    /// Destructor.
    ~ModifiedExponentialDamageHardeningLaw();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Hardening functions
     */

    virtual double& CalculateHardening(const PlasticDataType& rVariables, double &rHardening) override;
      
    /**
     * Calculate Hardening function derivatives
     */

    virtual double& CalculateDeltaHardening(const PlasticDataType& rVariables, double &rDeltaHardening) override;

    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    
    /// Turn back information as a string.
    virtual std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "ModifiedExponentialDamageHardeningLaw" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "ModifiedExponentialDamageHardeningLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "ModifiedExponentialDamageHardeningLaw Data";
    } 

    ///@}
    ///@name Friends
    ///@{


    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    
    ///@}
    ///@name Protected Operators
    ///@{

    
    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

  private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    
    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
    }
    
    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class ModifiedExponentialDamageHardeningLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// // input stream function

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODIFIED_EXPONENTIAL_DAMAGE_HARDENING_LAW_H_INCLUDED  defined 
