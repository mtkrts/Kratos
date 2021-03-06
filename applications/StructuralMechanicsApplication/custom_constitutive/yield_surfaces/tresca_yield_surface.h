// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//

#if !defined(KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED)
#define KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED

// System includes

// Project includes
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"

namespace Kratos
{
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
/**
 * @class TrescaYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TrescaYieldSurface
{
  public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of TrescaYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(TrescaYieldSurface);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    TrescaYieldSurface()
    {
    }

    /// Copy constructor
    TrescaYieldSurface(TrescaYieldSurface const &rOther)
    {
    }

    /// Assignment operator
    TrescaYieldSurface &operator=(TrescaYieldSurface const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~TrescaYieldSurface(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param StressVector The stress vector
     * @param StrainVector The StrainVector vector
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStress(
        const Vector &StressVector,
        const Vector &StrainVector,
        double &rEqStress,
        const Properties &rMaterialProperties)
    {
        double I1, J2, J3, lode_angle;
        Vector deviator = ZeroVector(6);

        ConstitutiveLawUtilities::CalculateI1Invariant(StressVector, I1);
        ConstitutiveLawUtilities::CalculateJ2Invariant(StressVector, I1, deviator, J2);
        ConstitutiveLawUtilities::CalculateJ3Invariant(deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, lode_angle);

        rEqStress = 2.0 * std::cos(lode_angle) * std::sqrt(J2);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rMaterialProperties The material properties
     */
    static void GetInitialUniaxialThreshold(const Properties &rMaterialProperties, double &rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_TENSION]); // TODO Check
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param AParameter The damage parameter
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameter(
        const Properties &rMaterialProperties,
        double &AParameter,
        const double CharacteristicLength)
    {
        const double Gf = rMaterialProperties[FRACTURE_ENERGY];
        const double E = rMaterialProperties[YOUNG_MODULUS];
        const double sigma_c = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double sigma_t = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = sigma_c / sigma_t;

        if (rMaterialProperties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential))
        {
            AParameter = 1.00 / (Gf * n * n * E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
            KRATOS_ERROR_IF(AParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
        }
        else
        { // linear
            AParameter = -std::pow(sigma_c, 2) / (2.0 * E * Gf * n * n / CharacteristicLength);
        }
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param StressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rg The derivative of the plastic potential
     * @param rMaterialProperties The material properties
     */
    static void CalculatePlasticPotentialDerivative(
        const Vector &StressVector,
        const Vector &Deviator,
        const double J2,
        Vector &rg,
        const Properties &rMaterialProperties)
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, rg, rMaterialProperties);
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param StressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rFFlux The derivative of the yield surface
     * @param rMaterialProperties The material properties
     */
    static void CalculateYieldSurfaceDerivative(
        const Vector &StressVector,
        const Vector &Deviator,
        const double J2,
        Vector &rFFlux,
        const Properties &rMaterialProperties)
    {
        Vector first_vector, second_vector, third_vector;

        ConstitutiveLawUtilities::CalculateFirstVector(first_vector);
        ConstitutiveLawUtilities::CalculateSecondVector(Deviator, J2, second_vector);
        ConstitutiveLawUtilities::CalculateThirdVector(Deviator, J2, third_vector);

        double J3, lode_angle;
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, lode_angle);

        const double checker = std::abs(lode_angle * 180 / Globals::Pi);

        double c1, c2, c3;
        c1 = 0.0;

        if (checker < 29.0)
        {
            c2 = 2.0 * (std::cos(lode_angle) + std::sin(lode_angle) * std::tan(3.0 * lode_angle));
            c3 = std::sqrt(3.0) * std::sin(lode_angle) / (J2 * std::cos(3.0 * lode_angle));
        }
        else
        {
            c2 = std::sqrt(3.0);
            c3 = 0.0;
        }

        noalias(rFFlux) = c1 * first_vector + c2 * second_vector + c3 * third_vector;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const
    {
    }

    void load(Serializer &rSerializer)
    {
    }

    ///@}

}; // Class TrescaYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
