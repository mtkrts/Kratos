// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined (KRATOS_SPRING_LAW_H_INCLUDED)
#define  KRATOS_SPRING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

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
 * @class PythonConstitutiveLawFunction
 * @ingroup StructuralMechanicsApplication
 * @brief This function allows to call a function method of type f(x, y, z, t) implemented in python.
 * @details Uses python functions to evaluate the bahaviour of the CL
 * The functions can be constructed by providing a python-defined method of the type
 *
 *  class aux_object_cpp_callback:
 *    def __init__(self, function_string ):
 *        #TODO: check python version
 *        self.compiled_function = compile(function_string, '', 'eval', optimize=2)
 *
 *    def f(self,x,y,z,t):
 *        return  eval(self.compiled_function)
 *
 * the object is then insantiated as
 * aux_function = PythonConstitutiveLawFunction(aux_object_cpp_callback(self.function_string))
 * @note Based on python_function_callback_utility
 * @note This makes this file to depend on python
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class PythonConstitutiveLawFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// The node defintion
    typedef Node<3> NodeType;

    /// Counted pointer of PythonConstitutiveLawFunction
    KRATOS_CLASS_POINTER_DEFINITION(PythonConstitutiveLawFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rFunctionBody The text containing the function
     */
    PythonConstitutiveLawFunction( const std::string& rFunctionBody)
    {
        // Compile the function starting from the string function body
        try {
            mMainModule = pybind11::module::import("__main__");
            mMainNameSpace = mMainModule.attr("__dict__");
            pybind11::exec("from math import *", mMainNameSpace);
            mFunctionBody = rFunctionBody;
        } catch(pybind11::error_already_set const&) {
            PyErr_Print();
        }
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This methods returns the value of the function
     * @param ThisNode The node whre to evaluate
     * @param Time The current time
     * @return The resulting value of the function
     */
    double CallFunction(
        NodeType& ThisNode,
        const double Time
        )
    {
        /// Position
        mMainNameSpace["x"] = ThisNode.X();
        mMainNameSpace["y"] = ThisNode.Y();
        mMainNameSpace["z"] = ThisNode.Z();
        mMainNameSpace["X"] = ThisNode.X0();
        mMainNameSpace["Y"] = ThisNode.Y0();
        mMainNameSpace["Z"] = ThisNode.Z0();

        /// Values
        const array_1d<double, 3>& disp = ThisNode.FastGetSolutionStepValue(DISPLACEMENT);
        mMainNameSpace["disp0"] = disp[0];
        mMainNameSpace["disp1"] = disp[1];
        mMainNameSpace["disp2"] = disp[2];
        const array_1d<double, 3>& vel = ThisNode.FastGetSolutionStepValue(VELOCITY);
        mMainNameSpace["vel0"] = vel[0];
        mMainNameSpace["vel1"] = vel[1];
        mMainNameSpace["vel2"] = vel[2];
        const array_1d<double, 3>& accel = ThisNode.FastGetSolutionStepValue(ACCELERATION);
        mMainNameSpace["accel0"] = accel[0];
        mMainNameSpace["accel1"] = accel[1];
        mMainNameSpace["accel2"] = accel[2];
        const array_1d<double, 3>& theta = ThisNode.FastGetSolutionStepValue(ROTATION);
        mMainNameSpace["theta0"] = theta[0];
        mMainNameSpace["theta1"] = theta[1];
        mMainNameSpace["theta2"] = theta[2];
        const array_1d<double, 3>& vang = ThisNode.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        mMainNameSpace["vang0"] = vang[0];
        mMainNameSpace["vang1"] = vang[1];
        mMainNameSpace["vang2"] = vang[2];
        const array_1d<double, 3>& aang = ThisNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION);
        mMainNameSpace["aang0"] = aang[0];
        mMainNameSpace["aang1"] = aang[1];
        mMainNameSpace["aang2"] = aang[2];

        /// Increment of values
        const array_1d<double, 3>& dispDelta = ThisNode.FastGetSolutionStepValue(DISPLACEMENT) - ThisNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
        mMainNameSpace["dispDelta0"] = dispDelta[0];
        mMainNameSpace["dispDelta1"] = dispDelta[1];
        mMainNameSpace["dispDelta2"] = dispDelta[2];
        const array_1d<double, 3>& velDelta = ThisNode.FastGetSolutionStepValue(VELOCITY) - ThisNode.FastGetSolutionStepValue(VELOCITY, 1);
        mMainNameSpace["velDelta0"] = velDelta[0];
        mMainNameSpace["velDelta1"] = velDelta[1];
        mMainNameSpace["velDelta2"] = velDelta[2];
        const array_1d<double, 3>& accelDelta = ThisNode.FastGetSolutionStepValue(ACCELERATION) - ThisNode.FastGetSolutionStepValue(ACCELERATION, 1);
        mMainNameSpace["accelDelta0"] = accelDelta[0];
        mMainNameSpace["accelDelta1"] = accelDelta[1];
        mMainNameSpace["accelDelta2"] = accelDelta[2];
        const array_1d<double, 3>& thetaDelta = ThisNode.FastGetSolutionStepValue(ROTATION) - ThisNode.FastGetSolutionStepValue(ROTATION, 1);
        mMainNameSpace["thetaDelta0"] = thetaDelta[0];
        mMainNameSpace["thetaDelta1"] = thetaDelta[1];
        mMainNameSpace["thetaDelta2"] = thetaDelta[2];
        const array_1d<double, 3>& vangDelta = ThisNode.FastGetSolutionStepValue(ANGULAR_VELOCITY) - ThisNode.FastGetSolutionStepValue(ANGULAR_VELOCITY, 1);
        mMainNameSpace["vangDelta0"] = vangDelta[0];
        mMainNameSpace["vangDelta1"] = vangDelta[1];
        mMainNameSpace["vangDelta2"] = vangDelta[2];
        const array_1d<double, 3>& aangDelta = ThisNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION) - ThisNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 1);
        mMainNameSpace["aangDelta0"] = aangDelta[0];
        mMainNameSpace["aangDelta1"] = aangDelta[1];
        mMainNameSpace["aangDelta2"] = aangDelta[2];

        /// Time
        mMainNameSpace["t"] = Time;

        return pybind11::eval(mFunctionBody, mMainNameSpace).cast<double>();
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    pybind11::object mMainModule;    /// The python main module
    pybind11::object mMainNameSpace; /// The variables that generate the dependence of the function
    std::string mFunctionBody;       /// The text definting the function

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

    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
};

/**
 * @class SpringConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief Spring/damper/mass/inertia... constitutive law for 3D and 2D points
 * @details Uses python functions to evaluate the bahaviour of the CL
 * @tparam TDim The size of the space
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SpringConstitutiveLaw
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// The process info type
    typedef ProcessInfo       ProcessInfoType;

    /// The constitutive law definiton
    typedef ConstitutiveLaw          BaseType;

    /// The size definition
    typedef std::size_t              SizeType;

    /// The index type
    typedef std::size_t             IndexType;

    /// Counted pointer of SpringConstitutiveLaw
    KRATOS_CLASS_POINTER_DEFINITION( SpringConstitutiveLaw );

    /**
     * @brief Flags related to constitutive law computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( NULL_MASS );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_INERTIA_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_INERTIA_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_INERTIA_Z );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_STIFFNESS_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_STIFFNESS_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_STIFFNESS_Z );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATTIONAL_STIFFNESS_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATTIONAL_STIFFNESS_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATTIONAL_STIFFNESS_Z );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SpringConstitutiveLaw();

    /**
     * @brief Default constructor. Using parameters instead of the one by default
     * @param NewParameters The configuration parameters of the new constitutive law
     */
    SpringConstitutiveLaw(Kratos::Parameters NewParameters);

    /**
     * @brief Copy constructor.
     */
    SpringConstitutiveLaw (const SpringConstitutiveLaw& rOther);


    /**
     * @brief Destructor.
     */
    ~SpringConstitutiveLaw() override;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * @note implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
     *      return p_clone;
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief Creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    Pointer Create(Kratos::Parameters NewParameters) override
    {
        return Kratos::make_shared<SpringConstitutiveLaw<TDim>>(NewParameters);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 1;
    }
    
    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (array of 3 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * @note Fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    bool Has(const Variable<array_1d<double, 3 > >& rThisVariable) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (array of 3 components)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    array_1d<double, 3 > & CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue
        ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

protected:

    ///@name Protected static Member Variables
    ///@{
    
    ///@}
    ///@name Protected member Variables
    ///@{
    
    /* The definition of the functions */ // TODO: Add to serializer
    PythonConstitutiveLawFunction::Pointer mMassFuntion = nullptr;                      /// The function of the mass
    array_1d<PythonConstitutiveLawFunction::Pointer, TDim> mInertiaFunction;            /// The function of the inertia
    array_1d<PythonConstitutiveLawFunction::Pointer, TDim> mStiffnessunction;           /// The function of the stiffness
    array_1d<PythonConstitutiveLawFunction::Pointer, TDim> mRotationalStiffnessunction; /// The function of the rotational stiffness

    Flags mConstitutiveLawFlags;       /// Constitutive flags

    array_1d<double, 2> mTimeInterval; /// The time interval

    ///@}
    ///@name Protected Operators
    ///@{
    
    ///@}
    ///@name Protected Operations
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

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.save("ConstitutiveLawFlags",mConstitutiveLawFlags);
        rSerializer.save("TimeInterval",mTimeInterval);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("ConstitutiveLawFlags",mConstitutiveLawFlags);
        rSerializer.load("TimeInterval",mTimeInterval);
    }


}; // Class SpringConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_SPRING_LAW_H_INCLUDED  defined
