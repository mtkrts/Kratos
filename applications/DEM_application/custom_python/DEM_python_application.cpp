/*
==============================================================================
KratosDEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author:  $
//   Date:                $Date: $
//   Revision:            $Revision: 1.3 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "DEM_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosDEMApplication)
{

    class_<KratosDEMApplication,
           KratosDEMApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosDEMApplication")
           ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();

    //registering variables in python

    // OPTIONS AND FLAGS

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(BOUNDING_BOX_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ROTATION_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CRITICAL_TIME_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(VIRTUAL_MASS_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(UNIFORM_MATERIAL_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ACTIVATE_SEARCH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LIMIT_SURFACE_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CLEAN_INDENT_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_VARIABLES_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TRIHEDRON_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTINUUM_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FORCE_CALCULATION_TYPE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DAMP_TYPE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ROTA_DAMP_TYPE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGH_INITIALIZED)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GROUP_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EXPORT_GROUP_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EXPORT_ID)

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DELTA_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CASE_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SKIN_SPHERE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ROTATION_SPRING_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NON_LINEAR_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_MESH_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FAILURE_CRITERION_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCRETE_TEST_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_CONTINUUM)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_COHESION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_TENSION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(plot_OPTIONS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SURFACE_FRICC)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(IF_BOUNDARY_ELEMENT)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(IF_BOUNDARY_FACE)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_CONTACT_FAILURE_ID)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_FAILURE_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EXPORT_PARTICLE_FAILURE_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PREDEFINED_SKIN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EXPORT_SKIN_SPHERE)
    // *************** Continuum only END ***************

    // *************** Swimming only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NON_NEWTONIAN_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MANUALLY_IMPOSED_DRAG_LAW_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DRAG_MODIFIER_TYPE)
    // *************** Swimming only END ***************

    // COUNTING

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NUM_OF_NEIGH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TOTAL_CONTACTS)
    // KRATOS_REGISTER_IN_PYTHON_VARIABLE(POTENTIAL_NEIGHBOURS)

    // TOLERANCES

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DISTANCE_TOLERANCE)

    // *************** Swimming only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(POWER_LAW_TOLERANCE)
    // *************** Swimming only END ***************

    // TIME RELATED

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TIME_INCREASING_RATIO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FINAL_SIMULATION_TIME)

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INITIAL_PRESSURE_TIME)
    // *************** Continuum only END ***************

    // MATERIAL PARAMETERS

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_MASS_COEFF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_MAGIC_FACTOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_MOMENT_OF_INERTIA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ROLLING_FRICTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_DENSITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_YOUNG_MODULUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_POISSON_RATIO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_COHESION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_ROLLING_FRICTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_STATIC_FRICTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_DYNAMIC_FRICTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GENERAL_RESTITUTION_COEFF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(HISTORICAL_MIN_K)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_KN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_KT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_INERTIA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_DENSITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_FRICTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_STATIC_FRICTION_COEF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_DYNAMIC_FRICTION_COEF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RESTITUTION_COEFF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_ROTATION_DAMP_RATIO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SQRT_OF_MASS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LN_OF_RESTITUTION_COEFF)

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_KR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_RN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_RT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_RR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GLOBAL_FRI_ANG)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NON_ELASTIC_STAGE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SLOPE_FRACTION_N1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SLOPE_FRACTION_N2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SLOPE_LIMIT_COEFF_C1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SLOPE_LIMIT_COEFF_C2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_FAILURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_SIGMA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_TAU)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FAILURE_CRITERION_STATE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_SIGMA_MAX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_SIGMA_MIN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_TAU_ZERO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTACT_INTERNAL_FRICC)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LOW_POISSON_FORCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(HIGH_POISSON_FORCE)
    // *************** Continuum only END *************

    // *************** Swimming only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GEL_STRENGTH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(POWER_LAW_N)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(POWER_LAW_K)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INIT_DRAG_FORCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DRAG_LAW_SLOPE)
    // *************** Swimming only END ***************

    // GEOMETRIC PARAMETERS

    // *************** Swimming only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SPHERICITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CYLINDER_RADIUS)
    // *************** Swimming only END ***************

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MAX_INDENTATION)

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LOCAL_CONTACT_AREA_HIGH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LOCAL_CONTACT_AREA_LOW)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MEAN_CONTACT_AREA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(REPRESENTATIVE_VOLUME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_INITIAL_DELTA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_CONTACT_DELTA)

    // *************** Continuum only END ***************

    // KINEMATICS

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_ROTATION_ANGLE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(EULER_ANGLES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ORIENTATION_REAL)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ORIENTATION_IMAG)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(INITIAL_AXES_TRACKING)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL_DIR)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SURFACE_POINT_COOR)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DELTA_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DELTA_ROTA_DISPLACEMENT)

    // FORCE AND MOMENTUM

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(TOTAL_FORCES)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_SURFACE_CONTACT_FORCES)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_MOMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MAX_ROTA_MOMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(APPLIED_FORCE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ANGULAR_MOMENTUM)

    // ENERGY

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(KINETIC_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELASTIC_ENERGY_OF_CONTACTS)

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(INITIAL_ROTA_MOMENT)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_BLOCK_CONTACT_FORCE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_SURFACE_ROTATE_SPRING_MOMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(EXTERNAL_APPLIED_FORCE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(LOCAL_CONTACT_FORCE_LOW)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(LOCAL_CONTACT_FORCE_HIGH)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_CONTACT_FORCES )

    // CONCRETE TEST

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FIXED_VEL_TOP)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FIXED_VEL_BOT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIAL_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AREA_VERTICAL_TAPA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AREA_VERTICAL_CENTRE)

    // TENSION

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_XX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_XY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_XZ)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_YX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_YY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_YZ)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_ZX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_ZY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_STRESS_ZZ)
    // *************** Continuum only END *************

    // Possible future blocks (no FEM) interaction

    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_BLOCK_CONTACT_FAILURE_ID)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTICLE_BLOCK_IF_INITIAL_CONTACT)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

    // DUMMIES INT AND DOUBLE VARIABLES

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_SWITCH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_3)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_4)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_5)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_6)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_7)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_8)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_9)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INT_DUMMY_10)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_FORCES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_LOCAL_AXES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_3)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_4)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_5)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_6)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_7)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_8)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DOUBLE_DUMMY_9)

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
