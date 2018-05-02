//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes


// External includes


// Project includes
#include "custom_conditions/meshless_force_interface_condition.h"
#include "utilities/math_utils.h"
#include "includes/define.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"


namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	MeshlessForceInterfaceCondition::MeshlessForceInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseCondition(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}


	//************************************************************************************
	//************************************************************************************
	MeshlessForceInterfaceCondition::MeshlessForceInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: MeshlessBaseCondition(NewId, pGeometry, pProperties)
	{
	}


	Condition::Pointer MeshlessForceInterfaceCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return MeshlessBaseCondition::Pointer(new MeshlessForceInterfaceCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}


	// Destructor
	MeshlessForceInterfaceCondition::~MeshlessForceInterfaceCondition()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessForceInterfaceCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		MatrixType temp(0, 0);
		CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
	}

	//************************************************************************************
	//************************************************************************************
	/**
	* CalculateLocalSystem
	* Returns only rRightHandSideVector as there is no impact on the 
	* siffness due to loads
	*/
	void MeshlessForceInterfaceCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
			std::cout << "start MeshlessForceInterfaceCondition" << std::endl;
		const unsigned int number_of_points = GetGeometry().size();

		if (rLeftHandSideMatrix.size1() != number_of_points * 3)
			rLeftHandSideMatrix.resize(number_of_points * 3, number_of_points * 3, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points * 3, number_of_points * 3); //resetting LHS

		if (rRightHandSideVector.size() != number_of_points * 3)
			rRightHandSideVector.resize(number_of_points * 3, false);
		rRightHandSideVector = ZeroVector(number_of_points * 3); //resetting RHS

		const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
		const Vector& force_vector = this->GetValue(EXTERNAL_FORCES_VECTOR);
		KRATOS_WATCH(N)
		KRATOS_WATCH(force_vector)
		Vector fLoads = ZeroVector(number_of_points * 3);

		for (unsigned int i = 0; i < number_of_points; i++)
		{
			int index = 3 * i;
			fLoads[index]	  = - force_vector[0] * N[i];
			fLoads[index + 1] = - force_vector[1] * N[i];
			fLoads[index + 2] = - force_vector[2] * N[i];
		}
		noalias(rRightHandSideVector) -= fLoads;

		KRATOS_WATCH(rRightHandSideVector)
		KRATOS_CATCH("")
	}
	/***********************************************************************************/
	/***********************************************************************************/

	int  MeshlessForceInterfaceCondition::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		if (!Has(SHAPE_FUNCTION_VALUES))
			KRATOS_ERROR << "No SHAPE_FUNCTION_VALUES assigned!" << std::endl;
		if (!this->Has(EXTERNAL_FORCES_VECTOR))
			KRATOS_ERROR << "EXTERNAL_FORCES_VECTOR not assigned!" << std::endl;
		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MeshlessForceInterfaceCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = number_of_nodes * 3;

		if (rResult.size() != dim)
			rResult.resize(dim);

		for (unsigned int i = 0; i < number_of_nodes; i++)
		{
			int index = i * 3;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
		}

		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void MeshlessForceInterfaceCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
	{

		ElementalDofList.resize(0);

		for (unsigned int i = 0; i < GetGeometry().size(); i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessForceInterfaceCondition::GetFirstDerivativesVector(Vector& values, int Step)
	{
		if (values.size() != 3)
			values.resize(3, false);

		const int& number_of_nodes = GetGeometry().size();
		Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

		for (SizeType i = 0; i < number_of_nodes; i++)
		{
			const NodeType & iNode = GetGeometry()[i];
			const array_1d<double, 3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

			values[0] += N[i] * vel[0];
			values[1] += N[i] * vel[1];
			values[2] += N[i] * vel[2];
		}
	}

} // Namespace Kratos