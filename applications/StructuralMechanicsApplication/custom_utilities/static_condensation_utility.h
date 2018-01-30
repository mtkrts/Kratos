//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter
//                   Philipp Bucher
//

#if !defined(KRATOS_STATIC_CONDENSATION_UTILITY_H_INCLUDED )
#define  KRATOS_STATIC_CONDENSATION_UTILITY_H_INCLUDED


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"


namespace Kratos
{
	namespace StaticCondensationUtility
	{
        typedef Element ElementType;
        typedef std::size_t SizeType;
        typedef Matrix MatrixType;

		void CondenseLeftHandSide(
            const ElementType& rTheElement,
			MatrixType& rLeftHandSideMatrix,
			const std::vector<int> & rDofList);

		std::vector<MatrixType> CalculateSchurComplements(
            const ElementType& rTheElement,
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int> & rDofList);

		std::vector<int> CreateRemainingDofList(
            const ElementType& rTheElement,
			const std::vector<int> & rDofList);

		void FillSchurComplements(
			MatrixType& Submatrix,
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int>& rVecA,
			const std::vector<int>& rVecB,
			const SizeType& rSizeA,
			const SizeType& rSizeB); //maybe inline

		void ConvertingCondensation(
            const ElementType& rTheElement,
			Vector& rValues,
			const std::vector<int>& rDofList,
			const MatrixType& rLeftHandSideMatrix);

	}  // namespace StaticCondensationUtility
  
}  // namespace Kratos.

#endif // KRATOS_STATIC_CONDENSATION_UTILITY_H_INCLUDED  defined 


