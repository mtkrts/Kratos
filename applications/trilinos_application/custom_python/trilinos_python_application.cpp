//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-12-09 20:20:55 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>

#include "custom_python/add_trilinos_space_to_python.h" 
#include "custom_python/add_trilinos_convergence_criterias_to_python.h" 
#include "custom_python/add_trilinos_schemes_to_python.h" 
#include "custom_python/add_trilinos_linear_solvers_to_python.h" 
#include "custom_python/add_trilinos_strategies_to_python.h" 

//Trilinos includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"


// Project includes 
#include "includes/define.h"
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
// #include "add_trilinos_linear_solvers_to_python.h"
#include "includes/model_part.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
// #include "solving_strategies/schemes/scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_lagrangian_monolithic_scheme.h"
// #include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
// #include "custom_strategies/schemes/trilinos_predictorcorrector_velocity_bossak_scheme.h"

//convergence criterias
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"
// #include "solving_strategies/convergencecriterias/displacement_criteria.h"

//Builder And Solver
// #include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
/*#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"*/
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_vec.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_mixed.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

#include "external_includes/aztec_solver.h"
#include "external_includes/amesos_solver.h"
#include "external_includes/ml_solver.h"

//configuration files
#include "../../incompressible_fluid_application/custom_strategies/strategies/solver_configuration.h"
#include "custom_strategies/strategies/trilinos_fractionalstep_configuration.h"
#include "../../incompressible_fluid_application/custom_strategies/strategies/fractional_step_strategy.h"
#include "../../incompressible_fluid_application/incompressible_fluid_application.h"



namespace Kratos {

    namespace Python {

        //************************************************************************************************
        typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
        typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;



        BOOST_PYTHON_MODULE(KratosTrilinosApplication) {

 	class_<KratosTrilinosApplication,
                    KratosTrilinosApplication::Pointer,
                    bases<KratosApplication>, boost::noncopyable > ("KratosTrilinosApplication")
                    ;

		AddBasicOperations();
		AddConvergenceCriterias();
		AddSchemes();
		AddLinearSolvers();
		AddStrategies();


//             class_<TrilinosLinearSolverType, TrilinosLinearSolverType::Pointer > ("TrilinosLinearSolver")
                    /*		.def("Initialize",&TrilinosLinearSolverType::Initialize)
                                    .def("Solve",pointer_to_solve)
                                            //.def("",&LinearSolverType::)
                                            .def(self_ns::str(self))*/
                    ;



            typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;
            typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
            typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;



//             //********************************************************************
//             //********************************************************************
//             class_< TrilinosBaseSchemeType, boost::noncopyable >
//                     ("TrilinosScheme", init< >())
//                     .def("Initialize", &TrilinosBaseSchemeType::Initialize)
//                     .def("SchemeIsInitialized", &TrilinosBaseSchemeType::SchemeIsInitialized)
//                     .def("ElementsAreInitialized", &TrilinosBaseSchemeType::ElementsAreInitialized)
//                     .def("InitializeElements", &TrilinosBaseSchemeType::InitializeElements)
//                     .def("InitializeSolutionStep", &TrilinosBaseSchemeType::InitializeSolutionStep)
//                     .def("FinalizeSolutionStep", &TrilinosBaseSchemeType::FinalizeSolutionStep)
//                     .def("InitializeNonLinIteration", &TrilinosBaseSchemeType::InitializeNonLinIteration)
//                     .def("FinalizeNonLinIteration", &TrilinosBaseSchemeType::FinalizeNonLinIteration)
//                     .def("Predict", &TrilinosBaseSchemeType::Predict)
//                     .def("Update", &TrilinosBaseSchemeType::Update)
//                     .def("CalculateOutputData", &TrilinosBaseSchemeType::CalculateOutputData)
//                     .def("Clean", &TrilinosBaseSchemeType::Clean)
//                     .def("MoveMesh", MoveMesh)
//                     ;
// 
//             class_< TrilinosResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
//                     bases< TrilinosBaseSchemeType >, boost::noncopyable >
//                     (
//                     "TrilinosResidualBasedIncrementalUpdateStaticScheme", init< >()
//                     );
// 
//             class_< TrilinosResidualBasedLagrangianMonolithicScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
//                     bases< TrilinosBaseSchemeType >, boost::noncopyable >
//                     (
//                     "TrilinosResidualBasedLagrangianMonolithicScheme", init<int >()
//                     );
// 
//             class_< TrilinosResidualBasedLagrangianMonolithicScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
//                     bases< TrilinosBaseSchemeType >, boost::noncopyable >
//                     (
//                     "TrilinosResidualBasedLagrangianMonolithicScheme", init<int >()
//                     );
// 
//             typedef ResidualBasedPredictorCorrectorVelocityBossakScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme;
// 
//             class_< TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme,
//                     bases< TrilinosBaseSchemeType >, boost::noncopyable >
//                     (
//                     "TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme", init<double, double >()
//                     );
// 
//             class_< TrilinosPredictorCorrectorVelocityBossakScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
//                     bases< TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme >, boost::noncopyable >
//                     (
//                     "TrilinosPredictorCorrectorVelocityBossakScheme", init<double, double >()
//                     );

//             //********************************************************************
//             //********************************************************************
//             //convergence criteria base class
//             typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosConvergenceCriteria;
//             class_< TrilinosConvergenceCriteria, boost::noncopyable > ("TrilinosConvergenceCriteria", init<>())
//                     .def("SetActualizeRHSFlag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::SetActualizeRHSFlag)
//                     .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
//                     .def("PreCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PreCriteria)
//                     .def("PostCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PostCriteria)
//                     .def("Initialize", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Initialize)
//                     .def("InitializeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::InitializeSolutionStep)
//                     .def("FinalizeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::FinalizeSolutionStep)
//                     ;
// 
//             class_< TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
//                     bases< TrilinosConvergenceCriteria >,
//                     boost::noncopyable >
//                     ("TrilinosDisplacementCriteria", init< double, double, Epetra_MpiComm& >());
// 
//             class_< TrilinosUPCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
//                     bases< TrilinosConvergenceCriteria >,
//                     boost::noncopyable >
//                     ("TrilinosUPCriteria", init< double, double, double, double, Epetra_MpiComm& >());

//             //********************************************************************
//             //********************************************************************
//             //Builder and Solver
// 
//             class_< TrilinosBuilderAndSolverType::DofsArrayType, boost::noncopyable > ("DofsArrayType", init<>());
// 
//             class_< TrilinosBuilderAndSolverType, boost::noncopyable >
//                     ("TrilinosResidualBasedBuilderAndSolver", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > ())
//                     .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::SetCalculateReactionsFlag)
//                     .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::GetCalculateReactionsFlag)
//                     .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::SetDofSetIsInitializedFlag)
//                     .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::GetDofSetIsInitializedFlag)
//                     .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::SetReshapeMatrixFlag)
//                     .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::GetReshapeMatrixFlag)
//                     .def("GetEquationSystemSize", &TrilinosBuilderAndSolverType::GetEquationSystemSize)
//                     .def("BuildLHS", &TrilinosBuilderAndSolverType::BuildLHS)
//                     .def("BuildRHS", &TrilinosBuilderAndSolverType::BuildRHS)
//                     .def("Build", &TrilinosBuilderAndSolverType::Build)
//                     .def("SystemSolve", &TrilinosBuilderAndSolverType::SystemSolve)
//                     .def("BuildAndSolve", &TrilinosBuilderAndSolverType::BuildAndSolve)
//                     .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverType::BuildRHSAndSolve)
//                     .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverType::ApplyDirichletConditions)
//                     .def("SetUpDofSet", &TrilinosBuilderAndSolverType::SetUpDofSet)
//                     .def("GetDofSet", &TrilinosBuilderAndSolverType::GetDofSet, return_internal_reference<>())
//                     .def("SetUpSystem", &TrilinosBuilderAndSolverType::SetUpSystem)
//                     .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverType::ResizeAndInitializeVectors)
//                     .def("InitializeSolutionStep", &TrilinosBuilderAndSolverType::InitializeSolutionStep)
//                     .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverType::FinalizeSolutionStep)
//                     .def("CalculateReactions", &TrilinosBuilderAndSolverType::CalculateReactions)
//                     .def("Clear", &TrilinosBuilderAndSolverType::Clear)
//                     .def("SetEchoLevel", &TrilinosBuilderAndSolverType::SetEchoLevel)
//                     .def("GetEchoLevel", &TrilinosBuilderAndSolverType::GetEchoLevel)
//                     ;
// 
// 
//             typedef TrilinosBuilderAndSolverML< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverMLtype;
// 
//             class_< TrilinosBuilderAndSolverMLtype, bases<TrilinosBuilderAndSolverType>, boost::noncopyable >
//                     ("TrilinosBuilderAndSolverML", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > ())
//                     .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLtype::SetCalculateReactionsFlag)
//                     .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLtype::GetCalculateReactionsFlag)
//                     .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLtype::SetDofSetIsInitializedFlag)
//                     .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLtype::GetDofSetIsInitializedFlag)
//                     .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLtype::SetReshapeMatrixFlag)
//                     .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLtype::GetReshapeMatrixFlag)
//                     .def("GetEquationSystemSize", &TrilinosBuilderAndSolverMLtype::GetEquationSystemSize)
//                     .def("BuildLHS", &TrilinosBuilderAndSolverMLtype::BuildLHS)
//                     .def("BuildRHS", &TrilinosBuilderAndSolverMLtype::BuildRHS)
//                     .def("Build", &TrilinosBuilderAndSolverMLtype::Build)
//                     .def("SystemSolve", &TrilinosBuilderAndSolverMLtype::SystemSolve)
//                     .def("BuildAndSolve", &TrilinosBuilderAndSolverMLtype::BuildAndSolve)
//                     .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverMLtype::BuildRHSAndSolve)
//                     .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverMLtype::ApplyDirichletConditions)
//                     .def("SetUpDofSet", &TrilinosBuilderAndSolverMLtype::SetUpDofSet)
//                     .def("GetDofSet", &TrilinosBuilderAndSolverMLtype::GetDofSet, return_internal_reference<>())
//                     .def("SetUpSystem", &TrilinosBuilderAndSolverMLtype::SetUpSystem)
//                     .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverMLtype::ResizeAndInitializeVectors)
//                     .def("InitializeSolutionStep", &TrilinosBuilderAndSolverMLtype::InitializeSolutionStep)
//                     .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverMLtype::FinalizeSolutionStep)
//                     .def("CalculateReactions", &TrilinosBuilderAndSolverMLtype::CalculateReactions)
//                     .def("Clear", &TrilinosBuilderAndSolverMLtype::Clear)
//                     .def("SetEchoLevel", &TrilinosBuilderAndSolverMLtype::SetEchoLevel)
//                     .def("GetEchoLevel", &TrilinosBuilderAndSolverMLtype::GetEchoLevel)
//                     ;
// 
// 
// 
//             typedef TrilinosBuilderAndSolverML2D< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverML2Dtype;
// 
//             class_< TrilinosBuilderAndSolverML2Dtype, boost::noncopyable >
//                     ("TrilinosBuilderAndSolverML2D", init<Epetra_MpiComm&, int, int, TrilinosLinearSolverType::Pointer > ())
//                     .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverML2Dtype::SetCalculateReactionsFlag)
//                     .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverML2Dtype::GetCalculateReactionsFlag)
//                     .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverML2Dtype::SetDofSetIsInitializedFlag)
//                     .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverML2Dtype::GetDofSetIsInitializedFlag)
//                     .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverML2Dtype::SetReshapeMatrixFlag)
//                     .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverML2Dtype::GetReshapeMatrixFlag)
//                     .def("GetEquationSystemSize", &TrilinosBuilderAndSolverML2Dtype::GetEquationSystemSize)
//                     .def("BuildLHS", &TrilinosBuilderAndSolverML2Dtype::BuildLHS)
//                     .def("BuildRHS", &TrilinosBuilderAndSolverML2Dtype::BuildRHS)
//                     .def("Build", &TrilinosBuilderAndSolverML2Dtype::Build)
//                     .def("SystemSolve", &TrilinosBuilderAndSolverML2Dtype::SystemSolve)
//                     .def("BuildAndSolve", &TrilinosBuilderAndSolverML2Dtype::BuildAndSolve)
//                     .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverML2Dtype::BuildRHSAndSolve)
//                     .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverML2Dtype::ApplyDirichletConditions)
//                     .def("SetUpDofSet", &TrilinosBuilderAndSolverML2Dtype::SetUpDofSet)
//                     .def("GetDofSet", &TrilinosBuilderAndSolverML2Dtype::GetDofSet, return_internal_reference<>())
//                     .def("SetUpSystem", &TrilinosBuilderAndSolverML2Dtype::SetUpSystem)
//                     .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverML2Dtype::ResizeAndInitializeVectors)
//                     .def("InitializeSolutionStep", &TrilinosBuilderAndSolverML2Dtype::InitializeSolutionStep)
//                     .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverML2Dtype::FinalizeSolutionStep)
//                     .def("CalculateReactions", &TrilinosBuilderAndSolverML2Dtype::CalculateReactions)
//                     .def("Clear", &TrilinosBuilderAndSolverML2Dtype::Clear)
//                     .def("SetEchoLevel", &TrilinosBuilderAndSolverML2Dtype::SetEchoLevel)
//                     .def("GetEchoLevel", &TrilinosBuilderAndSolverML2Dtype::GetEchoLevel)
//                     ;
// 
// 
//             typedef TrilinosBuilderAndSolverMLmixed< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverMLmixedType;
// 
//             class_< TrilinosBuilderAndSolverMLmixedType, boost::noncopyable >
//                     ("TrilinosBuilderAndSolverMLmixed", init<Epetra_MpiComm&, int, int, TrilinosLinearSolverType::Pointer > ())
//                     .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLmixedType::SetCalculateReactionsFlag)
//                     .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLmixedType::GetCalculateReactionsFlag)
//                     .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLmixedType::SetDofSetIsInitializedFlag)
//                     .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLmixedType::GetDofSetIsInitializedFlag)
//                     .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLmixedType::SetReshapeMatrixFlag)
//                     .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLmixedType::GetReshapeMatrixFlag)
//                     .def("GetEquationSystemSize", &TrilinosBuilderAndSolverMLmixedType::GetEquationSystemSize)
//                     .def("BuildLHS", &TrilinosBuilderAndSolverMLmixedType::BuildLHS)
//                     .def("BuildRHS", &TrilinosBuilderAndSolverMLmixedType::BuildRHS)
//                     .def("Build", &TrilinosBuilderAndSolverMLmixedType::Build)
//                     .def("SystemSolve", &TrilinosBuilderAndSolverMLmixedType::SystemSolve)
//                     .def("BuildAndSolve", &TrilinosBuilderAndSolverMLmixedType::BuildAndSolve)
//                     .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverMLmixedType::BuildRHSAndSolve)
//                     .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverMLmixedType::ApplyDirichletConditions)
//                     .def("SetUpDofSet", &TrilinosBuilderAndSolverMLmixedType::SetUpDofSet)
//                     .def("GetDofSet", &TrilinosBuilderAndSolverMLmixedType::GetDofSet, return_internal_reference<>())
//                     .def("SetUpSystem", &TrilinosBuilderAndSolverMLmixedType::SetUpSystem)
//                     .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverMLmixedType::ResizeAndInitializeVectors)
//                     .def("InitializeSolutionStep", &TrilinosBuilderAndSolverMLmixedType::InitializeSolutionStep)
//                     .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverMLmixedType::FinalizeSolutionStep)
//                     .def("CalculateReactions", &TrilinosBuilderAndSolverMLmixedType::CalculateReactions)
//                     .def("Clear", &TrilinosBuilderAndSolverMLmixedType::Clear)
//                     .def("SetEchoLevel", &TrilinosBuilderAndSolverMLmixedType::SetEchoLevel)
//                     .def("GetEchoLevel", &TrilinosBuilderAndSolverMLmixedType::GetEchoLevel)
//                     ;

            /*		typedef TrilinosBuilderAndSolverML< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverMLtype;

                            class_< TrilinosBuilderAndSolverMLtype,bases<TrilinosBuilderAndSolverType>, boost::noncopyable >
                                    ("TrilinosBuilderAndSolverML",init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer >() )
                                    ;


                            typedef TrilinosBuilderAndSolverML2D< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverML2Dtype;

                            class_< TrilinosBuilderAndSolverML2Dtype,bases<TrilinosBuilderAndSolverMLtype>, boost::noncopyable >
                                    ("TrilinosBuilderAndSolverML2D",init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer >() )
                                    ;*/



//             typedef AztecSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AztecSolverType;
//             class_<AztecSolverType, bases<TrilinosLinearSolverType>, boost::noncopyable >
//                     ("AztecSolver",
//                     init< Teuchos::ParameterList&, std::string, Teuchos::ParameterList&, double, int, int >());
// 
//             typedef AmesosSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmesosSolverType;
//             class_<AmesosSolverType, bases<TrilinosLinearSolverType>, boost::noncopyable >
//                     ("AmesosSolver",
//                     init<const std::string&, Teuchos::ParameterList& >());
// 
//             typedef MultiLevelSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > MLSolverType;
//             class_<MLSolverType, bases<TrilinosLinearSolverType>, boost::noncopyable >
//                     ("MultiLevelSolver",
//                     init<Teuchos::ParameterList&, double, int >());


//             //********************************************************************************************
//             class_< SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >,
//                     boost::noncopyable >
//                     ("SolverConfiguration", init< ModelPart&, unsigned int>())
//                     .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
//                     .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
//                     ;
// 
//             class_< TrilinosFractionalStepConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >,
//                     bases< SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > >,
//                     boost::noncopyable >
//                     ("TrilinosFractionalStepConfiguration", init < Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, TrilinosLinearSolverType::Pointer,
//                     unsigned int, unsigned int, bool >());
// 
// 
//             //********************************************************************************************
//             class_< FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, boost::noncopyable >
//                     ("TrilinosFractionalStepStrategy",
//                     init < ModelPart&,
//                     SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >&,
//                     bool,
//                     double, double,
//                     int, int,
//                     unsigned int, unsigned int,
//                     bool
//                     >())
//                     .def("Solve", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Solve)
//                     .def("SolveStep1", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep1)
//                     .def("SolveStep2", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep2)
//                     .def("SolveStep3", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep3)
//                     .def("SolveStep4", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep4)
//                     .def("ActOnLonelyNodes", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ActOnLonelyNodes)
//                     .def("Clear", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Clear)
//                     .def("FractionalVelocityIteration", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::FractionalVelocityIteration)
//                     .def("ConvergenceCheck", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ConvergenceCheck)
//                     .def("InitializeFractionalStep", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::InitializeFractionalStep)
//                     .def("PredictVelocity", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::PredictVelocity)
//                     .def("InitializeProjections", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::InitializeProjections)
//                     .def("AssignInitialStepValues", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::AssignInitialStepValues)
//                     .def("IterativeSolve", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::IterativeSolve)
//                     .def("SavePressureIteration", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SavePressureIteration)
//                     .def("ApplyFractionalVelocityFixity", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ApplyFractionalVelocityFixity)
//                     .def("SetEchoLevel", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SetEchoLevel )
//             ;

        }


    } // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
