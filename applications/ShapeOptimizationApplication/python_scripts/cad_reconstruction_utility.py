# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import json as json
import time
import os

# ======================================================================================================================================
# Class definition
# ======================================================================================================================================

class CADReconstrutionUtilities():

    # --------------------------------------------------------------------------
    def __init__( self, reconstruction_parameters ):
        self.Parameters = reconstruction_parameters

    # --------------------------------------------------------------------------
    def Initialize( self ):
        self.__ReadFEData()
        self.__RefineFEModel()
        self.__ReadCADData()
        self.__InitializeConstraints()
        self.__CreateReconstructionDataBase()
        self.__CreateReconstructionOutputWriter()
    
    # --------------------------------------------------------------------------
    def SetDisplacementCouplingOnAllCouplingPoints( self ):
        self.IsDisplacementCouplingSpecifiedForAllCouplingPoints = True

    # --------------------------------------------------------------------------
    def SetRotationCouplingOnAllCouplingPoints( self ):
        self.IsRotationCouplingSpecifiedForAllCouplingPoints = True
        
    # --------------------------------------------------------------------------
    def SetDirichletConstraints( self, edges_with_dirichlet_conditions ):
        self.AreDirichletConstraintsSpecified = True
        self.EdgesWithDirichletConditions = edges_with_dirichlet_conditions

    # --------------------------------------------------------------------------
    def PerformReconstruction( self ):
        self.__CreateReconstructionConditions()
        self.__CreateSolverForReconstruction()
        self.__RunSolutionAlorithm()
        
    # --------------------------------------------------------------------------
    def OutputFEData( self ):
        fem_input_filename = self.Parameters["inpute_parameters"]["fem_filename"].GetString()
        output_folder = self.Parameters["result_output_parameters"]["results_output_folder"].GetString()

        from gid_output import GiDOutput
        fem_output_filename = output_folder + "/" + fem_input_filename + "_as_used_for_reconstruction"
        nodal_results=["SHAPE_CHANGE_ABSOLUTE"]
        gauss_points_results=[]
        VolumeOutput = True
        GiDPostMode = "Binary"
        GiDWriteMeshFlag = False
        GiDWriteConditionsFlag = True
        GiDWriteParticlesFlag = False
        GiDMultiFileFlag = "Single"

        gig_io = GiDOutput(fem_output_filename, VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
        gig_io.initialize_results(self.FEModelPart)
        gig_io.write_results(1, self.FEModelPart, nodal_results, gauss_points_results)
        gig_io.finalize_results()

    # --------------------------------------------------------------------------
    def OutputCADSurfacePoints( self, file_to_write ):
        self.OutputWriter.OutputCADSurfacePoints( file_to_write, self.Parameters )

    # --------------------------------------------------------------------------
    def OutputGaussPointsOfFEMesh( self, file_to_write ):
        self.OutputWriter.OutputGaussPointsOfFEMesh( file_to_write, self.FEMGaussIntegrationDegree )

    # --------------------------------------------------------------------------
    def __ReadFEData( self ):
        print("\n> Start importing FE data")
        fem_input_filename = self.Parameters["inpute_parameters"]["fem_filename"].GetString()
        self.FEModelPart = ModelPart("name_of_empty_mdpa")
        self.FEModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        model_part_io = ModelPartIO(fem_input_filename)
        model_part_io.ReadModelPart(self.FEModelPart)
        print("> Importing FE data finished.")        

    # --------------------------------------------------------------------------
    def __RefineFEModel( self ):       
        # Assign pseudo material to elements (required by refinement)
        prop_id = 1
        prop = self.FEModelPart.Properties[prop_id]
        mat = LinearElasticPlaneStress2DLaw()
        prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

        refinement_level = self.Parameters["inpute_parameters"]["fe_refinement_level"].GetInt()
        for refinement_level in range(0,refinement_level):

            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            nodal_neighbour_search = FindNodalNeighboursProcess(self.FEModelPart, number_of_avg_elems, number_of_avg_nodes)
            neighbour_calculator = FindElementalNeighboursProcess(self.FEModelPart,2,10)
            nodal_neighbour_search.Execute()
            neighbour_calculator.Execute()

            for elem in self.FEModelPart.Elements:
                elem.SetValue(SPLIT_ELEMENT,True)

            refine_on_reference = False
            interpolate_internal_variables = True
            Refine = LocalRefineTriangleMesh(self.FEModelPart)
            Refine.LocalRefineMesh(refine_on_reference, interpolate_internal_variables)

    # --------------------------------------------------------------------------
    def __ReadCADData( self ):
        print("\n> Start importing CAD data")
        cad_geometry_filename = self.Parameters["inpute_parameters"]["cad_geometry_filename"].GetString()        
        self.CADGeometry = {}
        with open(cad_geometry_filename) as cad_data1:
            self.CADGeometry = json.load(cad_data1)

        cad_integration_data_filename = self.Parameters["inpute_parameters"]["cad_integration_data_filename"].GetString()                    
        self.CADIntegrationData = {}
        with open(cad_integration_data_filename) as cad_data2:
            self.CADIntegrationData = json.load(cad_data2)
        print("> Importing CAD data finished.")

    # --------------------------------------------------------------------------
    def __InitializeConstraints( self ):
        self.IsDisplacementCouplingSpecifiedForAllCouplingPoints = False
        self.IsRotationCouplingSpecifiedForAllCouplingPoints = False
        self.AreDirichletConstraintsSpecified = False
    
    # --------------------------------------------------------------------------
    def __CreateReconstructionDataBase( self ):
        self.DataBase = ReconstructionDataBase(self.FEModelPart, self.CADGeometry, self.CADIntegrationData)
        self.DataBase.Create()

    # --------------------------------------------------------------------------
    def __CreateReconstructionOutputWriter( self ):
        output_folder = self.Parameters["result_output_parameters"]["results_output_folder"].GetString()
        if not os.path.exists( output_folder ):
            os.makedirs( output_folder )    
        self.OutputWriter = ReconstructionOutputWriter( self.DataBase, self.Parameters )    

    # --------------------------------------------------------------------------
    def __CreateReconstructionConditions( self ):
        # Container to store all conditions (including constraints and reguarlization )
        self.ConditionsContainer = ReconstructionConditionContainer( self.DataBase, self.Parameters )

        # Basic reconstruction condition
        reconstruction_strategy = self.Parameters["solution_parameters"]["strategy"].GetString()
        if reconstruction_strategy == "displacement_mapping":
            self.ConditionsContainer.CreateDisplacementMappingConditions()
        elif reconstruction_strategy == "distance_minimization":
            self.ConditionsContainer.CreateDistanceMinimizationConditions()                                                                          
        else:
            raise ValueError( "The following reconstruction strategy does not exist: ", reconstruction_strategy )

        # Reconstruction constraints
        if self.IsDisplacementCouplingSpecifiedForAllCouplingPoints: 
            self.ConditionsContainer.CreateDisplacementCouplingConstraintsOnAllCouplingPoints()
        if self.IsRotationCouplingSpecifiedForAllCouplingPoints: 
            self.ConditionsContainer.CreateRotationCouplingConstraintsOnAllCouplingPoints()            
        if self.AreDirichletConstraintsSpecified:
            self.ConditionsContainer.CreateDirichletConstraints( self.DirichletConstraints )

        # Regularization
        minimize_control_point_distance_to_surface = self.Parameters["solution_parameters"]["regularization_parameters"]["minimize_control_point_distance_to_surface"].GetBool()
        minimize_control_point_displacement = self.Parameters["solution_parameters"]["regularization_parameters"]["minimize_control_point_displacement"].GetBool()
        if minimize_control_point_distance_to_surface: 
            self.ConditionsContainer.CreateMinimalControlPointDistanceToSurfaceCondition()                               
        if minimize_control_point_displacement: 
            self.ConditionsContainer.CreateMinimalControlPointDisplacementCondition()

    # --------------------------------------------------------------------------
    def __CreateSolverForReconstruction( self ):
        linear_solver_name = self.Parameters["solution_parameters"]["linear_solver_name"].GetString()
        self.LinearSolver = None
        if linear_solver_name == "SuperLU":
            self.LinearSolver = SuperLUSolver()
        elif linear_solver_name == "BICGSTAB":
            DiagPrecond = DiagonalPreconditioner()
            self.LinearSolver =  BICGSTABSolver(1e-9, 5000, DiagPrecond)
        elif linear_solver_name == "AMGCL":
            self.LinearSolver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)        
        else:
            raise NameError("Linear solver not implemented!")              
        self.ReconstructionSolver = CADReconstructionSolver( self.DataBase, self.ConditionsContainer, self.LinearSolver, self.Parameters )

    # --------------------------------------------------------------------------
    def __RunSolutionAlorithm( self ):
        solution_iterations = self.Parameters["solution_parameters"]["general_parameters"]["solution_iterations"].GetInt()        
        penalty_multiplier =  self.Parameters["solution_parameters"]["general_parameters"]["penalty_multiplier"].GetDouble()
        self.ReconstructionSolver.InitializeEquationSystem()

        for iteration in range(1,solution_iterations+1): 
            print("\n===========================================")
            print("Starting reconstruction iteration ", iteration,"...")
            print("===========================================")            

            if iteration > 1:  
                self.ReconstructionSolver.MultiplyAllPenaltyFactorsByInputFactor( self.PenaltyMultiplier )               

            self.ReconstructionSolver.ComputeLHS()
            self.ReconstructionSolver.ComputeRHS()
            self.ReconstructionSolver.SolveEquationSystem()

            self.ReconstructionSolver.UpdateControlPointsAccordingReconstructionStrategy()
            
            self.OutputWriter.OutputResultsInRhinoFormat( iteration ) 

        print("\n===========================================")
        print("Finished reconstruction loop.")
        print("===========================================")                    