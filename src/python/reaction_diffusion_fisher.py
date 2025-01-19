#!/usr/bin/env python

import sys

# Intialise OpenCMISS
from opencmiss.opencmiss import OpenCMISS_Python as oc

#-----------------------------------------------------------------------------------------------------------
# SET PROBLEM PARAMETERS
#-----------------------------------------------------------------------------------------------------------

LINEAR_LAGRANGE = 1
QUADRATIC_LAGRANGE = 2
CUBIC_LAGRANGE = 3

NOSPLIT_SOLUTION = 1
GUDUNOV_SPLIT_SOLUTION = 2
STRANG_SPLIT_SOLUTION = 3

LENGTH = 1.0

ANALYTIC_FUNCTION_TYPE = oc.ReactionDiffusionAnalyticFunctionTypes.OC_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_ONE_DIM_1
A = 1.0
c = 2.1
aParam = 1.0
bParam = -c
cParam = -c
DParam = -1.0

TIME_START = 0.0
TIME_STOP = 1.0
PDE_TIME_STEP = 0.1
ODE_TIME_STEP = 0.01

OUTPUT_FREQUENCY = 1

(CONTEXT_USER_NUMBER,
 COORDINATE_SYSTEM_USER_NUMBER,
 REGION_USER_NUMBER,
 BASIS_USER_NUMBER,
 GENERATED_MESH_USER_NUMBER,
 MESH_USER_NUMBER,
 DECOMPOSITION_USER_NUMBER,
 DECOMPOSER_USER_NUMBER,
 GEOMETRIC_FIELD_USER_NUMBER,
 EQUATIONS_SET_FIELD_USER_NUMBER,
 EQUATIONS_SET_USER_NUMBER,
 DEPENDENT_FIELD_USER_NUMBER,
 MATERIALS_FIELD_USER_NUMBER,
 ANALYTIC_FIELD_USER_NUMBER,
 CELLML_USER_NUMBER,
 CELLML_MODELS_FIELD_USER_NUMBER,
 CELLML_STATE_FIELD_USER_NUMBER,
 CELLML_PARAMETERS_FIELD_USER_NUMBER,
 PROBLEM_USER_NUMBER) = range(1,20)

NUMBER_OF_GAUSS_XI = 2

numberOfElements = 10
interpolationOrder = LINEAR_LAGRANGE
#solutionMethod = GUDUNOV_SPLIT_SOLUTION
solutionMethod = NOSPLIT_SOLUTION

# Override with command line arguments if need be
if len(sys.argv) > 1:
    if len(sys.argv) > 4:
        sys.exit('ERROR: too many arguments- currently only accepting up to 3 options: numberOfElements interpolationOrder solutionMethod')
        numberOfElements = int(sys.argv[1])
    if len(sys.argv) > 2:
        interpolationOrder = int(sys.argv[2])
    if len(sys.argv) > 3:
        solutionMethod = int(sys.argv[3])
        
if (numberOfElements <= 0):
    sys.exit('ERROR: number of elements must be greater than 0.')
    
if ((interpolationOrder != LINEAR_LAGRANGE) and \
    (interpolationOrder != QUADRATIC_LAGRANGE) and \
    (interpolationOrder != CUBIC_LAGRANGE)):
    sys.exit('ERROR: invalid interpolation order.')
    
if ((solutionMethod != NOSPLIT_SOLUTION) and \
    (solutionMethod != GUDUNOV_SPLIT_SOLUTION) and \
    (solutionMethod != STRANG_SPLIT_SOLUTION)):
    sys.exit('ERROR: invalid solution method.')

#-----------------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#-----------------------------------------------------------------------------------------------------------

oc.OutputSetOn("FishersEquation")

#-----------------------------------------------------------------------------------------------------------
# CONTEXT
#-----------------------------------------------------------------------------------------------------------

# Create a context for the example
context = oc.Context()
context.Create(CONTEXT_USER_NUMBER)

# Get the world region
worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

# Get the computational nodes information
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

#-----------------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
#-----------------------------------------------------------------------------------------------------------

coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context)
coordinateSystem.DimensionSet(1)
coordinateSystem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# REGION
#-----------------------------------------------------------------------------------------------------------
region = oc.Region()
region.CreateStart(REGION_USER_NUMBER,worldRegion)
region.LabelSet("FishersEquation")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BASIS
#-----------------------------------------------------------------------------------------------------------

basis = oc.Basis()
basis.CreateStart(BASIS_USER_NUMBER,context)
basis.TypeSet(oc.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(1)
if (interpolationOrder == LINEAR_LAGRANGE):
    basis.InterpolationXiSet([oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE])
elif (interpolationOrder == QUADRATIC_LAGRANGE):
    basis.InterpolationXiSet([oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE])
elif (interpolationOrder == CUBIC_LAGRANGE):
    basis.InterpolationXiSet([oc.BasisInterpolationSpecifications.CUBIC_LAGRANGE])
else:
    sys.exit('ERROR: invalid interpolation order.')    
basis.QuadratureNumberOfGaussXiSet([NUMBER_OF_GAUSS_XI])
basis.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# MESH
#-----------------------------------------------------------------------------------------------------------
generatedMesh = oc.GeneratedMesh()
generatedMesh.CreateStart(GENERATED_MESH_USER_NUMBER,region)
generatedMesh.TypeSet(oc.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
generatedMesh.ExtentSet([LENGTH])
generatedMesh.NumberOfElementsSet([numberOfElements])
mesh = oc.Mesh()
generatedMesh.CreateFinish(MESH_USER_NUMBER,mesh)

#-----------------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
#-----------------------------------------------------------------------------------------------------------

decomposition = oc.Decomposition()
decomposition.CreateStart(DECOMPOSITION_USER_NUMBER,mesh)
decomposition.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# DECOMPOSER
#-----------------------------------------------------------------------------------------------------------

decomposer = oc.Decomposer()
decomposer.CreateStart(DECOMPOSER_USER_NUMBER,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
#-----------------------------------------------------------------------------------------------------------

geometricField = oc.Field()
geometricField.CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region)
geometricField.DecompositionSet(decomposition)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

#-----------------------------------------------------------------------------------------------------------
# EQUATION SETS
#-----------------------------------------------------------------------------------------------------------

# Create standard Laplace equations set
equationsSetField = oc.Field()
equationsSet = oc.EquationsSet()
if (solutionMethod == NOSPLIT_SOLUTION):
    equationsSetSpecification = [oc.EquationsSetClasses.CLASSICAL_FIELD,
                                 oc.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                 oc.EquationsSetSubtypes.GEN_FISHERS_NOSPLIT_REACT_DIFF]
elif ((solutionMethod == GUDUNOV_SPLIT_SOLUTION) or
      (solutionMethod == STRANG_SPLIT_SOLUTION)):
    equationsSetSpecification = [oc.EquationsSetClasses.CLASSICAL_FIELD,
                                 oc.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                 oc.EquationsSetSubtypes.CELLML_SPLIT_GEN_REACT_DIFF]
else:
    sys.exit('ERROR: invalid solution method.')
equationsSet.CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,
        equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField)
equationsSet.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# DEPENDENT FIELD
#-----------------------------------------------------------------------------------------------------------

dependentField = oc.Field()
equationsSet.DependentCreateStart(DEPENDENT_FIELD_USER_NUMBER,dependentField)
dependentField.DOFOrderTypeSet(oc.FieldVariableTypes.U,oc.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(oc.FieldVariableTypes.DELUDELN,oc.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,0.5)

#-----------------------------------------------------------------------------------------------------------
# MATERIALS FIELD
#-----------------------------------------------------------------------------------------------------------

materialsField = oc.Field()
equationsSet.MaterialsCreateStart(MATERIALS_FIELD_USER_NUMBER,materialsField)
if (solutionMethod == NOSPLIT_SOLUTION):
    materialsField.ComponentInterpolationSet(oc.FieldVariableTypes.U,1,oc.FieldInterpolationTypes.CONSTANT)
    materialsField.ComponentInterpolationSet(oc.FieldVariableTypes.U,2,oc.FieldInterpolationTypes.CONSTANT)
    materialsField.ComponentInterpolationSet(oc.FieldVariableTypes.U,3,oc.FieldInterpolationTypes.CONSTANT)
    materialsField.ComponentInterpolationSet(oc.FieldVariableTypes.U,4,oc.FieldInterpolationTypes.CONSTANT)
elif ((solutionMethod == GUDUNOV_SPLIT_SOLUTION) or
      (solutionMethod == STRANG_SPLIT_SOLUTION)):
    materialsField.ComponentInterpolationSet(oc.FieldVariableTypes.U,1,oc.FieldInterpolationTypes.CONSTANT)
else:
    sys.exit('ERROR: invalid solution method.')
equationsSet.MaterialsCreateFinish()

# Initialise materials field
if (solutionMethod == NOSPLIT_SOLUTION):
    materialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,aParam)
    materialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,bParam)
    materialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3,cParam)
    materialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,4,DParam)
elif ((solutionMethod == GUDUNOV_SPLIT_SOLUTION) or
      (solutionMethod == STRANG_SPLIT_SOLUTION)):
    materialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,DParam)
else:
    sys.exit('ERROR: invalid solution method.')

#-----------------------------------------------------------------------------------------------------------
# ANALYTIC FIELD
#-----------------------------------------------------------------------------------------------------------

analyticField = oc.Field()
equationsSet.AnalyticCreateStart(ANALYTIC_FUNCTION_TYPE,ANALYTIC_FIELD_USER_NUMBER,analyticField)
analyticField.ComponentInterpolationSet(oc.FieldVariableTypes.U,1,oc.FieldInterpolationTypes.CONSTANT)
analyticField.ComponentInterpolationSet(oc.FieldVariableTypes.U,2,oc.FieldInterpolationTypes.CONSTANT)
equationsSet.AnalyticCreateFinish()

# Use defaults for analytic field
analyticField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,A)
analyticField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,c)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS
#-----------------------------------------------------------------------------------------------------------

equations = oc.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(oc.EquationsSparsityTypes.SPARSE)
#equations.OutputTypeSet(oc.EquationsOutputTypes.NONE)
#equations.OutputTypeSet(oc.EquationsOutputTypes.MATRIX)
equations.OutputTypeSet(oc.EquationsOutputTypes.ELEMENT_MATRIX)
equationsSet.EquationsCreateFinish()

if ((solutionMethod == GUDUNOV_SPLIT_SOLUTION) or (solutionMethod == STRANG_SPLIT_SOLUTION)):
    
    #-----------------------------------------------------------------------------------------------------------
    # CELLML
    #-----------------------------------------------------------------------------------------------------------
    
    # Create the CellML environment
    cellML = oc.CellML()
    cellML.CreateStart(CELLML_USER_NUMBER,region)
    # Import the cell model from a file
    cellModel = cellML.ModelImport('./Fishers.cellml')
    # Now we have imported the model we are able to specify which variables from the model we want to set from OpenCMISS
    cellML.VariableSetAsKnown(cellModel,"main/a")
    cellML.VariableSetAsKnown(cellModel,"main/b")
    cellML.VariableSetAsKnown(cellModel,"main/c")
    # and variables to get from the CellML
    cellML.CreateFinish()

    #-----------------------------------------------------------------------------------------------------------
    # CELLML FIELD MAPS
    #-----------------------------------------------------------------------------------------------------------
    
    # Start the creation of CellML <--> OpenCMISS field maps
    cellML.FieldMapsCreateStart()
    # Now we can set up the field variable component <--> CellML model variable mappings.
    # Map u
    cellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U,1,oc.FieldParameterSetTypes.VALUES,
                                  cellModel,"main/u", oc.FieldParameterSetTypes.VALUES)
    cellML.CreateCellMLToFieldMap(cellModel,"main/u", oc.FieldParameterSetTypes.VALUES,
                                  dependentField,oc.FieldVariableTypes.U,1,oc.FieldParameterSetTypes.VALUES)
    #Finish the creation of CellML <--> OpenCMISS field maps
    cellML.FieldMapsCreateFinish()

    #-----------------------------------------------------------------------------------------------------------
    # CELLML MODELS FIELD
    #-----------------------------------------------------------------------------------------------------------
    
    cellMLModelsField = oc.Field()
    cellML.ModelsFieldCreateStart(CELLML_MODELS_FIELD_USER_NUMBER,cellMLModelsField)
    cellML.ModelsFieldCreateFinish()
    
    #-----------------------------------------------------------------------------------------------------------
    # CELLML STATE FIELD
    #-----------------------------------------------------------------------------------------------------------
    
    cellMLStateField = oc.Field()
    cellML.StateFieldCreateStart(CELLML_STATE_FIELD_USER_NUMBER,cellMLStateField)
    cellML.StateFieldCreateFinish()

    #-----------------------------------------------------------------------------------------------------------
    # CELLML PARAMETERS FIELD
    #-----------------------------------------------------------------------------------------------------------
    
    cellMLParametersField = oc.Field()
    cellML.ParametersFieldCreateStart(CELLML_PARAMETERS_FIELD_USER_NUMBER,cellMLParametersField)
    cellML.ParametersFieldCreateFinish()
    
    # Set the initial Cm values
    aComponent = cellML.FieldComponentGet(cellModel,oc.CellMLFieldTypes.PARAMETERS, "main/a")
    bComponent = cellML.FieldComponentGet(cellModel,oc.CellMLFieldTypes.PARAMETERS, "main/b")
    cComponent = cellML.FieldComponentGet(cellModel,oc.CellMLFieldTypes.PARAMETERS, "main/c")
    cellMLParametersField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,aComponent,aParam)
    cellMLParametersField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,bComponent,bParam)
    cellMLParametersField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,cComponent,cParam)
        
#-----------------------------------------------------------------------------------------------------------
# PROBLEM
#-----------------------------------------------------------------------------------------------------------

problem = oc.Problem()
if (solutionMethod == NOSPLIT_SOLUTION):
    problemSpecification = [oc.ProblemClasses.CLASSICAL_FIELD,
                            oc.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                            oc.ProblemSubtypes.NONLINEAR_NOSPLIT_REACT_DIFF]
elif (solutionMethod == GUDUNOV_SPLIT_SOLUTION):
    problemSpecification = [oc.ProblemClasses.CLASSICAL_FIELD,
                            oc.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                            oc.ProblemSubtypes.CELLML_GUDUNOV_SPLIT_REACT_DIFF]
elif (solutionMethod == STRANG_SPLIT_SOLUTION):
    problemSpecification = [oc.ProblemClasses.CLASSICAL_FIELD,
                            oc.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                            oc.ProblemSubtypes.CELLML_STRANG_SPLIT_REACT_DIFF]
else:
    sys.exit('ERROR: invalid solution method.')
problem.CreateStart(PROBLEM_USER_NUMBER,context,problemSpecification)
problem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# CONTROL LOOPS
#-----------------------------------------------------------------------------------------------------------

# Create control loops
problem.ControlLoopCreateStart()
controlLoop = oc.ControlLoop()
problem.ControlLoopGet([oc.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.TimesSet(TIME_START,TIME_STOP,PDE_TIME_STEP)
controlLoop.OutputTypeSet(oc.ControlLoopOutputTypes.TIMING)
controlLoop.TimeOutputSet(OUTPUT_FREQUENCY)
problem.ControlLoopCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVERS
#-----------------------------------------------------------------------------------------------------------

# Create problem solvers
if (solutionMethod == NOSPLIT_SOLUTION):
    dynamicSolver = oc.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,dynamicSolver)
    #pdeSolver.OutputTypeSet(oc.SolverOutputTypes.SOLVER)
    dynamicSolver.OutputTypeSet(oc.SolverOutputTypes.MATRIX)
    #dynamicSolver.LinearTypeSet(oc.LinearSolverTypes.ITERATIVE)
    #dynamicSolver.LinearIterativeAbsoluteToleranceSet(1.0E-12)
    #dynamicSolver.LinearIterativeRelativeToleranceSet(1.0E-12)
    problem.SolversCreateFinish()
elif (solutionMethod == GUDUNOV_SPLIT_SOLUTION):
    odeSolver1 = oc.Solver()
    dynamicSolver = oc.Solver()
    problem.SolversCreateStart()
    # Get the first ODE solver
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,odeSolver1)
    odeSolver1.DAETimeStepSet(ODE_TIME_STEP)
    odeSolver1.OutputTypeSet(oc.SolverOutputTypes.NONE)
    # Get the second dynamic solver for the parabolic problem
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],2,dynamicSolver)
    #dynamicSolver.OutputTypeSet(oc.SolverOutputTypes.NONE)
    dynamicSolver.OutputTypeSet(oc.SolverOutputTypes.MATRIX)
    problem.SolversCreateFinish()
elif (solutionMethod == STRANG_SPLIT_SOLUTION):
    odeSolver1 = oc.Solver()
    odeSolver2 = oc.Solver()
    dynamicSolver = oc.Solver()
    problem.SolversCreateStart()
    # Get the first ODE solver
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,odeSolver1)
    odeSolver1.DAETimeStepSet(ODE_TIME_STEP)
    odeSolver1.OutputTypeSet(oc.SolverOutputTypes.NONE)
    # Get the second dynamic solver for the parabolic problem
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],2,dynamicSolver)
    #dynamicSolver.OutputTypeSet(oc.SolverOutputTypes.NONE)
    dynamicSolver.OutputTypeSet(oc.SolverOutputTypes.MATRIX)
    # Get the second ODE solver
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,odeSolver2)
    odeSolver2.DAETimeStepSet(ODE_TIME_STEP)
    odeSolver2.OutputTypeSet(oc.SolverOutputTypes.NONE)
    problem.SolversCreateFinish()
else:
    sys.exit('ERROR: invalid solution method.')

#-----------------------------------------------------------------------------------------------------------
# SOLVER EQUATIONS
#-----------------------------------------------------------------------------------------------------------

# Create solver equations and add equations set to solver equations
dynamicSolver = oc.Solver()
solverEquations = oc.SolverEquations()
problem.SolverEquationsCreateStart()
if (solutionMethod == NOSPLIT_SOLUTION):
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,dynamicSolver)
elif ((solutionMethod == GUDUNOV_SPLIT_SOLUTION) or
    (solutionMethod == STRANG_SPLIT_SOLUTION)):    
    problem.SolverGet([oc.ControlLoopIdentifiers.NODE],2,dynamicSolver)
else:
    sys.exit('ERROR: invalid solution method.')
dynamicSolver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(oc.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# CELLML EQUATIONS
#-----------------------------------------------------------------------------------------------------------

if ((solutionMethod == GUDUNOV_SPLIT_SOLUTION) or
    (solutionMethod == STRANG_SPLIT_SOLUTION)):    

    # Create CellML equations and add the CellML environment
    problem.CellMLEquationsCreateStart()
    cellMLEquations1 = oc.CellMLEquations()
    odeSolver1.CellMLEquationsGet(cellMLEquations1)
    cellmlIndex1 = cellMLEquations1.CellMLAdd(cellML)
    if (solutionMethod == STRANG_SPLIT_SOLUTION):
        cellMLEquations2 = oc.CellMLEquations()
        odeSolver2.CellMLEquationsGet(cellMLEquations2)
        cellmlIndex2 = cellMLEquations2.CellMLAdd(cellML)    
    problem.CellMLEquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------------------------

# Create analytic boundary conditions
boundaryConditions = oc.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
solverEquations.BoundaryConditionsAnalytic()
solverEquations.BoundaryConditionsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVE
#-----------------------------------------------------------------------------------------------------------

problem.Solve()

#-----------------------------------------------------------------------------------------------------------
# ANALYTIC ANALYSIS
#-----------------------------------------------------------------------------------------------------------

dependentField.AnalyticAnalysis_Output('FishersEquation')

#-----------------------------------------------------------------------------------------------------------
# OUTPUT
#-----------------------------------------------------------------------------------------------------------

# Export results
fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("FishersEquation","FORTRAN")
fields.ElementsExport("FishersEquation","FORTRAN")
fields.Finalise()

# Destroy the context
context.Destroy()
# Finalise OpenCMISS
oc.Finalise()
