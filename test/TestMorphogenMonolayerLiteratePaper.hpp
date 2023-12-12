#ifndef TESTMORPHOGENMONOLAYERLITERATEPAPER_HPP_
#define TESTMORPHOGENMONOLAYERLITERATEPAPER_HPP_


/*
 * = Long-range Signalling Example =
 *
 * On this wiki page we describe in detail the code that is used to run this example from the paper.
 *
 * The easiest way to visualize these simulations is with Paraview.
 * 
 * [[EmbedYoutube(Yl2GT2x2ohc)]]
 *
 * == Code overview ==
 *
 * The first thing to do is to include the necessary header files.
 */

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "DefaultCellProliferativeType.hpp"

#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "TissueWidthWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "MorphogenCellwiseSourceParabolicPde.hpp"
#include "VolumeTrackingModifier.hpp"

#include "FixedDurationCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForceWithMinDistanceItem.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GrowthInhibitionModifier.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "CaBasedCellPopulation.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "AdhesionCaSwitchingUpdateRule.hpp"

#include "RandomNumberGenerator.hpp"

#include "PetscSetupAndFinalize.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 *
 *  The first block (commented out) are the original parameter values.
 *  The second block are parameters for a much shorter simulation, and are used for continuous testing with Chaste.
 */

//static const double M_TIME_FOR_SIMULATION = 100; //100
//static const double M_NUM_CELLS_ACROSS = 10; // 10
//static const double M_UPTAKE_RATE = 0.01; // S in paper
//static const double M_DIFFUSION_CONSTANT = 1e-4; // D in paper
//static const double M_DUDT_COEFFICIENT = 1.0; // Not used in paper so 1

static const double M_TIME_FOR_SIMULATION = 1440.0;
static const double M_NUM_CELLS_ACROSS = 57;
static const double M_UPTAKE_RATE = 0.01; // S in paper
static const double M_DIFFUSION_CONSTANT = 1e-4; // D in paper
static const double M_DUDT_COEFFICIENT = 1.0;

class TestMorphogenMonolayerLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
private:

    /*
     * This is a helper method to generate cells and is used in all simulations.
     */ 
    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        for (unsigned i=0; i<num_cells; i++)
        {
            //UniformlyDistributedCellCycleModel* p_cycle_model = new UniformlyDistributedCellCycleModel();
            FixedDurationCellCycleModel* p_cycle_model = new FixedDurationCellCycleModel();
            p_cycle_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            // Note the first few recorded ages will be too short as cells start with some mass.
            //const double birth_time = -RandomNumberGenerator::Instance()->ranf() * 18.0;
            const double birth_time = -20;
            p_cell->SetBirthTime(birth_time);
            p_cycle_model->SetPhaseTimer(birth_time);


            p_cell->InitialiseCellCycleModel();

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            p_cell->GetCellData()->SetItem("growth inhibited", 0.0);
            p_cell->GetCellData()->SetItem("Radius", 0.1);
            p_cell->GetCellData()->SetItem("cell age", birth_time);
            rCells.push_back(p_cell);
        }
     }

public:

    /*
     * == OS ==
     *
     * Simulate reaction diffusion on a growing a population of cells in the
     * Overlapping Spheres model.
     */
    void TestNodeBasedMorphogenMonolayer()
    {
        HoneycombMeshGenerator generator(2.0 * M_NUM_CELLS_ACROSS, 2.0 * M_NUM_CELLS_ACROSS,0);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();

        p_generating_mesh->Translate(-M_NUM_CELLS_ACROSS / 2.0, -M_NUM_CELLS_ACROSS / 2.0);

        //Remove all elements outside the specified initial radius
        for (AbstractMesh<2, 2>::NodeIterator node_iter = p_generating_mesh->GetNodeIteratorBegin();
             node_iter != p_generating_mesh->GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            c_vector<double,2> node_location = node_iter->rGetLocation();

            if (norm_2(node_location)>0.5*M_NUM_CELLS_ACROSS + 1e-5)
            {
                p_generating_mesh->DeleteNodePriorToReMesh(node_index);
            }
        }
        p_generating_mesh->ReMesh();

        double cut_off_length = 1.5; //this is the default

        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, cut_off_length);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddPopulationWriter<TissueWidthWriter>();
        cell_population.SetUseVariableRadii(true);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MorphogenMonolayer/Node");
        simulator.SetDt(0.05);
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForceWithMinDistanceItem<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(5.00);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);
        
        MAKE_PTR(GrowthInhibitionModifier<2>, p_growth_inhibition_modifier);
        simulator.AddSimulationModifier(p_growth_inhibition_modifier);


        simulator.Solve();

        delete p_mesh; // to stop memory leaks

    }

    /*
     * == VT ==
     *
     * Simulate reaction diffusion on a growing a population of cells in the
     * Voronoi Tesselation model.
     */

    void TestMeshBasedMorphogenMonolayer()
    {
    }

    /*
     * == VM ==
     *
     * Simulate reaction diffusion on a growing a population of cells in the
     * Cell Vertex model.
     */
    void TestVertexBasedMorphogenMonolayer()
    {
    //    // Create Mesh
    //    HoneycombVertexMeshGenerator generator(2.0*M_NUM_CELLS_ACROSS, 3.0*M_NUM_CELLS_ACROSS);
    //    boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
    //    p_mesh->SetCellRearrangementThreshold(0.1);

    //    p_mesh->Translate(-M_NUM_CELLS_ACROSS,-sqrt(3.0)*M_NUM_CELLS_ACROSS+ sqrt(3.0)/6.0);

    //    //Remove all elements outside the specified initial radius
    //    for (VertexMesh<2,2>::VertexElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
    //             elem_iter != p_mesh->GetElementIteratorEnd();
    //             ++elem_iter)
    //    {
    //        unsigned elem_index = elem_iter->GetIndex();
    //        c_vector<double,2> element_centre = p_mesh->GetCentroidOfElement(elem_index);

    //        if (norm_2(element_centre)>0.5*M_NUM_CELLS_ACROSS + 1e-5)
    //        {
    //            p_mesh->DeleteElementPriorToReMesh(elem_index);
    //        }
    //    }
    //    p_mesh->ReMesh();

    //    // Create Cells
    //    std::vector<CellPtr> cells;
    //    GenerateCells(p_mesh->GetNumElements(),cells);

    //    // Create Population
    //    VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
    //    cell_population.AddCellWriter<CellIdWriter>();
    //    cell_population.AddCellWriter<CellAgesWriter>();
    //    cell_population.AddCellWriter<CellMutationStatesWriter>();
    //    //Make cell data writer so can pass in variable name
    //    boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
    //    cell_population.AddCellWriter(p_cell_data_item_writer);


    //    // Create Simulation
    //    OffLatticeSimulation<2> simulator(cell_population);
    //    simulator.SetOutputDirectory("MorphogenMonolayer/Vertex");
    //    simulator.SetDt(1.0/200.0);
    //    simulator.SetSamplingTimestepMultiple(200);
    //    simulator.SetEndTime(M_TIME_FOR_SIMULATION);

    //    simulator.SetOutputDivisionLocations(true);

    //    // Create Forces and pass to simulation NOTE: these are not the default ones and chosen to give a stable growing monolayer
    //    MAKE_PTR(NagaiHondaForce<2>, p_force);
    //    p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
    //    p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
    //    p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
    //    p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
    //    simulator.AddForce(p_force);

    //    // Create Modifiers and pass to simulation

    //    // Create a pde modifier and pass it to the simulation 

    //    // Make the Pde and BCS
    //    MAKE_PTR_ARGS(MorphogenCellwiseSourceParabolicPde<2>, p_pde, (cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE));
    //    MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

    //    // Create a PDE Modifier object using this pde and bcs object
    //    MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true));
    //    p_pde_modifier->SetDependentVariableName("morphogen");
    //    simulator.AddSimulationModifier(p_pde_modifier);

    //    simulator.Solve();
    }
};

#endif /* TESTMORPHOGENMONOLAYERLITERATEPAPER_HPP_ */
