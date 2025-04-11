/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTT3SWAPFIX_HPP_
#define TESTT3SWAPFIX_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GammaG1CellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellLocationWriter.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "DiffusionForce.hpp"
#include "VolumeTrackingModifier.hpp"
#include "OffLatticeSimulation.hpp"

#include "VertexBoundaryRefinementModifier.hpp"

#include "CellSrnModel.hpp"


/**
 * These tests check and demonstrate simulation of vertex based models with edge  Srn models
 */
class TestT3SwapFix : public AbstractCellBasedTestSuite
{
  
public:

    /*
     * Test running vertex based model with edge based SRN.
     */
    void TestTemporaryT3SwapFix()
    {
     // make the simulation
        std::vector<CellPtr> cells;

        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh;
        VoronoiVertexMeshGenerator mesh_generator = VoronoiVertexMeshGenerator(8,8,0,1.0);
        p_mesh = mesh_generator.GetMesh();
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;

        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        //p_mesh->SetT2Threshold(0.0);
        p_mesh->SetCheckForT3Swaps(true);
        p_mesh->SetDistanceForT3SwapChecking(5.0);
        p_mesh->SetCheckForInternalIntersections(0.0);
        p_mesh->SetT2Threshold(0.00);

        for (auto node_iter = p_mesh->GetNodeIteratorBegin();
              node_iter != p_mesh->GetNodeIteratorEnd();
              ++node_iter)
              {
                node_iter->SetRadius(1.0);
              }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellLocationWriter>();

        double initial_target_area = 1.0;

        for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            // target areas
            cell_iter->GetCellData()->SetItem("target area", initial_target_area);

            // The following code can be removed and is simply due to system original test was running on
            boost::shared_ptr<AbstractCellProperty> p_vec_data(CellPropertyRegistry::Instance()->Get<CellVecData>());
            cell_iter->AddCellProperty(p_vec_data);
            TS_ASSERT(cell_iter->HasCellVecData());

            CellPropertyCollection parent_cell_property_collection = cell_iter->rGetCellPropertyCollection().GetPropertiesType<CellVecData>();
            boost::shared_ptr<CellVecData> p_parent_cell_vec_data = boost::static_pointer_cast<CellVecData>(parent_cell_property_collection.GetProperty());

            Vec item_1 = PetscTools::CreateAndSetVec(2, -17.3); // <-17.3, -17.3>
            cell_iter->GetCellVecData()->SetItem("item 1", item_1);
            // Remove above code
        } 

      OffLatticeSimulation<2> simulator(cell_population);

       simulator.SetNoBirth(true);
       simulator.SetOutputCellVelocities(true);
       simulator.SetOutputDivisionLocations(true);
       simulator.SetSamplingTimestepMultiple(10);
       simulator.SetDt(0.001);
       simulator.SetEndTime(50.00);

       cell_population.SetWriteCellVtkResults(true);
       
       // Make the diffusion force
        MAKE_PTR(DiffusionForce<2>, std_diffusion_force);
       // Make the Farhadifar force
        MAKE_PTR(FarhadifarForce<2>, p_force);

        // before passing the force to the simulation
        p_force->SetAreaElasticityParameter(1.0);
        p_force->SetPerimeterContractilityParameter(1.0 / 0.001);
        p_force->SetLineTensionParameter(-2.0 * 3.0 / 0.001);
        p_force->SetBoundaryLineTensionParameter(2.0);
        p_force->SetTargetAreaParameter(1.00);
        simulator.AddForce(p_force);

        // We need a FarhadifarType target area modifier
        MAKE_PTR(TargetAreaLinearGrowthModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR(VertexBoundaryRefinementModifier<2>, refinement_modifier);
        simulator.AddSimulationModifier(refinement_modifier);

        // Add random motion
        std_diffusion_force->SetAbsoluteTemperature(296.0);
        simulator.AddForce(std_diffusion_force);

        simulator.SetOutputDirectory("T3SwapFix");

        simulator.Solve();
    }

};


#endif /*TESTT3SWAPFIX_HPP_*/