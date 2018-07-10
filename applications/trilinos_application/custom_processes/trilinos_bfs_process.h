//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#pragma once

// System includes
#include <list>

// External includes
#include "mpi.h"

// Project includes
#include "includes/node.h"
#include "includes/model_part.h"

// Processes
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos {

/**
 * @brief Probably we will have to change the name
 * 
 */
class TrilinosBfsProcess : public Process {
public:

    /// Pointer definition of TrilinosBfsProcess
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosBfsProcess);

    TrilinosBfsProcess(ModelPart & modelPart) 
    : mrModelPart(modelPart) {}

    void SetInitialCoord();

    void SetInitialNode();

    void SetInitialElement();

    void SetInitialCondition();

    void Execute() {
        
        Node<3>::WeakPointer originNodePtr = *(mrModelPart.NodesArray().begin());

        std::list<Node<3>::WeakPointer> modelNodes;
        std::list<Node<3>::WeakPointer> componentNodes;

        // Locate the elements for every node 
        auto findNeighbourProcesses = FindNodalNeighboursProcess(mrModelPart);
        findNeighbourProcesses.Execute();

        // Push the first node to initilize the modelNodes
        modelNodes.push_back(originNodePtr);

        std::size_t clusterID = 1;
        std::size_t visitId = 1;
        bool updateClusterID = false;

        while(!modelNodes.empty()) {
            auto & itModelNodePtr = modelNodes.front();
            modelNodes.pop_front();
            int visited = itModelNodePtr.lock()->FastGetSolutionStepValue(LOAD_RESTART);
            if(visited < 1) {
                std::cout << "New non-empty component" << std::endl;
                componentNodes.push_back(itModelNodePtr);
                while(!componentNodes.empty()) {
                    auto & itCompNodePtr = componentNodes.front();
                    componentNodes.pop_front();
                    itCompNodePtr.lock()->FastGetSolutionStepValue(PARTITION_INDEX) = visitId++;


                    int isXTOKEN1 = itCompNodePtr.lock()->FastGetSolutionStepValue(STATIONARY);
                    if(isXTOKEN1 < 1) {
                        itCompNodePtr.lock()->FastGetSolutionStepValue(RIGID_BODY_ID) = clusterID;
                        updateClusterID = true;
                    } else {
                        itCompNodePtr.lock()->FastGetSolutionStepValue(RIGID_BODY_ID) = 0;
                    }
                    auto & neighbourNodes = itCompNodePtr.lock()->GetValue(NEIGHBOUR_NODES);

                    for(auto neigItr = neighbourNodes.ptr_begin(); neigItr != neighbourNodes.ptr_end(); neigItr++) {
                        if((*neigItr).lock()->FastGetSolutionStepValue(LOAD_RESTART) < 1) {
                            int isXTOKEN2 = (*neigItr).lock()->FastGetSolutionStepValue(STATIONARY);

                            // It doesn't matter at which connected comp we are, just that is the same 
                            if(isXTOKEN1 != isXTOKEN2) {
                                modelNodes.push_back((*neigItr));
                            } else {
                                (*neigItr).lock()->FastGetSolutionStepValue(LOAD_RESTART) = 2;
                                componentNodes.push_back((*neigItr));
                            }
                        }
                    }

                }
            }

            // Update the cluster Id
            if(updateClusterID) {
                updateClusterID = false;
                clusterID++;
            }
        }
    }

private:

    ModelPart & mrModelPart;

};


}