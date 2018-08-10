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
template<class VariableType>
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

    void CalculateLocalClustering(Variable<int> NODEDOMAIN, VariableType NODEDATA) {

        Node<3>::WeakPointer originNodePtr = *(mrModelPart.NodesArray().begin());

        std::list<Node<3>::WeakPointer> modelNodes;         // List of nodes in the queue
        std::list<Node<3>::WeakPointer> componentNodes;     // List of nodes in the current set

        // Locate the elements for every node 
        auto findNeighbourProcesses = FindNodalNeighboursProcess(mrModelPart);
        findNeighbourProcesses.Execute();

        // Push the first node to the nodes queue
        modelNodes.push_back(originNodePtr);

        std::size_t clusterID = 1;
        std::size_t visitId = 1;
        std::size_t nodeCluster = 0;

        bool updateClusterID = false;

        // Guess the size?
        std::vector<double> cluster_info(100, 0);
        std::vector<int> cluster_size(100, 0);

        // Use a traversing hash to avoid marking nodes every repetition
        int traverse_fingerptint = rand();

        while(!modelNodes.empty()) {
            auto & itModelNodePtr = modelNodes.front();
            modelNodes.pop_front();
            int fingerprint = itModelNodePtr.lock()->FastGetSolutionStepValue(FINGERPRINT);
            if(fingerprint != traverse_fingerptint) {
                std::cout << "New non-empty component" << std::endl;
                componentNodes.push_back(itModelNodePtr);
                while(!componentNodes.empty()) {
                    auto & itCompNodePtr = componentNodes.front();
                    componentNodes.pop_front();
                    itCompNodePtr.lock()->FastGetSolutionStepValue(FINGERPRINT) = traverse_fingerptint;

                    int currentDomainType = itCompNodePtr.lock()->FastGetSolutionStepValue(NODEDOMAIN);

                    // IsSolid?
                    if(currentDomainType != 1) {
                        nodeCluster = clusterID;
                        updateClusterID = true;
                    } else {
                        nodeCluster = 0;
                    }

                    itCompNodePtr.lock()->FastGetSolutionStepValue(CLUSTER_ID) = nodeCluster;

                    cluster_info[nodeCluster] += itCompNodePtr.lock()->FastGetSolutionStepValue(NODEDATA);
                    cluster_size[nodeCluster] += 1;

                    auto & neighbourNodes = itCompNodePtr.lock()->GetValue(NEIGHBOUR_NODES);

                    for(auto neigItr = neighbourNodes.ptr_begin(); neigItr != neighbourNodes.ptr_end(); neigItr++) {
                        if((*neigItr).lock()->FastGetSolutionStepValue(FINGERPRINT) != traverse_fingerptint) {
                            int neighbourDomainType = (*neigItr).lock()->FastGetSolutionStepValue(NODEDOMAIN);

                            // If the node is in another cluster add it to the model queue, otherwise add it to the set queue
                            if(currentDomainType != neighbourDomainType) {
                                modelNodes.push_back((*neigItr));
                            } else {
                                (*neigItr).lock()->FastGetSolutionStepValue(FINGERPRINT) = traverse_fingerptint;
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
                cluster_info.resize(clusterID);
                cluster_size.resize(clusterID);
                cluster_info[clusterID-1] = 0.0;
                cluster_info[clusterID-1] = 0;
            }
        }

        mLocalClusterNum = clusterID;
    }

    void Execute() {
        
        CalculateLocalClustering(STATIONARY, NODAL_VOLUME);
    }

private:

    ModelPart & mrModelPart;

    int mTotalClusterNum;
    int mLocalClusterNum;

};


}