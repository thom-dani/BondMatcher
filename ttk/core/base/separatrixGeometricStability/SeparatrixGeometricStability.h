/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::SeparatrixGeometricStability
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %SeparatrixGeometricStability class that computes
/// for each vertex of a triangulation the average scalar value of itself and
/// its direct neighbors.
///
/// \b Related \b publication: \n
/// 'SeparatrixGeometricStability'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The SeparatrixGeometricStability class provides methods to compute for each
   * vertex of a triangulation the average scalar value of itself and its direct
   * neighbors.
   */
  class SeparatrixGeometricStability : virtual public Debug {

  public:
    SeparatrixGeometricStability();

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeAverages(dataType *outputData,
                        const dataType *inputData,
                        const triangulationType *triangulation) const {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Vertex Averages
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing Averages",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

        // compute the average of each vertex in parallel
        size_t const nVertices = triangulation->getNumberOfVertices();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(size_t i = 0; i < nVertices; i++) {
          // initialize average
          outputData[i] = inputData[i];

          // add neighbor values to average
          size_t const nNeighbors = triangulation->getVertexNeighborNumber(i);
          ttk::SimplexId neighborId{-1};
          for(size_t j = 0; j < nNeighbors; j++) {
            triangulation->getVertexNeighbor(i, j, neighborId);
            outputData[i] += inputData[neighborId];
          }

          // divide by neighbor number
          outputData[i] /= (nNeighbors + 1);
        }

        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computing Averages",
                       1, // progress
                       localTimer.getElapsedTime(), this->threadNumber_);
      }

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

  }; // SeparatrixGeometricStability class

} // namespace ttk
