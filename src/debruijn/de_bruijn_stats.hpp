/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * de_bruijn_stats.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DE_BRUIJN_STATS_HPP_
#define DE_BRUIJN_STATS_HPP_

#include <vector>
#include <mxx/comm.hpp>
#include <mxx/reduction.hpp>

namespace bliss {

  namespace de_bruijn {

    template <typename Kmer, typename Edge>
    void print_dbg_edge_histogram(::std::vector<::std::pair<Kmer, Edge> > const & nodes, mxx::comm const & comm) {

      // local computation
      // 5 possibilities for number of in or out edges
      std::vector<size_t> histogram((Edge::maxEdgeCount + 1) * (Edge::maxEdgeCount + 1), 0UL);   // in count x out count
      size_t x, y;

      for (auto t : nodes) {
        y = t.second.get_in_edge_count();
        x = t.second.get_out_edge_count();

        ++histogram[y * (Edge::maxEdgeCount + 1) + x];
      }

      // finally, print
      if (comm.rank() == 0) {
        printf("LOCAL Edge Existence Histogram: \n");
        for (size_t j = 0; j <= Edge::maxEdgeCount; ++j) {
          printf("\t%lu", j);
        }
        printf("\n");

        for (size_t j = 0; j <= Edge::maxEdgeCount; ++j) {
          printf("%lu", j);
          for (size_t k = 0; k <= Edge::maxEdgeCount; ++k) {
            printf("\t%ld", histogram[j * (Edge::maxEdgeCount + 1) + k]);
          }
          printf("\n");
        }
      }

      // then global reduction
      std::vector<size_t> complete = ::mxx::reduce(histogram, 0, comm);



      // finally, print
      if (comm.rank() == 0) {
        printf("TOTAL Edge Existence Histogram: \n");
        for (size_t j = 0; j <= Edge::maxEdgeCount; ++j) {
          printf("\t%lu", j);
        }
        printf("\n");

        for (size_t j = 0; j <= Edge::maxEdgeCount; ++j) {
          printf("%lu", j);
          for (size_t k = 0; k <= Edge::maxEdgeCount; ++k) {
            printf("\t%ld", complete[j * (Edge::maxEdgeCount + 1) + k]);
          }
          printf("\n");
        }
      }
    }



  } // namespace de_bruijn

} // namespace bliss

#endif // DE_BRUIJN_STATS_HPP_
