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
 * debruijn_stats.hpp
 *
 *  Created on: June 10, 2016
 *      Author: Tony Pan
 */

#ifndef DEBRUIJN_STATS_HPP_
#define DEBRUIJN_STATS_HPP_

#include <vector>
#include <utility>   // pair

#include <mxx/comm.hpp>
#include <mxx/reduction.hpp>

namespace bliss {

  namespace debruijn {

    namespace graph {

      template <typename Iterator>
      void print_compact_multi_biedge_histogram(Iterator first, Iterator last, mxx::comm const & comm) {

    	  using Edge = typename std::iterator_traits<Iterator>::value_type::second_type;

        std::vector<size_t> complete;
        {
          // local computation
          // 5 possibilities for number of in or out edges
          std::vector<size_t> histogram((Edge::maxEdgeCount + 1) * (Edge::maxEdgeCount + 1), 0UL);   // in count x out count

          for (auto it = first; it != last; ++it) {
            ++histogram[it->second.get_in_edge_count() * (Edge::maxEdgeCount + 1) + it->second.get_out_edge_count()];
          }

          // then global reduction
          ::mxx::reduce(histogram, 0, comm).swap(complete);
        }

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

    } // ns graph

  } // namespace debruijn

} // namespace bliss

#endif // DEBRUIJN_STATS_HPP_
