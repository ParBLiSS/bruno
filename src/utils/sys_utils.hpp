/*
 * Copyright 2017 Georgia Institute of Technology
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

/**
 * @file    sys_utils.h
 * @ingroup utils
 * @brief   Defines helper utilities for getting various system properties.  linux only for now.
 * @author  Tony Pan <tpan7@gatech.edu>
 */

#ifndef SYS_UTILS_H
#define SYS_UTILS_H

#include <sys/types.h>
#include <sys/sysinfo.h>
#include <unistd.h> // for sysconf

#include <mxx/comm.hpp>
#include <mxx/reduction.hpp>

namespace utils
{

inline unsigned long get_free_mem() {
	struct sysinfo memInfo;
	sysinfo(&memInfo);

	return memInfo.freeram * memInfo.mem_unit;
}

inline long get_num_cpus() {
	return sysconf(_SC_NPROCESSORS_CONF);
}


inline unsigned long get_free_mem_per_proc(const mxx::comm & comm) {
	unsigned long free_mem = utils::get_free_mem();
	int local_procs = 1;
	{
		mxx::comm shared_comm = comm.split_shared();
		local_procs = shared_comm.size();
	}
	if (local_procs < 1) throw std::logic_error("evaluated number of local processes on same node to be 0 or less.");
	if (free_mem < 1) throw std::logic_error("evaluated free memory on the node to be 0 or less.");

	size_t mem_per_p = (free_mem / local_procs);  // number of elements that can be held in freemem
	// find the minimum across all processors.
	mem_per_p = mxx::allreduce(mem_per_p, mxx::min<size_t>(), comm);

	if (comm.rank() == 0) std::cout << "estimate available mem=" << free_mem << " bytes, p=" << local_procs <<
			", alloc " << mem_per_p << " elements" << std::endl;

	return mem_per_p;
}

inline unsigned long get_free_mem_per_proc() {
	unsigned long free_mem = utils::get_free_mem();
	int local_procs = utils::get_num_cpus();

	if (local_procs < 1) throw std::logic_error("evaluated number of local processes on same node to be 0 or less.");
	if (free_mem < 1) throw std::logic_error("evaluated free memory on the node to be 0 or less.");

	size_t mem_per_p = (free_mem / local_procs);  // number of elements that can be held in freemem
	// find the minimum across all processors.

//	std::cout << "estimate available mem=" << free_mem << " bytes, p=" << local_procs <<
//			", alloc " << mem_per_p << " elements" << std::endl;
	return mem_per_p;
}



}// namespace utils

#endif // SYS_UTILS_H
