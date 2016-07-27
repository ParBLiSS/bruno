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

/**
 * @file    kmer_index.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief	High level kmer indexing API
 * @details	3 primary Kmer Index classes are currently provided:
 * 			Kmer Count Index,
 * 			Kmer Position Index, and
 * 			Kmer Position + Quality score index.
 *
 *
 *		Currently, KmerIndex classes only supports FASTQ files, but the intent is to
 *		  allow replaceable module for file reading.
 *
 *		Data distribution is via replaceable Hash functions.  currently,
 *		  PrefixHasher is used to distribute to MPI processes based on the first log_2(P) bits of the kmer
 *		  InfixHasher is used to distribute to threads within a process, based on the next log_2(p) bits of the kmer
 *		  SuffixHasher is used for threadlocal unordered_map hashing
 *
 *		All are deterministic to allow simple lookup processes.
 *
 */
#ifndef KMER_INDEX_HPP_
#define KMER_INDEX_HPP_

#include "bliss-config.hpp"

#if defined(USE_MPI)
#include "mpi.h"
#endif

//#if defined(USE_OPENMP)
//#include "omp.h"
//#endif

#include <unistd.h>     // sysconf
#include <sys/stat.h>   // block size.
#include <tuple>        // tuple and utility functions
#include <utility>      // pair and utility functions.
#include <type_traits>
#include <cctype>       // tolower.

#include "io/file.hpp"
#include "io/fastq_loader.hpp"
#include "io/fasta_loader.hpp"
//#include "io/fasta_iterator.hpp"

#include "utils/logging.h"
#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "utils/kmer_utils.hpp"
#include "index/kmer_hash.hpp"
#include "common/kmer_transform.hpp"

#include "io/mxx_support.hpp"
//#include "containers/distributed_hashed_vec.hpp"
#include "containers/distributed_unordered_map.hpp"
#include "containers/distributed_sorted_map.hpp"
//#include "containers/distributed_map.hpp"
#include "containers/distributed_densehash_map.hpp"

#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"
#include "iterators/transform_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "iterators/zip_iterator.hpp"
#include "iterators/unzip_iterator.hpp"
#include "iterators/constant_iterator.hpp"
#include "index/quality_score_iterator.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/file_utils.hpp"

#include <fstream> // debug only
#include <iostream>  // debug only
#include <sstream>  // debug only

namespace bliss
{
namespace index
{
namespace kmer
{

/**
 * @tparam MapType  	container type
 * @tparam KmerParser		functor to generate kmer (tuple) from input.  specified here so we specialize for different index.  note KmerParser needs to be supplied with a data type.
 */
template <typename MapType, typename KmerParser>
class Index {
protected:

	MapType map;

	const mxx::comm& comm;

public:
	using KmerType = typename MapType::key_type;
	using ValueType   = typename MapType::mapped_type;
	using TupleType =      std::pair<KmerType, ValueType>;
	using Alphabet = typename KmerType::KmerAlphabet;

	using KmerParserType = KmerParser;

	Index(const mxx::comm& _comm) : map(_comm), comm(_comm) {
	}

	virtual ~Index() {};

	MapType & get_map() {
		return map;
	}
	MapType const & get_map() const {
		return map;
	}


	/**
	 * @brief  generate kmers or kmer tuples for 1 block of raw data.
	 * @note   requires that SeqParser be passed in and operates on the Block's Iterators.
	 *          Mostly, this is because we need to broadcast the state of SeqParser to procs on the same node (L1 seq info) and recreate on child procs a new SeqParser (for L2)
	 * @tparam SeqParser		parser type for extracting sequences.  supports FASTQ and FASTA.  template template parameter, param is iterator
	 * @tparam KP           parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
	 * @tparam BlockType		input partition type, supports in memory (vector) vs memmapped.
	 * @param partition
	 * @param result        output vector.  should be pre allocated.
	 */
	template <typename KP, template <typename> class SeqParser, typename BlockType>
	size_t read_block(BlockType const & partition,
			SeqParser<typename BlockType::iterator> const &seq_parser,
			std::vector<typename KP::value_type>& result) const {

		// from FileLoader type, get the block iter type and range type
		using BlockIterType = typename BlockType::const_iterator;

		using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, SeqParser >;
		//		using SeqType = typename ::std::iterator_traits<SeqIterType>::value_type;

		//== sequence parser type
		KP kmer_parser;

		//== process the chunk of data

		//==  and wrap the chunk inside an iterator that emits Reads.
		SeqIterType seqs_start(seq_parser, partition.cbegin(), partition.in_mem_cend(), partition.getRange().start);
		SeqIterType seqs_end(partition.in_mem_cend());

		::fsc::back_emplace_iterator<std::vector<typename KP::value_type> > emplace_iter(result);

		size_t before = result.size();

		//== loop over the reads
		for (; seqs_start != seqs_end; ++seqs_start)
		{
			if ((*seqs_start).seq_size() == 0) continue;
			//		  std::cout << "** seq: " << (*seqs_start).id.id << ", ";
			//		  ostream_iterator<typename std::iterator_traits<typename SeqType::IteratorType>::value_type> osi(std::cout);
			//		  std::copy((*seqs_start).seq_begin, (*seqs_start).seq_end, osi);
			//		  std::cout << std::endl;

			emplace_iter = kmer_parser(*seqs_start, emplace_iter);

			//	    std::cout << "Last: pos - kmer " << result.back() << std::endl;

		}

		return result.size() - before;
	}


	template <template <typename> class SeqParser, typename KP, typename BlockType>
	size_t parse_file_data(const BlockType & partition,
	                       std::vector<typename KP::value_type>& result, const mxx::comm & _comm) const {
		 size_t before = result.size();

	    BL_BENCH_INIT(file);
	    {
	      // not reusing the SeqParser in loader.  instead, reinitializing one.
	      BL_BENCH_START(file);
	      SeqParser<typename BlockType::const_iterator> l1parser;
	      l1parser.init_parser(partition.in_mem_cbegin(), partition.parent_range_bytes, partition.in_mem_range_bytes, partition.getRange(), _comm);
	      BL_BENCH_END(file, "mark_seqs", partition.getRange().size());

	      //== reserve
	      BL_BENCH_START(file);
	      // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
	      // index reserve internally sends a message to itself.

	      // call after getting first L1Block to ensure that file is loaded.  (rank 0 reads and broadcast)
	      size_t record_size = 0;
	      size_t seq_len = 0;
	      std::tie(record_size, seq_len) = l1parser.get_record_size(partition.cbegin(), partition.parent_range_bytes, partition.getRange(), partition.getRange(), _comm, 10);
	      size_t est_size = (record_size == 0) ? 0 : (partition.getRange().size() + record_size - 1) / record_size;  // number of records
	      est_size *= (seq_len < KmerType::size) ? 0 : (seq_len - KmerType::size + 1) ;  // number of kmers in a record
	      result.reserve(est_size);
	      BL_BENCH_END(file, "reserve", est_size);

	      BL_BENCH_START(file);
	      //=== copy into array
	      if (partition.getRange().size() > 0) {
	        read_block<KP, SeqParser>(partition, l1parser, result);
	      }
	      BL_BENCH_END(file, "read", result.size());
	      // std::cout << "Last: pos - kmer " << result.back() << std::endl;
	    }

	    BL_BENCH_REPORT_MPI_NAMED(file, "index:read_file_data", _comm);
	    return result.size() - before;

	}

	template <typename FileType>
	::bliss::io::file_data open_file(const std::string & filename, const mxx::comm & _comm) const {
	      // file extension determines SeqParserType
	      std::string extension = ::bliss::utils::file::get_file_extension(filename);
	      std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
	      if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0) && (extension.compare("fa") != 0)) {
	        throw std::invalid_argument("input filename extension is not supported.");
	      }

	      ::bliss::io::file_data partition;

		    BL_BENCH_INIT(file);
		    {  // ensure that fileloader is closed at the end.

		      BL_BENCH_START(file);
		      //==== create file Loader
		//      FileLoaderType loader(filename, _comm, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
		//      typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();

		      FileType fobj(filename, KmerType::size - 1, _comm);
		      partition = fobj.read_file();
		      BL_BENCH_END(file, "open", partition.getRange().size());
		    }
		      BL_BENCH_REPORT_MPI_NAMED(file, "index:open_file", _comm);
		return partition;
	}


	/**
	 * @brief read a file's content and generate kmers, place in a vector as return result.
	 * @note  static so can be used wihtout instantiating a internal map.
	 * @tparam SeqParser		parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
	 * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
	 */
	template <template <typename> class SeqParser, typename KP>
	size_t read_file_mpiio(const std::string & filename,
	                       std::vector<typename KP::value_type>& result,
	                       const mxx::comm & _comm) const {

	    size_t before = result.size();

	    using FileType = ::bliss::io::parallel::mpiio_file<SeqParser > ;

	    BL_BENCH_INIT(file);
	    {  // ensure that fileloader is closed at the end.

	      BL_BENCH_START(file);
	      ::bliss::io::file_data partition = open_file<FileType>(filename, _comm);
	      BL_BENCH_END(file, "open", partition.getRange().size());

	      //std::cout << "rank " << _comm.rank() << " mpiio " << partition.getRange() << " in mem " << partition.in_mem_range_bytes << std::endl;


	      // not reusing the SeqParser in loader.  instead, reinitializing one.
	      BL_BENCH_START(file);
	      parse_file_data<SeqParser, KP>(partition, result, _comm);
	      BL_BENCH_END(file, "read", result.size());
	      // std::cout << "Last: pos - kmer " << result.back() << std::endl;
	    }


      BL_BENCH_REPORT_MPI_NAMED(file, "index:read:mpiio", _comm);
	    return result.size() - before;
	}


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <template <typename> class SeqParser, typename KP>
  size_t read_file_mmap(const std::string & filename,
                        std::vector<typename KP::value_type>& result,
                        const mxx::comm & _comm) const {

//      // file extension determines SeqParserType
//      std::string extension = ::bliss::utils::file::get_file_extension(filename);
//      std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
//      if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0)) {
//        throw std::invalid_argument("input filename extension is not supported.");
//      }
//
//      // check to make sure that the file parser will work
//      if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
//        throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
//      } else if ((extension.compare("fasta") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
//        throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
//      }


      // partitioned file with mmap or posix do not seem to be much faster than mpiio and may result in more jitter when congested.
      using FileType = ::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, SeqParser >;
	    //====  now process the file, one L1 block (block partition by MPI Rank) at a time

      size_t before = result.size();

      BL_BENCH_INIT(file);
      {  // ensure that fileloader is closed at the end.

        BL_BENCH_START(file);
	      ::bliss::io::file_data partition = open_file<FileType>(filename, _comm);
        BL_BENCH_END(file, "open", partition.getRange().size());


	      // not reusing the SeqParser in loader.  instead, reinitializing one.
	      BL_BENCH_START(file);
	      parse_file_data<SeqParser, KP>(partition, result, _comm);
	      BL_BENCH_END(file, "read", result.size());
	      // std::cout << "Last: pos - kmer " << result.back() << std::endl;
	  }

      BL_BENCH_REPORT_MPI_NAMED(file, "index:read:mmap_file", _comm);
      return result.size() - before;
  }


  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <template <typename> class SeqParser, typename KP>
  size_t read_file_posix(const std::string & filename,
                         std::vector<typename KP::value_type>& result,
                         const mxx::comm & _comm) const {

//      // file extension determines SeqParserType
//      std::string extension = ::bliss::utils::file::get_file_extension(filename);
//      std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
//      if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0)) {
//        throw std::invalid_argument("input filename extension is not supported.");
//      }
//
//      // check to make sure that the file parser will work
//      if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
//        throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
//      } else if ((extension.compare("fasta") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
//        throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
//      }


      // partitioned file with mmap or posix do not seem to be much faster than mpiio and may result in more jitter when congested.
      using FileType = ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, SeqParser >;
	    //====  now process the file, one L1 block (block partition by MPI Rank) at a time

      size_t before = result.size();

      BL_BENCH_INIT(file);
      {  // ensure that fileloader is closed at the end.

        BL_BENCH_START(file);
	      ::bliss::io::file_data partition = open_file<FileType>(filename, _comm);
        BL_BENCH_END(file, "open", partition.getRange().size());

	      // not reusing the SeqParser in loader.  instead, reinitializing one.
	      BL_BENCH_START(file);
	      parse_file_data<SeqParser, KP>(partition, result, _comm);
	      BL_BENCH_END(file, "read", result.size());
	      // std::cout << "Last: pos - kmer " << result.back() << std::endl;
	  }

      BL_BENCH_REPORT_MPI_NAMED(file, "index:read:posix_file", _comm);
      return result.size() - before;
  }





	std::vector<TupleType> find_overlap(std::vector<KmerType> &query) const {
		return map.find_overlap(query);
	}
	std::vector<TupleType> find(std::vector<KmerType> &query) const {
		return map.find(query);
	}
	std::vector<TupleType> find_collective(std::vector<KmerType> &query) const {
		return map.find_collective(query);
	}
  std::vector<TupleType> find_sendrecv(std::vector<KmerType> &query) const {
    return map.find_sendrecv(query);
  }
	std::vector< std::pair<KmerType, size_t> > count(std::vector<KmerType> &query) const {
		return map.count(query);
	}

	void erase(std::vector<KmerType> &query) {
		map.erase(query);
	}


	template <typename Predicate>
	std::vector<TupleType> find_if_overlap(std::vector<KmerType> &query, Predicate const &pred) const {
		return map.find_overlap(query, false, pred);
	}
	template <typename Predicate>
	std::vector<TupleType> find_if(std::vector<KmerType> &query, Predicate const &pred) const {
		return map.find(query, false, pred);
	}
	template <typename Predicate>
	std::vector<TupleType> find_if_collective(std::vector<KmerType> &query, Predicate const &pred) const {
		return map.find_collective(query, false, pred);
	}
  template <typename Predicate>
  std::vector<TupleType> find_if_sendrecv(std::vector<KmerType> &query, Predicate const &pred) const {
    return map.find_sendrecv(query, false, pred);
  }
	template <typename Predicate>
	std::vector<TupleType> find_if(Predicate const &pred) const {
		return map.find(pred);
	}

	template <typename Predicate>
	std::vector<std::pair<KmerType, size_t> > count_if(std::vector<KmerType> &query, Predicate const &pred) const {
		return map.count(query, false, pred);
	}

	template <typename Predicate>
	std::vector<std::pair<KmerType, size_t> > count_if(Predicate const &pred) const {
		return map.count(pred);
	}


	template <typename Predicate>
	void erase_if(std::vector<KmerType> &query, Predicate const &pred) {
		map.erase(query, false, pred);
	}

	template <typename Predicate>
	void erase_if(Predicate const &pred) {
		map.erase(pred);
	}


	/**
	 * @tparam T 	input type may not be same as map's value types, so map need to provide overloads (and potentially with transform operators)
	 */
	 template <typename T>
	void insert(std::vector<T> &temp) {
		BL_BENCH_INIT(insert);

		// do not reserve until insertion - less transient memory used.
//		BL_BENCH_START(build);
//		this->map.reserve(this->map.size() + temp.size());
//		BL_BENCH_END(build, "reserve", temp.size());

		// distribute
		BL_BENCH_START(insert);
		this->map.insert(temp);  // COLLECTIVE CALL...
		BL_BENCH_END(insert, "map_insert", this->map.local_size());

#if (BL_BENCHMARK == 1)
		BL_BENCH_START(insert);
		size_t m = 0;  // here because sortmap needs it.
		m = this->map.get_multiplicity();
		BL_BENCH_END(insert, "multiplicity", m);
#else
		auto result = this->map.get_multiplicity();
		BLISS_UNUSED(result);
#endif
		BL_BENCH_REPORT_MPI_NAMED(insert, "index:insert", this->comm);

	 }

	 // Note that KmerParserType may depend on knowing the Sequence Parser Type (e.g. provide quality score iterators)
	 //	Output type of KmerParserType may not match Map value type, in which case the map needs to do its own transform.
	 //     since Kmer template parameter is not explicitly known, we can't hard code the return types of KmerParserType.

	 //============= THESE ARE TO BE DEPRECATED


	 // Note that KmerParserType may depend on knowing the Sequence Parser Type (e.g. provide quality score iterators)
	 //	Output type of KmerParserType may not match Map value type, in which case the map needs to do its own transform.
	 //     since Kmer template parameter is not explicitly known, we can't hard code the return types of KmerParserType.

	 /// convenience function for building index.
	 template <template <typename> class SeqParser>
	 void build_mpiio(const std::string & filename, MPI_Comm comm) {

		 // file extension determines SeqParserType
		 std::string extension = ::bliss::utils::file::get_file_extension(filename);
		 std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
		 if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0) && (extension.compare("fa") != 0)) {
			 throw std::invalid_argument("input filename extension is not supported.");
		 }

		 // check to make sure that the file parser will work
		 if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
			 throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
		 } else if (((extension.compare("fasta") == 0) || (extension.compare("fa") == 0)) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
			 throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
		 }
     BL_BENCH_INIT(build);

		 // proceed
     BL_BENCH_START(build);
		 ::std::vector<typename KmerParser::value_type> temp;
		 this->read_file_mpiio<SeqParser, KmerParser >(filename, temp, comm);
     BL_BENCH_END(build, "read", temp.size());


		 //        // dump the generated kmers to see if they look okay.
		 //         std::stringstream ss;
		 //         ss << "test." << commRank << ".log";
		 //         std::ofstream ofs(ss.str());
		 //         for (int i = 0; i < temp.size(); ++i) {
		 //          ofs << "item: " << i << " value: " << temp[i] << std::endl;
		 //         }
		 //         ofs.close();

     BL_BENCH_START(build);
		 this->insert(temp);
     BL_BENCH_END(build, "insert", temp.size());


     BL_BENCH_REPORT_MPI_NAMED(build, "index:build_mpiio", this->comm);

	 }


	  /// convenience function for building index.
	   template <template <typename> class SeqParser>
	   void build_mmap(const std::string & filename, MPI_Comm comm) {

	     // file extension determines SeqParserType
	     std::string extension = ::bliss::utils::file::get_file_extension(filename);
	     std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
	     if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0) && (extension.compare("fa") != 0)) {
	       throw std::invalid_argument("input filename extension is not supported.");
	     }

	     // check to make sure that the file parser will work
	     if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
	       throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
	     } else if (((extension.compare("fasta") == 0) || (extension.compare("fa") == 0)) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
	       throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
	     }
	     BL_BENCH_INIT(build);

	     // proceed
	     BL_BENCH_START(build);
	     ::std::vector<typename KmerParser::value_type> temp;
	     this->read_file_mmap<SeqParser, KmerParser >(filename, temp, comm);
	      BL_BENCH_END(build, "read", temp.size());


	     //        // dump the generated kmers to see if they look okay.
	     //         std::stringstream ss;
	     //         ss << "test." << commRank << ".log";
	     //         std::ofstream ofs(ss.str());
	     //         for (int i = 0; i < temp.size(); ++i) {
	     //          ofs << "item: " << i << " value: " << temp[i] << std::endl;
	     //         }
	     //         ofs.close();

	     BL_BENCH_START(build);
	     this->insert(temp);
	      BL_BENCH_END(build, "insert", temp.size());


	     BL_BENCH_REPORT_MPI_NAMED(build, "index:build_mmap", this->comm);

	   }




		 /// convenience function for building index.
		 template <template <typename> class SeqParser>
		 void build_posix(const std::string & filename, MPI_Comm comm) {

			 // file extension determines SeqParserType
			 std::string extension = ::bliss::utils::file::get_file_extension(filename);
			 std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
			 if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0) && (extension.compare("fa") != 0)) {
				 throw std::invalid_argument("input filename extension is not supported.");
			 }

			 // check to make sure that the file parser will work
			 if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
				 throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
			 } else if (((extension.compare("fasta") == 0) || (extension.compare("fa") == 0)) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser<char*> >::value)) {
				 throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
			 }
	     BL_BENCH_INIT(build);

			 // proceed
	     BL_BENCH_START(build);
			 ::std::vector<typename KmerParser::value_type> temp;
			 this->read_file_posix<SeqParser, KmerParser >(filename, temp, comm);
	     BL_BENCH_END(build, "read", temp.size());


			 //        // dump the generated kmers to see if they look okay.
			 //         std::stringstream ss;
			 //         ss << "test." << commRank << ".log";
			 //         std::ofstream ofs(ss.str());
			 //         for (int i = 0; i < temp.size(); ++i) {
			 //          ofs << "item: " << i << " value: " << temp[i] << std::endl;
			 //         }
			 //         ofs.close();

	     BL_BENCH_START(build);
			 this->insert(temp);
	     BL_BENCH_END(build, "insert", temp.size());


	     BL_BENCH_REPORT_MPI_NAMED(build, "index:build_posix", this->comm);

		 }




   typename MapType::const_iterator cbegin() const
   {
     return map.cbegin();
   }

   typename MapType::const_iterator cend() const {
     return map.cend();
   }

   size_t size() const {
     return map.size();
   }

   size_t local_size() const {
     return map.local_size();
   }

};


struct NotEOL {
  template <typename CharType>
  bool operator()(CharType const & x) {
	return (x != '\n') && (x != '\r' );
  }

  template <typename CharType, typename MDType>
  bool operator()(std::pair<CharType, MDType> const & x) {
	  return (x.first != '\n') && (x.first != '\r');
  }
};

template <typename Iter>
using NonEOLIter = bliss::iterator::filter_iterator<NotEOL, Iter>;



/**
 * @tparam KmerType       output value type of this parser.  not necessarily the same as the map's final storage type.
 */
template <typename KmerType>
struct KmerParser {

    /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = KmerType;
	using kmer_type = KmerType;


  /**
   * @brief generate kmers from 1 sequence.  result inserted into output_iter, which may be preallocated.
   * @param read          sequence object, which has pointers to the raw byte array.
   * @param output_iter   output iterator pointing to insertion point for underlying container.
   * @return new position for output_iter
   * @tparam SeqType      type of sequence.  inferred.
   * @tparam OutputIt     output iterator type, inferred.
   */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {

		static_assert(std::is_same<KmerType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
						"output type and output container value type are not the same");

		using Alphabet = typename KmerType::KmerAlphabet;

		/// converter from ascii to alphabet values
		using BaseCharIterator = bliss::iterator::transform_iterator<NonEOLIter<typename SeqType::IteratorType>, bliss::common::ASCII2<Alphabet> >;

		/// kmer generation iterator
		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIterType>::value_type,
				KmerType>::value,
				"input iterator and output iterator's value types differ");

		// then compute and store into index (this will generate kmers and insert into index)
    if (::std::distance(read.seq_begin, read.seq_end) < KmerType::size) return output_iter;

		//== filtering iterator
		NotEOL neol;
		NonEOLIter<typename SeqType::IteratorType> eolstart(neol, read.seq_begin, read.seq_end);
		NonEOLIter<typename SeqType::IteratorType> eolend(neol, read.seq_end);

		//== set up the kmer generating iterators.
		KmerIterType start(BaseCharIterator(eolstart, bliss::common::ASCII2<Alphabet>()), true);
		KmerIterType end(BaseCharIterator(eolend, bliss::common::ASCII2<Alphabet>()), false);

//    printf("First: pos %lu kmer %s\n", read.id.id, bliss::utils::KmerUtils::toASCIIString(*start).c_str());

		return ::std::copy(start, end, output_iter);

	}
};


/**
 * @tparam TupleType       output value type of this parser.  not necessarily the same as the map's final storage type.
 */
template <typename TupleType>
struct KmerPositionTupleParser {

    /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = TupleType;
	using kmer_type = typename ::std::tuple_element<0, value_type>::type;


  /**
   * @brief generate kmer-position pairs from 1 sequence.  result inserted into output_iter, which may be preallocated.
   * @param read          sequence object, which has pointers to the raw byte array.
   * @param output_iter   output iterator pointing to insertion point for underlying container.
   * @return new position for output_iter
   * @tparam SeqType      type of sequence.  inferred.
   * @tparam OutputIt     output iterator type, inferred.
   */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {

		static_assert(std::is_same<TupleType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
						"output type and output container value type are not the same");
		static_assert(::std::tuple_size<TupleType>::value == 2, "kmer-pos-qual index data type should be a pair");

		using KmerType = typename std::tuple_element<0, TupleType>::type;
		using Alphabet = typename KmerType::KmerAlphabet;

		// filter out EOL characters
		using CharIter = NonEOLIter<typename SeqType::IteratorType>;
    // converter from ascii to alphabet values
    using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
    // kmer generation iterator
    using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

		//== next figure out starting positions for the kmers, accounting for EOL char presenses.
		using IdType = typename std::tuple_element<1, TupleType>::type;
		// kmer position iterator type
		using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

    // use zip iterator to tie together the iteration of sequence raw data and id.
    using PairedIter = bliss::iterator::ZipIterator<typename SeqType::IteratorType, IdIterType>;
    using CharPosIter = NonEOLIter<PairedIter>;
		// now use 2 unzip iterators to access the values.  one of them advances.  all subsequent wrapping
		// iterators trickle back to the zip iterator above, again, one of the 2 unzip iterator will call operator++ on the underlying zip iterator
		using IdIter = bliss::iterator::AdvancingUnzipIterator<CharPosIter, 1>;

		// rezip the results
		using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIter, IdIter>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
				TupleType>::value,
				"input iterator and output iterator's value types differ");


		// then compute and store into index (this will generate kmers and insert into index)
		if (::std::distance(read.seq_begin, read.seq_end) < KmerType::size) return output_iter;  // if too short...

    //== set up the kmer generating iterators.
    NotEOL neol;
    KmerIter start(BaseCharIterator(CharIter(neol, read.seq_begin, read.seq_end), bliss::common::ASCII2<Alphabet>()), true);
    KmerIter end(BaseCharIterator(CharIter(neol, read.seq_end), bliss::common::ASCII2<Alphabet>()), false);


		//== set up the position iterators
    IdType seq_begin_id(read.id);
    seq_begin_id += read.seq_begin_offset;  // change id to point to start of sequence (in file coord)
		IdType seq_end_id(seq_begin_id);
		seq_end_id += read.seq_size();

		// tie chars and id together
		PairedIter pp_begin(read.seq_begin, IdIterType(seq_begin_id));
		PairedIter pp_end(read.seq_end, IdIterType(seq_end_id));

    // filter eol
		CharPosIter cp_begin(neol, pp_begin, pp_end);
    CharPosIter cp_end(neol, pp_end);

		// ==== extract new id and rezip iterators
		KmerIndexIterType index_start(start, IdIter(cp_begin));
		KmerIndexIterType index_end(end, IdIter(cp_end));


//		for (; index_start != index_end; ++index_start) {
//		  auto tp = *index_start;
//
//		  printf("TCP id = %lu, pos = %lu, kmer = %s\n", tp.second.get_id(), tp.second.get_pos(), bliss::utils::KmerUtils::toASCIIString(tp.first).c_str());
//
//		  *output_iter = tp;
//
//		}
		return ::std::copy(index_start, index_end, output_iter);
	}

};


/**
 * @tparam TupleType       output value type of this parser.  not necessarily the same as the map's final storage type.
 */
template <typename TupleType, template<typename> class QualityEncoder = bliss::index::Illumina18QualityScoreCodec>
struct KmerPositionQualityTupleParser {

    /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = TupleType;
	using kmer_type = typename ::std::tuple_element<0, value_type>::type;

  /**
   * @brief generate kmer-position-quality pairs from 1 sequence.  result inserted into output_iter, which may be preallocated.
   * @param read          sequence object, which has pointers to the raw byte array.
   * @param output_iter   output iterator pointing to insertion point for underlying container.
   * @return new position for output_iter
   * @tparam SeqType      type of sequence.  inferred.
   * @tparam OutputIt     output iterator type, inferred.
   */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {


		static_assert(SeqType::has_quality(), "Sequence Parser needs to support quality scores");

		static_assert(std::is_same<TupleType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
				"output type and output container value type are not the same");

		static_assert(::std::tuple_size<TupleType>::value == 2, "kmer-pos-qual index data type should be a pair");

		using KmerType = typename std::tuple_element<0, TupleType>::type;
		using Alphabet = typename KmerType::KmerAlphabet;

		static_assert(::std::tuple_size<typename std::tuple_element<1, TupleType>::type>::value == 2, "pos-qual index data type should be a pair");


    // filter out EOL characters
    using CharIter = NonEOLIter<typename SeqType::IteratorType>;
    // converter from ascii to alphabet values
    using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
    // kmer generation iterator
    using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

    //== next figure out starting positions for the kmers, accounting for EOL char presenses.
    using IdType = typename std::tuple_element<0, typename std::tuple_element<1, TupleType>::type >::type;
    // kmer position iterator type
    using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

    // use zip iterator to tie together the iteration of sequence raw data and id.
    using PairedIter = bliss::iterator::ZipIterator<typename SeqType::IteratorType, IdIterType>;
    using CharPosIter = NonEOLIter<PairedIter>;
    // now use 2 unzip iterators to access the values.  one of them advances.  all subsequent wrapping
    // iterators trickle back to the zip iterator above, again, one of the 2 unzip iterator will call operator++ on the underlying zip iterator
    using IdIter = bliss::iterator::AdvancingUnzipIterator<CharPosIter, 1>;


		using QualType = typename std::tuple_element<1, typename std::tuple_element<1, TupleType>::type>::type;
		//static_assert(::std::is_same<typename SeqType::IdType, IdType>::value, "position type does not match for input and output iterators" );

		// also remove eol from quality score
		using QualIterType =
				bliss::index::QualityScoreGenerationIterator<NonEOLIter<typename SeqType::IteratorType>, KmerType::size, QualityEncoder<QualType> >;

		/// combine kmer iterator and position iterator to create an index iterator type.
		using KmerInfoIterType = bliss::iterator::ZipIterator<IdIter, QualIterType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerInfoIterType>::value_type,
				typename std::tuple_element<1, TupleType>::type >::value,
				"kmer info input iterator and output iterator's value types differ");


		using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIter, KmerInfoIterType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
				TupleType>::value,
				"input iterator and output iterator's value types differ");


		// then compute and store into index (this will generate kmers and insert into index)
		if ((::std::distance(read.seq_begin, read.seq_end) < KmerType::size) ||
		    (::std::distance(read.qual_begin, read.qual_end) < KmerType::size)) return output_iter;
		assert(::std::distance(read.seq_begin, read.seq_end) <= ::std::distance(read.qual_begin, read.qual_end));


    //== set up the kmer generating iterators.
    NotEOL neol;
    KmerIter start(BaseCharIterator(CharIter(neol, read.seq_begin, read.seq_end), bliss::common::ASCII2<Alphabet>()), true);
    KmerIter end(BaseCharIterator(CharIter(neol, read.seq_end), bliss::common::ASCII2<Alphabet>()), false);


    //== set up the position iterators
    IdType seq_begin_id(read.id);
    seq_begin_id += read.seq_begin_offset;  // change id to point to start of sequence (in file coord)
    IdType seq_end_id(seq_begin_id);
    seq_end_id += read.seq_size();

    // tie chars and id together
    PairedIter pp_begin(read.seq_begin, IdIterType(seq_begin_id));
    PairedIter pp_end(read.seq_end, IdIterType(seq_end_id));

    // filter eol
    CharPosIter cp_begin(neol, pp_begin, pp_end);
    CharPosIter cp_end(neol, pp_end);

    // ==== quality scoring
		// filter eol and generate quality scores
		QualIterType qual_start(CharIter(neol, read.qual_begin, read.qual_end));
		QualIterType qual_end(CharIter(neol, read.qual_end));

		KmerInfoIterType info_start(IdIter(cp_begin), qual_start);
		KmerInfoIterType info_end(IdIter(cp_end), qual_end);


		// ==== set up the zip iterators
		KmerIndexIterType index_start(start, info_start);
		KmerIndexIterType index_end(end, info_end);

		//    printf("First: pos %lu kmer %s\n", read.id.id, bliss::utils::KmerUtils::toASCIIString(*start).c_str());

		return ::std::copy(index_start, index_end, output_iter);
	}
};



/**
 * @tparam TupleType       output value type of this parser.  not necessarily the same as the map's final storage type.
 */
template <typename TupleType>
struct KmerCountTupleParser {

  /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = TupleType;
	using kmer_type = typename ::std::tuple_element<0, value_type>::type;

	/**
	 * @brief generate kmer-count pairs from 1 sequence.  result inserted into output_iter, which may be preallocated.
	 * @param read          sequence object, which has pointers to the raw byte array.
	 * @param output_iter   output iterator pointing to insertion point for underlying container.
	 * @return new position for output_iter
	 * @tparam SeqType      type of sequence.  inferred.
	 * @tparam OutputIt     output iterator type, inferred.
	 */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {

		static_assert(std::is_same<TupleType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
				"output type and output container value type are not the same");
		static_assert(::std::tuple_size<TupleType>::value == 2, "count data type should be a pair");

		using KmerType = typename std::tuple_element<0, TupleType>::type;
		using CountType = typename std::tuple_element<1, TupleType>::type;

		using Alphabet = typename KmerType::KmerAlphabet;

		/// converter from ascii to alphabet values
		using BaseCharIterator = bliss::iterator::transform_iterator<NonEOLIter<typename SeqType::IteratorType>, bliss::common::ASCII2<Alphabet> >;

		/// kmer generation iterator
		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
		using CountIterType = bliss::iterator::ConstantIterator<CountType>;

		using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, CountIterType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type, std::pair<KmerType, CountType> >::value,
				"count zip iterator not producing the right type.");

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
				TupleType>::value,
				"count: input iterator and output iterator's value types differ");

		CountIterType count_start(1);

		// then compute and store into index (this will generate kmers and insert into index)
    if (::std::distance(read.seq_begin, read.seq_end) < KmerType::size) return output_iter;

		//== set up the kmer generating iterators.
		NotEOL neol;
		NonEOLIter<typename SeqType::IteratorType> eolstart(neol, read.seq_begin, read.seq_end);
		NonEOLIter<typename SeqType::IteratorType> eolend(neol, read.seq_end);

		KmerIterType start(BaseCharIterator(eolstart, bliss::common::ASCII2<Alphabet>()), true);
		KmerIterType end(BaseCharIterator(eolend, bliss::common::ASCII2<Alphabet>()), false);

		KmerIndexIterType istart(start, count_start);
		KmerIndexIterType iend(end, count_start);

		//    printf("First: pos %lu kmer %s\n", read.id.id, bliss::utils::KmerUtils::toASCIIString(*start).c_str());

		return ::std::copy(istart, iend, output_iter);
	}

};

// TODO: the types of Map that is used should be restricted.  (perhaps via map traits)
template <typename MapType>
using KmerIndex = Index<MapType, KmerParser<typename MapType::key_type> >;

template <typename MapType>
using PositionIndex = Index<MapType, KmerPositionTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;

template <typename MapType>
using PositionQualityIndex = Index<MapType, KmerPositionQualityTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;

template <typename MapType>
using CountIndex = Index<MapType, KmerCountTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;
template <typename MapType>
using CountIndex2 = Index<MapType, KmerParser<typename MapType::key_type> >;

// template aliases for hash to be used as distribution hash
template <typename Key>
using DistHashFarm = ::bliss::kmer::hash::farm<Key, true>;
template <typename Key>
using DistHashMurmur = ::bliss::kmer::hash::murmur<Key, true>;
template <typename Key>
using DistHashStd = ::bliss::kmer::hash::cpp_std<Key, true>;
template <typename Key>
using DistHashIdentity = ::bliss::kmer::hash::identity<Key, true>;


template <typename Key>
using StoreHashFarm = ::bliss::kmer::hash::farm<Key, false>;
template <typename Key>
using StoreHashMurmur = ::bliss::kmer::hash::murmur<Key, false>;
template <typename Key>
using StoreHashStd = ::bliss::kmer::hash::cpp_std<Key, false>;
template <typename Key>
using StoreHashIdentity = ::bliss::kmer::hash::identity<Key, false>;

// =================  Partially defined aliases for MapParams, for distributed_xxx_maps.
// NOTE: when using this, need to further alias so that only Key param remains.
// =================
template <typename Key,
			template <typename> class DistHash  = DistHashMurmur,
			template <typename> class StoreHash = StoreHashMurmur,
			template <typename> class DistTrans = ::bliss::kmer::transform::identity
			>
using SingleStrandHashMapParams = ::dsc::HashMapParams<
		Key,
		::bliss::kmer::transform::identity,  // precanonalizer
		 DistTrans,  				// could be iden, xor, lex_less
		  DistHash,
		  ::std::equal_to,
		   ::bliss::kmer::transform::identity,  // only one that makes sense given InputTransform
		    StoreHash,
		    ::std::equal_to
		  >;


template <typename Key,
	template <typename> class DistHash  = DistHashMurmur,
	template <typename> class StoreHash = StoreHashMurmur
>
using CanonicalHashMapParams = ::dsc::HashMapParams<
		Key,
		::bliss::kmer::transform::lex_less,  // precanonalizer
		 ::bliss::kmer::transform::identity,  // only one that makes sense given InputTransform
		  DistHash,
		  ::std::equal_to,
		   ::bliss::kmer::transform::identity,
		    StoreHash,
		    ::std::equal_to
		  >;

template <typename Key,
	template <typename> class DistHash  = DistHashMurmur,
	template <typename> class StoreHash = StoreHashMurmur
	>
using BimoleculeHashMapParams = ::dsc::HashMapParams<
		Key,
		::bliss::kmer::transform::identity,  // precanonalizer - only one that makes sense for bimole
		 ::bliss::kmer::transform::lex_less,  // only one that makes sense for bimole
		  DistHash,
		  ::std::equal_to,
		   ::bliss::kmer::transform::lex_less,
		    StoreHash,
		    ::std::equal_to
		  >;

//template <typename Key,
//			template <typename> class DistHash = DistHashMurmur,
//			template <typename> class StoreLess = ::std::less,
//			template <typename> class DistTrans = ::bliss::kmer::transform::identity >
//using SingleStrandOrderedMapParams = ::dsc::OrderedMapParams<
//		Key,
//		::bliss::kmer::transform::identity,  // precanonalizer
//		 DistTrans,  				// could be iden, xor, lex_less
//		  DistHash,
//		  ::std::equal_to,
//		   ::bliss::kmer::transform::identity,  // only one that makes sense given InputTransform
//		    StoreLess,
//		    ::std::equal_to
//		  >;
//
//
//template <typename Key,
//			template <typename> class DistHash = DistHashMurmur,
//			template <typename> class StoreLess = ::std::less
//			>
//using CanonicalOrderedMapParams = ::dsc::OrderedMapParams<
//		Key,
//		::bliss::kmer::transform::lex_less,  // precanonalizer
//		 ::bliss::kmer::transform::identity,  // only one that makes sense given InputTransform
//		  DistHash,
//		  ::std::equal_to,
//		   ::bliss::kmer::transform::identity,
//		    StoreLess,
//		    ::std::equal_to
//		  >;
//
//template <typename Key,
//			template <typename> class DistHash = DistHashMurmur,
//			template <typename> class StoreLess = ::std::less
//			>
//using BimoleculeOrderedMapParams = ::dsc::OrderedMapParams<
//		Key,
//		::bliss::kmer::transform::identity,  // precanonalizer - only one that makes sense for bimole
//		 ::bliss::kmer::transform::lex_less,  // only one that makes sense for bimole
//		  DistHash,
//		  ::std::equal_to,
//		   ::bliss::kmer::transform::lex_less,
//		    StoreLess,
//		    ::std::equal_to
//		  >;


template <typename Key,
  template <typename> class Less = ::std::less
>
using SingleStrandSortedMapParams = ::dsc::SortedMapParams<
		Key,
		::bliss::kmer::transform::identity,  // precanonalizer
		   ::bliss::kmer::transform::identity,  // only one that makes sense given InputTransform
		    Less,
		    ::std::equal_to
		  >;


template <typename Key,
			template <typename> class Less = ::std::less
>
using CanonicalSortedMapParams = ::dsc::SortedMapParams<
		Key,
		::bliss::kmer::transform::lex_less,  // precanonalizer
		 ::bliss::kmer::transform::identity,  // only one that makes sense given InputTransform
		    Less,
		    ::std::equal_to
		  >;

template <typename Key,
			template <typename> class Less = ::std::less
>
using BimoleculeSortedMapParams = ::dsc::SortedMapParams<
		Key,
		::bliss::kmer::transform::identity,  // precanonalizer - only one that makes sense for bimole
		 ::bliss::kmer::transform::lex_less,  // only one that makes sense for bimole
		    Less,
		    ::std::equal_to
		  >;
} /* namespace kmer */

} /* namespace index */
} /* namespace bliss */

#endif /* KMER_INDEX_HPP_ */
