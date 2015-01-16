
#include <mpi.h>
#include <omp.h>
#include <unistd.h> // for sleep!

#include <iostream>
#include <functional>

#include <io/CommunicationLayer.hpp>

//#define DEBUG(msg) std::cerr << msg << std::endl;

int my_rank;
std::atomic<int> msgs_received(0);
std::atomic<int> lookup_received(0);
std::atomic<int> answers_received(0);
int elems;

template <bool ThreadLocal = false>
struct Tester
{
  const int ANSWER_TAG = 12;
  const int FIRST_TAG = 1;
  const int LOOKUP_TAG = 13;

  int generate_message(int srcRank, int dstRank)
  {
    return (srcRank + 1) * 100000 + (dstRank + 1);
  }

  void receivedCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    //DEBUG("Rank " << my_rank << " received " << count << " message from process: " << fromRank);

    // first: deserialize
    int* msgs = reinterpret_cast<int*>(msg);
    int msg_count = count / sizeof(int);

    // for all received requests, send the value from the lookup
    for (int i = 0; i < msg_count; ++i)
    {
      // check that the message is as expected
      if(msgs[i] != generate_message(fromRank, my_rank))
      {
//        ERROR("ERROR: message not as expected.  Expected: " << generate_message(fromRank, my_rank) << " Actual: "<< msgs[i] << "");
//        ERROR("my rank: " << my_rank << " from rank " << fromRank);
        std::cerr << "ERROR: message not as expected.  Expected: " << generate_message(fromRank, my_rank) << " Actual: "<< msgs[i] << "" << std::endl;
        std::cerr << "my rank: " << my_rank << " from rank " << fromRank << " message id = " << msgs_received.load() << std::endl;

        exit(EXIT_FAILURE);
      }
      else
      {
        // DEBUG("SUCCESS: message received");
        msgs_received.fetch_add(1);
      }
    }
  }

  void lookup_callback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // first: deserialize
    int* msgs = reinterpret_cast<int*>(msg);
    int msg_count = count / sizeof(int);

    // for all received requests, send the value from the lookup
    for (int i = 0; i < msg_count; ++i)
    {
      // check that the message is as expected
      if(msgs[i] != (fromRank+1)* (my_rank+1))
      {
        //ERROR("ERROR: LOOKUP message not as expected: " << msgs[i]);
        std::cerr << "ERROR: LOOKUP message not as expected: " << msgs[i] << std::endl;
        exit(EXIT_FAILURE);
      } else {
        // DEBUG("SUCCESS: message received");
        lookup_received.fetch_add(1);

      }
      int msg = msgs[i] + 13;
      commLayer.sendMessage(&msg, sizeof(int), fromRank, ANSWER_TAG);
    }
  }

  void answer_callback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // first: deserialize
    int* msgs = reinterpret_cast<int*>(msg);
    int msg_count = count / sizeof(int);

    for (int i = 0; i < msg_count; ++i)
    {
      // check that the message is as expected
      if(msgs[i] != (fromRank+1)* (my_rank+1) + 13)
      {
        //ERROR("ERROR: ANSWER message not as expected: " << msgs[i]);
        std::cerr << "ERROR: ANSWER message not as expected. val= " << msgs[i] << " expected = " << ((fromRank+1)* (my_rank+1) + 13) << std::endl;
        exit(EXIT_FAILURE);
      }
      answers_received.fetch_add(1);
    }
  }

  void test_comm_layer(int repeat_sends=10000)
  {
    DEBUG("Testing Comm Layer");
    DEBUG("Size: " << commLayer.getCommSize());
    DEBUG("Rank: " << commLayer.getCommRank());

    // set global rank
    my_rank = commLayer.getCommRank();
    using namespace std::placeholders;
    commLayer.addReceiveCallback(FIRST_TAG, std::bind(&Tester::receivedCallback, this, _1, _2, _3));
    commLayer.addReceiveCallback(LOOKUP_TAG, std::bind(&Tester::lookup_callback, this, _1, _2, _3));
    commLayer.addReceiveCallback(ANSWER_TAG, std::bind(&Tester::answer_callback, this, _1, _2, _3));

    commLayer.initCommunication();

    int iters = 2;
    int nthreads = numThreads;

    for (int it = 0; it < iters; ++it) {

      // start sending one message to each:
#pragma omp parallel for default(none) num_threads(nthreads) shared(repeat_sends, my_rank, it, stdout)
      for (int i = 0; i < repeat_sends; ++i)
      {

        for (int j = 0; j < commSize; ++j)
        {
          int msg = generate_message(my_rank, j);
          if (i % 1000 == 0) DEBUG("Iter " << it << " Rank " << my_rank << "." << omp_get_thread_num() << ", " << i << " of " << repeat_sends << " Sending " << msg << " to Rank " << j);
          commLayer.sendMessage(&msg, sizeof(int), j, FIRST_TAG);
        }
      }

//      if (commLayer.getCommRank() == 0) {
//        sleep(1);
//      }

      DEBUG("Iter " << it << " Rank " << commRank << " thread " << omp_get_thread_num() << " FLUSH. so far messages received = " << msgs_received.load());

      // call the flush function for this tag
      commLayer.flush(FIRST_TAG);

      //=== debug messages show that there are no messages waiting.  so where are the missing messages?
      DEBUG("Iter " << it << " Rank " << commRank << " thread " << omp_get_thread_num() << " flushed. messages received = " << msgs_received.load());


    }
    // call the finish function for this tag  //
    commLayer.finish(FIRST_TAG);
    INFO("Rank " << commRank << " COMM LAYER TAG 1 FINISHED.\n")
    // check that all messages have been received
    if (msgs_received.load() != repeat_sends * commSize * iters)
    {
//      ERROR("ERROR: wrong amount of messages received in phase 1");
//      ERROR("received: " << msgs_received.load() << ", should: " << repeat_sends * commLayer.getCommSize() * iters);
      std::cerr << "rank " << commRank << " FAIL: wrong amount of messages received in phase 1" << std::endl;
      std::cerr << "rank " << commRank << " received: " << msgs_received.load() << ", should be: " << repeat_sends * commSize * iters << std::endl;
      exit(EXIT_FAILURE);
    }
    //std::cerr << "INDEX: " << msgs_received << std::endl;


    /* phase 2 communication */

    // sending one message to each:
#pragma omp parallel for default(none) num_threads(nthreads) shared(repeat_sends, my_rank, stdout)
    for (int i = 0; i < repeat_sends; ++i)
    {

      for (int j = 0; j < commSize; ++j)
      {
        int msg = generate_message(my_rank, j);
        if (i % 1000 == 0) DEBUG("Rank " << my_rank << "." << omp_get_thread_num() << " " << i << " of " << repeat_sends << " Querying " << msg);
        commLayer.sendMessage(&msg, sizeof(int), j, LOOKUP_TAG);
      }
    }

    // flush both tags
    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());
    commLayer.flush(LOOKUP_TAG);
    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());

    // flush both tags
    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());
    commLayer.flush(LOOKUP_TAG);
    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());

    //======== call flush twice helps to avoid error below.  HACK
    // check that all messages have been received correctly
    if (lookup_received != repeat_sends * commSize)
    {
      std::cerr << "rank " << commRank << " FAIL: wrong amount of lookup messages received in phase 2" << std::endl;
      std::cerr << "rank " << commRank << " received: " << lookup_received.load() << ", should: " << repeat_sends * commSize << std::endl;
      exit(EXIT_FAILURE);
    }

    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());
    commLayer.flush(ANSWER_TAG);
    DEBUG("rank " << commRank << "thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());

    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());
    commLayer.flush(ANSWER_TAG);
    DEBUG("rank " << commRank << "thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());

    //======== call flush twice helps to avoid error below.  HACK
    // check that all messages have been received correctly
    if (answers_received != repeat_sends * commSize)
    {
      std::cerr << "rank " << commRank << " FAIL: wrong amount of answer messages received in phase 2" << std::endl;
      std::cerr << "rank " << commRank << " received: " << answers_received.load() << ", should: " << repeat_sends * commSize << std::endl;
      exit(EXIT_FAILURE);
    }

    INFO("Rank " << commRank << " DONE DONE DONE. ");



//    commLayer.finish(FIRST_TAG);
//    commLayer.finishTag(LOOKUP_TAG);
//    commLayer.finishTag(ANSWER_TAG);
    commLayer.finishCommunication();

    //std::cerr << "LOOKUP: " << lookup_received << " ANSWERS: " << answers_received << std::endl;

    DEBUG("This was a triumph.");
//    sleep(1);
    DEBUG("I'm making a note here: HUGE SUCCESS.");
//    sleep(1);
    DEBUG("It's hard to overstate my satisfaction.");

    INFO("DONE.");
  }

  Tester(MPI_Comm comm, int comm_size, int num_threads) :
    commLayer(comm, comm_size, num_threads), commSize(comm_size), numThreads(num_threads) {
    //commLayer.startThreads();
    MPI_Comm_rank(comm, &commRank);
  }

  bliss::io::CommunicationLayer<ThreadLocal> commLayer;

  int commSize;
  int commRank;
  int numThreads;
};

int main(int argc, char *argv[])
{
  int nthreads = 1;
  if (argc > 1) {
    nthreads = atoi(argv[1]);
  }

  elems = 1536 * nthreads;
  if (argc > 2) {
    elems = atoi(argv[2]);
  }

  // set up MPI
  MPI_Init(&argc, &argv);

  // get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int p, rank;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);

  /* code */
  {
    msgs_received.store(0);
    lookup_received.store(0);
    answers_received.store(0);

#if defined(THREADLOCAL)
    Tester<true> tester(comm, p, nthreads);
#else
    Tester<false> tester(comm, p, nthreads);
#endif
    tester.test_comm_layer(elems);

    MPI_Barrier(comm);
  }

  // finalize MPI
  MPI_Finalize();
  return 0;
}

