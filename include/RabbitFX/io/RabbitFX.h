#ifndef RABBITFX_H
#define RABBITFX_H

#include "./FastxStream.h"
#include "./FastxChunk.h"
#include "assert.h"

struct FM_NEOREF{  //formated neoReference data structure (for format result)
  std::vector<neoReference> vec;
  rabbit::fq::FastqDataChunk *chunk_p;
};

class FQ_SE{
typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> FqChunkQueue;  
typedef rabbit::fq::FastqDataPool FqDataPool;
private:
  rabbit::fq::FastqFileReader* fqFileReader;
  FqDataPool* dp_;
  FqChunkQueue* dq_;
public:
	typedef std::vector<Reference> FDTYPE; //formated type
	typedef FM_NEOREF FDTYPE_NCP;     //formated type no-copy
	std::thread* producer_;
  FQ_SE(std::string file1){
    dq_ = new  FqChunkQueue(128, 1); // data queue
    dp_ = new FqDataPool(256, 1 << 24); // data pool
		bool zipped = false;
		if(ends_with(file1, ".gz")) zipped = true;
    fqFileReader = new rabbit::fq::FastqFileReader(file1, *dp_, zipped);
  }
  void start_producer(){
    producer_ = new std::thread([&] {
      int n_chunks = 0;
      while (true) {
        rabbit::fq::FastqChunk  *fqchunk = new rabbit::fq::FastqChunk;
        fqchunk->chunk = fqFileReader->readNextChunk();
        if (fqchunk->chunk == NULL) break;
        n_chunks++;
        dq_->Push(n_chunks, fqchunk->chunk);
      }

      dq_->SetCompleted();
    });
  }
  FDTYPE get_formated_reads() {
    rabbit::fq::FastqChunk *fqchunk = new rabbit::fq::FastqChunk();  
    rabbit::int64 id = 0;
    int ref_num;
		std::vector<Reference> data;
    if (dq_->Pop(id, fqchunk->chunk)) {
      ref_num = rabbit::fq::chunkFormat(fqchunk->chunk, data, true);
      //------------------relaease-----------------//
      //release_chunk(fqchunk);
			dp_->Release(fqchunk->chunk);
    }else{
      ref_num = 0;
    }
		delete fqchunk;
    return data;
  }
	
  FM_NEOREF get_formated_reads_nocp() {    
    rabbit::fq::FastqChunk *fqchunk = new rabbit::fq::FastqChunk();   //TODO: it will become wild pointer!!!!
    rabbit::int64 id = 0;
    int ref_num = 0;
    FM_NEOREF data;
    if (dq_->Pop(id, fqchunk->chunk)) {
      ref_num = rabbit::fq::chunkFormat(fqchunk, data.vec);
      data.chunk_p = fqchunk->chunk;
    }
    return data;
  }
		//typedef FM_NEOREF FDTYPE_NCP;     //formated type no-copy
  void release_chunk(FDTYPE_NCP &data){
    dp_->Release(data.chunk_p);
  }

  ~FQ_SE(){
    delete fqFileReader;
    delete dp_;
    delete dq_;
  }
};
class FQ_PE{
  typedef rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> FqChunkQueue;  
  typedef rabbit::fq::FastqDataPool FqDataPool;
public:  
  std::thread* producer_;
  rabbit::fq::FastqFileReader* fqFileReader;
  FqDataPool* dp_;
  FqChunkQueue* dq_;
public:
  typedef std::pair<std::vector<Reference>, std::vector<Reference> > FDTYPE;
  typedef std::pair<FM_NEOREF, FM_NEOREF> FDTYPE_NCP;     //formated type no-copy
  FQ_PE(std::string file1, std::string file2){
    dq_ = new  FqChunkQueue(128, 1); // data queue
    dp_ = new FqDataPool(256, 1 << 24); // data pool
		bool zipped = false;
		if(ends_with(file1, ".gz")) zipped = true;
    fqFileReader = new rabbit::fq::FastqFileReader(file1, *dp_, zipped, file2);
    //start_producer();
  }
  void start_producer(){
    producer_ = new std::thread([&] {
      int n_chunks = 0;
      while (true) {
        rabbit::fq::FastqPairChunk *fqchunk = new rabbit::fq::FastqPairChunk;
        fqchunk->chunk = fqFileReader->readNextPairChunk1();
        if (fqchunk->chunk == NULL) break;
        n_chunks++;
        dq_->Push(n_chunks, fqchunk->chunk);
      }

      dq_->SetCompleted();
    }); 
  }
  FDTYPE get_formated_reads() {
    rabbit::fq::FastqPairChunk *fqchunk = new rabbit::fq::FastqPairChunk(); 
    rabbit::int64 id = 0;
    FDTYPE res_data;
    std::vector<Reference> &left_data = res_data.first;
    std::vector<Reference> &right_data = res_data.second;
    int ref_num_l, ref_num_r;
    if (dq_->Pop(id, fqchunk->chunk)) {
      ref_num_l = rabbit::fq::chunkFormat(fqchunk->chunk->left_part, left_data, true);
      ref_num_r = rabbit::fq::chunkFormat(fqchunk->chunk->right_part, right_data, true);
      assert(ref_num_l == ref_num_r); //TODO: add failed info message
      //------------------relaease-----------------//
      dp_->Release(fqchunk->chunk->left_part);
      dp_->Release(fqchunk->chunk->right_part);
    }else{
      ref_num_l = 0;
    }
		delete fqchunk;
    return res_data;
  }
  FDTYPE_NCP get_formated_reads_nocp() {
    rabbit::fq::FastqPairChunk *fqchunk = new rabbit::fq::FastqPairChunk; //TODO: it will become wild pointer!!!!
    rabbit::int64 id = 0;
    FDTYPE_NCP res_data;
    std::vector<neoReference> &left_data = res_data.first.vec;
    std::vector<neoReference> &right_data = res_data.second.vec;
    int ref_num_l, ref_num_r;
    if (dq_->Pop(id, fqchunk->chunk)) {
      ref_num_l = rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk*)fqchunk->chunk->left_part, left_data);
      ref_num_r = rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk*)fqchunk->chunk->right_part, right_data);
      assert(ref_num_l == ref_num_r); //TODO: add failed info message
      res_data.first.chunk_p = fqchunk->chunk->left_part;
      res_data.second.chunk_p = fqchunk->chunk->right_part;
    }else{
      res_data.first.chunk_p = NULL;
      res_data.second.chunk_p = NULL;
      ref_num_l = 0;
    }
    return res_data;
  }

  void release_chunk(FDTYPE_NCP &data){
    dp_->Release(data.first.chunk_p);
    dp_->Release(data.second.chunk_p);
  }

  ~FQ_PE(){
    delete fqFileReader;
    delete dp_;
    delete dq_;
    delete producer_;
  }
};

class FA{
typedef rabbit::core::TDataQueue<rabbit::fa::FastaChunk> FaChunkQueue;
typedef rabbit::fa::FastaDataPool FaDataPool;
public:
  typedef std::vector<Reference> FDTYPE;
  typedef std::vector<Reference> FDTYPE_NCP;
  rabbit::fa::FastaFileReader *faFileReader;
  FaDataPool* dp_;
  FaChunkQueue* dq_;
  std::thread* producer_;
  std::thread** consumers_;
  int consumer_th_;
public:
  FA(std::string &file){
    printf("in constuctor of FA!, fname: %s\n", file.c_str());
    dq_ = new  FaChunkQueue(128, 1); // data queue
    dp_ = new FaDataPool(256, 1 << 24); // data pool
		bool zipped = false;
		if(ends_with(file, ".gz")) zipped = true;
    faFileReader = new rabbit::fa::FastaFileReader(file, dp_, zipped);
  }

  void start_producer() {
    producer_ = new std::thread([&]{
      int n_chunks = 0;
      while (true) {
        rabbit::fa::FastaChunk *fachunk;// = new rabbit::fa::FastaChunk;
        fachunk = faFileReader->readNextChunkList();
        //fachunk = faFileReader->readNextChunk();
        if (fachunk == NULL) break;
        n_chunks++;
        dq_->Push(n_chunks, fachunk);
      }
      dq_->SetCompleted();
    });
  }
  template<typename T>
	void process_data_mt(const int tn, T (*func_ptr)(Reference &), std::vector<std::vector<T> > &v_res_data){ //tn means thread number
		consumers_ = new std::thread*[tn];
    assert(v_res_data.size() == tn);
    this->consumer_th_ = tn;
    for (int i = 0; i < tn; i++) {
      rabbit::int64 id = 0;
      std::vector<T> &res_data = v_res_data[i];

      consumers_[i] = new std::thread([&]() {
				//cout << "starting worker: " << i << endl;
        rabbit::fa::FastaChunk *fachunk;  // = new rabbit::fa::FastaChunk;
        while (dq_->Pop(id, fachunk)) {
          // rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
          std::vector<Reference> data;
          int ref_num = rabbit::fa::chunkFormat(*fachunk, data);
          for (Reference &r : data) {
            res_data.emplace_back(func_ptr(r));
          }
          //------------------relaease-----------------//
          rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
          do {
            dp_->Release(tmp);
            tmp = tmp->next;
          } while (tmp != NULL);
        }
      });
    }
  }

  std::vector<Reference> get_formated_reads(){
    rabbit::fa::FastaChunk *fachunk;  // = new rabbit::fa::FastaChunk;
    rabbit::int64 id = 0;
    int ref_num;
		std::vector<Reference> data;
    if (dq_->Pop(id, fachunk)) {
      // rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
      ref_num = rabbit::fa::chunkListFormat(*fachunk, data);
      //------------------relaease-----------------//
      release_chunk(fachunk);
    }else{
      ref_num = 0;
    }
    return data;
  }

  void release_chunk(rabbit::fa::FastaChunk *fachunk) {
    rabbit::fa::FastaDataChunk *tmp = fachunk->chunk;
    do {
      dp_->Release(tmp);
      tmp = tmp->next;
    } while (tmp != NULL);
  }

  ~FA(){
    delete faFileReader;
    delete dq_;
    delete dp_;
    delete producer_;
  }
};


// Just a wapper
template<class T>
class FXReader{
public:
  T reader_;
  FXReader(std::string fn)
    : reader_(fn)
  {
    reader_.start_producer();
  }
  FXReader(const std::string &fn1, const std::string &fn2)
    : reader_(fn1, fn2)
  {
    reader_.start_producer();
  }
 // bool pop_dq(rabbit::int64 &id, ){
 //   return reader_.dq_.Pop(id, chunk);
 // }
  typename T::FDTYPE get_formated_reads(){ // copy
    return reader_.get_formated_reads();
  }
  typename T::FDTYPE_NCP get_formated_reads_nocp(){ //no copy 
    return reader_.get_formated_reads_nocp();
  }
	void release_chunk(typename T::FDTYPE_NCP& data){
		reader_.release_chunk(data);
	}
  void join_producer(){
    reader_.producer_->join();
  }
  void join_consumers(){
    for(int i = 0; i < reader_.consumer_th_; i++){
      reader_.consumers_[i]->join();
    }
  }
};

#endif
