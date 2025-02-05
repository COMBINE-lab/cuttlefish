/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include <cstdio>
#include <vector>
#include <map>
#include <cstdio>

#include "Reference.h"
#include "Formater.h"
#include "Buffer.h"
#include "FastxStream.h" 
//#include "read.h"
//#include "Sequence.h"

using namespace std;

namespace rabbit
{

namespace fa
{
/*
FastaChunk* FastaReader::readNextChunk(){
	FastaDataChunk* part = NULL;
	recordsPool.Acquire(part);
	FastaChunk *dataPart = new FastaChunk;
	dataPart->chunk = part;
	if(fileReader.ReadNextChunk(dataPart, this->seqInfos))
	{
		return dataPart;
	}
	else
	{
		recordsPool.Release(part);
		return NULL;
	}
}
*/

//int chunkFormat(FastaDataChunk* &chunk, std::vector<Sequence*> &data, bool mHasQuality){
//	//format a whole chunk and return number of reads
//	int seq_count = 0;
//	int line_count = 0;
//	int pos_ = 0;
//
//	while(true){
//		//TODO rewrite to deal with part sequence
//		//string name = getLine(chunk, pos_);
//		//if(name.empty()) break;//dsrc guarantees that read are completed!
//		//std::cout << name << std::endl;
//
//		//string sequence = getSequence(chunk, pos_);
//		//std::cout << sequence << std::endl;
//
//		//data.push_back(new Read(name, sequence));
//		//seq_count++;
//
//	}
//
//	return seq_count;
//}


//FIXME:support enter = \n only
string getSequence(FastaDataChunk * &chunk, uint64 &pos)
{
	int start_pos = pos;
	char * data = (char *)chunk->data.Pointer();
	//cerr << "start pos: " << pos << endl << flush;
	//cerr << "chunk size: " << chunk->size << endl << flush;
	//cerr << "data[pos]: " << (int)data[pos] << endl << flush;
	string res="";

	while(pos < chunk->size - 1 )
	{
		if(data[pos] == '\n'){
			res += string(data+start_pos, pos-start_pos);
			pos++;
			start_pos = pos;
			if(data[pos] == '>')
				return res;		
		}
		else{
			pos++;
		}
	}

	//deal with last char	
	if(pos == chunk->size - 1)
	{
		if(data[pos] == '\n')
			res += string(data+start_pos, pos - start_pos);
		else
			res += string(data+start_pos, pos - start_pos + 1);

		return res;
	}	

		
	return "";
}

//only support uinx-like '\n'
string getLine(FastaDataChunk * &chunk, uint64 &pos)
{
	int start_pos = pos;
	char* data = (char *)chunk->data.Pointer();

	while(pos < chunk->size){
		if(data[pos] == '\n'){
			pos++;
			return string(data+start_pos, pos-start_pos - 1);
		}
		else{
			pos++;
		}
	}

	return "";
}

int chunkFormat(FastaChunk & fachunk, vector<Reference> & refs)
{
	uint64 pos = 0;
	bool done = false;
	//cerr << "into chunkFormat" << endl;
	while(true){
	
		Reference ref = getNextSeq(fachunk, done, pos);
		if(done) break;
		refs.push_back(ref);
	}

	ASSERT(refs.size() == fachunk.nseqs);

	return refs.size();
	
}

//design to filter out kmer apps
int chunkFormat(FastaChunk & fachunk, vector<Reference> & refs, int kmerSize)
{
	uint64 pos = 0;
	bool done = false;
	//cerr << "into chunkFormat" << endl;
	uint64 short_count = 0; //counter for seqs shorter than kmer
	while(true){
	
		Reference ref = getNextSeq(fachunk, done, pos);
		if(done) break;
		if(ref.seq.length() < static_cast<std::size_t>(kmerSize))
			short_count++;
		else
			refs.push_back(ref);
	}

	ASSERT(refs.size() + short_count == fachunk.nseqs);

	return refs.size();
	
}

Reference getNextSeq(FastaChunk & fachunk, bool & done, uint64 & pos)
{
	Reference ref;
	if(pos >= fachunk.chunk->size - 1){
		done = true;
		return ref;
	}

	char *data = (char *)fachunk.chunk->data.Pointer();	

	//while(data[pos] == '\n') pos++;//jump empty lines

	if(data[pos] != '>')
	{
		ref.seq = getSequence(fachunk.chunk, pos);	
		ref.length = ref.seq.size();
		ref.gid = fachunk.start;
		fachunk.start++;
	}else{
		string line = getLine(fachunk.chunk, pos);
		int str_pos = line.find_first_of(' ');
		ref.name = line.substr(1, str_pos-1);	//remove '>' and ' '
		if(static_cast<std::size_t>(str_pos) < line.size()) 
			ref.comment = line.substr(str_pos + 1);
		//cerr << "name: " << ref.name << endl << flush;
		ref.seq = getSequence(fachunk.chunk, pos);	
		//cerr << "seq: " << ref.seq << endl << flush;
		ref.length = ref.seq.size();
		ref.gid = fachunk.start;
		fachunk.start++;
	}

	return ref;
}

} // namesapce fa

namespace fq
{

/*
void FastqReader::readChunk()
{
	FastqDataChunk* part = NULL;

	recordsPool.Acquire(part);
	printf("fastqio: ready to into while\n");

	//while (!errorHandler.IsError() && fileReader.ReadNextChunk(part))
	while(fileReader.ReadNextChunk(part))
	{
		ASSERT(part->size > 0);

		//recordsQueue.Push(numParts, part); //[haoz:] Push:把<numparts, part>组成pair然后放到recordsQueue里面然后notifiy_one
		printf("numParts is %d\n", numParts);
		numParts++;

		recordsPool.Release(part);	
		recordsPool.Acquire(part);
	}

	ASSERT(part->size == 0);
	recordsPool.Release(part);		// the last empty part

	//recordsQueue.SetCompleted();
}

FastqDataChunk* FastqReader::readNextChunk(){
	FastqDataChunk* part = NULL;
	recordsPool.Acquire(part);
	if(fileReader.ReadNextChunk(part))
	{
		return part;
	}
	else
	{
		recordsPool.Release(part);
		return NULL;
	}
}
*/
//FastqDataChunk* FastqReader::readNextPairedChunk(){
//	FastqDataChunk* part = NULL;
//	recordsPool.Acquire(part);
//	if(fileReader.ReadNextPairedChunk(part))
//	{
//		//std::cout <<"read one chunk:   " << std::endl;
//		//std::cerr << (char*)part->data.Pointer() << endl;
//		//std::cerr << "chunk end" << endl;
//		return part;
//	}
//	else
//	{
//		//std::cout <<"read one chunk failed :   " << std::endl;
//		recordsPool.Release(part);
//		return NULL;
//	}
//}
void print_read(neoReference& ref){
	std::cout << std::string((char*)ref.base+ref.pname, ref.lname) << std::endl;
	std::cout << std::string((char*)ref.base+ref.pseq, ref.lseq) << std::endl;
	std::cout << std::string((char*)ref.base+ref.pstrand, ref.lstrand) << std::endl;
	std::cout << std::string((char*)ref.base+ref.pqual, ref.lqual) << std::endl;
}

int chunkFormat(FastqChunk* &fqChunk, std::vector<neoReference> &data, bool mHasQuality = true){
	(void)mHasQuality;
	FastqDataChunk * chunk = fqChunk->chunk;
	uint64_t seq_count = 0;
	// uint64_t line_count = 0;
	uint64_t pos_ = 0;
	neoReference ref;
	while(true){
		ref.base = chunk->data.Pointer();
		ref.pname = pos_;
		if(neoGetLine(chunk, pos_, ref.lname)){
			ref.pseq = pos_; 
		} 
		else{ break;}
		neoGetLine(chunk, pos_, ref.lseq); 
		ref.pstrand = pos_; 
		neoGetLine(chunk, pos_, ref.lstrand); 
		ref.pqual = pos_;  
		neoGetLine(chunk, pos_, ref.lqual);
		seq_count++;
		//print_read(ref);
		data.emplace_back(ref);
	}

	return seq_count;
}
/*
int chunkFormat(FastqChunk* &fqChunk, std::vector<Reference> &data, bool mHasQuality = true){
	//format a whole chunk and return number of reads
	FastqDataChunk * chunk = fqChunk->chunk;
	int seq_count = 0;
	int line_count = 0;
	int pos_ = 0;
	Reference ref;

	while(true){
		string name = getLine(chunk, pos_);
		if(name.empty()) break;//dsrc guarantees that read are completed!
		//std::cerr << name << std::endl;

		string sequence = getLine(chunk, pos_);
		//std::cerr<< sequence << std::endl;

		string strand = getLine(chunk, pos_);
		//std::cerr << strand << std::endl;
		ref.name = name;
		ref.seq = sequence;
		ref.strand = strand;

		if(!mHasQuality){
			string quality = string(sequence.length(), 'K');
			//std::cerr << quality << std::endl;
			//data.push_back(new Read(name, sequence, strand, quality));
			//data.emplace_back(Sketch::Reference(name, "", sequence, strand, quality));
			ref.quality = quality;
			data.push_back(ref);
			seq_count++;

		}else{
			string quality = getLine(chunk, pos_);
			//std::cerr << quality << std::endl;
			//data.push_back(new Read(name, sequence, strand, quality));
			//data.emplace_back(Sketch::Reference(name, "", sequence, strand, quality));
			ref.quality = quality;
			data.push_back(ref);
			seq_count++;
		}
	}

	return seq_count;
}
*/
//int pairedChunkFormat(FastqDataChunk* &chunk, std::vector<ReadPair*> &data, bool mHasQuality){
//	//format a whole chunk and return number of reads
//	int seq_count = 0;
//	int pos_ = 0;
//	//处理一开始的两条序列
//	Read* first = getOnePairedRead(chunk,pos_,mHasQuality);
//	Read* second = getOnePairedRead(chunk,pos_,mHasQuality);
//
//
//	if(first->mName.substr(0, first->mName.find(' ')) == second->mName.substr(0, second->mName.find(' ')) ){
//		//cerr << first->mName<<" + " << endl;
//
//		data.push_back(new ReadPair(first,second));
//		seq_count++;
//
//	}else{
//		//std::cout<<"Right"<<std::endl;
//		//cerr << first->mName<<" - "<<endl;
//		Read* third = getOnePairedRead(chunk,pos_,mHasQuality);
//		data.push_back(new ReadPair(second,third));
//		seq_count++;
//	}
//	while(true){
//		Read* t1=getOnePairedRead(chunk,pos_,mHasQuality);
//		if(t1 == NULL){
//			break;
//		}
//		Read* t2 =getOnePairedRead(chunk,pos_,mHasQuality);
//		if(t2 == NULL){
//			//cerr << t1->mName << " with + without -" << endl;
//			break;
//		}
//		data.push_back(new ReadPair(t1,t2));
//		seq_count++;
//	}
//
//	//for(int i =0;i<seq_count;i++){
//	//	cerr << data[i]->mLeft->mName<<endl;
//	//	cerr << data[i]->mLeft->mSeq.mStr<<endl;
//	//	cerr << data[i]->mLeft->mStrand<<endl;
//	//	cerr << data[i]->mLeft->mQuality<<endl;
//	//	cerr << data[i]->mRight->mName<<endl;
//	//	cerr << data[i]->mRight->mSeq.mStr<<endl;
//	//	cerr << data[i]->mRight->mStrand<<endl;
//	//	cerr << data[i]->mRight->mQuality<<endl;
//	//	
//	//}
//	//cerr << data[seq_count-1]->mLeft->mName<<endl;
//		
//	return seq_count;
//}
//Read* getOnePairedRead(FastqDataChunk* &chunk,int &pos_, bool mHasQuality){
//	while(true){
//		string name = getLine(chunk, pos_);
//		if(name.empty()) return NULL;
//		//std::cerr << name << std::endl;
//		string sequence = getLine(chunk, pos_);
//
//		string strand = getLine(chunk, pos_);
//		if(!mHasQuality){
//			string quality = string(sequence.length(), 'K');
//			return new Read(name, sequence, strand, quality);
//	
//		}else{
//			string quality = getLine(chunk, pos_);
//			//std::cerr << quality << std::endl;
//			return new Read(name, sequence, strand, quality);
//		}
//	}
//}


string getLine(FastqDataChunk* &chunk, int &pos){
	int start_pos = pos;
	char* data = (char *)chunk->data.Pointer();

	const int end = chunk->size + 1;
	while(pos <= end){
		if(data[pos] == '\n' || data[pos] == '\r' || pos == end){
			//find a line
			pos++;
			return string(data+start_pos, pos-start_pos - 1);
		}
		else{
			pos++;
		}
	}
	return "";
}
int neoGetLine(FastqDataChunk* &chunk, uint64_t &pos, uint64_t &len){
	int start_pos = pos;
	const char* data = (char *)chunk->data.Pointer();
	const uint64_t chunk_size = chunk->size + 1;
	/*
	while (pos <= chunk_size && data[pos] != '\n'){
		pos++;
	}
	len = pos - start_pos - 1;
	return len;
	*/
	//while(pos <= (chunk->size + 1)){
	while(pos <= chunk_size){
		//if(data[pos] == '\n' || data[pos] == '\r' || pos == (chunk->size + 1)){
		if(data[pos] == '\n'){
			//find a line
			pos++;
			//return string(data+start_pos, pos-start_pos - 1);
			//pos = start_pos; 
			len = pos - start_pos - 1;
			return len;
		}
		else{
			pos++;
		}
	}
	return 0;
}
/*
int neoGetLine(FastqDataChunk* &chunk, char* &pos, uint64_t &len){
	char* start_pos  = pos;
	char* data = (char *)chunk->data.Pointer();
	char* chunkend = data + chunk->size + 1; //the end address of this chunk(memory address)
	while(pos <= chunkend){
		if(*pos == '\n' || *pos == '\r' || pos == chunkend){
			//find a line
			pos++;
			len = (pos - start_pos)/sizeof(char); //-------if need the division??-------//
			//return string(data+start_pos, pos-start_pos - 1);
			return len;
		}
		else{
			pos++;
		}
	}
	return 0;
}
*/
} // namesapce fq

} // namespace rabbit
