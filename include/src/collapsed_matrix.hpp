#ifndef COLLAPSED_MATRIX_HPP
#define COLLAPSED_MATRIX_HPP

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>     

//modified from http://upcoder.com/2/efficient-vectors-of-vectors

template <typename T>
class CollapsedMatrix
{
  public:
    typedef typename std::vector<T>::size_type size_type;
    
    typedef typename std::vector<T>::iterator value_iterator;
    typedef typename std::vector<T>::const_iterator const_value_iterator;
    
    typedef typename std::vector<size_type>::iterator index_iterator;
    typedef typename std::vector<size_type>::const_iterator const_index_iterator;
    typedef typename std::vector<size_type>::reverse_iterator reverse_index_iterator;
    typedef typename std::vector<size_type>::const_reverse_iterator const_reverse_index_iterator;

  private:
    std::vector<T> _v;
    std::vector<size_type> _index;

    std::vector<T> sorted_v, unique_v;

  public:

    CollapsedMatrix() : _index(1){
      _index[0] = 0;
    }
    CollapsedMatrix(const std::vector<std::vector<T> >& buildFrom) : _index(buildFrom.size() + 1){
      _index.clear();
      _index.push_back(0);
      for(size_type i = 0; i != buildFrom.size(); ++i){
        _index.push_back(_index.back() + buildFrom[i].size());
      }
      _v.resize(_index.back());
      size_type vIndex = 0;
      for(size_type i = 0; i != buildFrom.size(); ++i){
        for(size_type j = 0; j != subVectorSize(i); ++j){
          _v[vIndex++] = buildFrom[i][j];
        }
      }
      assertD(vIndex == SizeL(_v));
    }

    void shrinkToFit(){
      if(_v.capacity() != _v.size()){
        std::vector<T>(_v).swap(_v);
      }
      if(_index.capacity() != _index.size()){
        std::vector<size_type>(_index).swap(_index);
      }
    }

    bool operator==(const CollapsedMatrix& rhs) const{
      return _index == rhs._index && _v == rhs._v;
    }

    void swap(CollapsedMatrix& rhs){
      _v.swap(rhs._v);
      _index.swap(rhs._index);
    }

    bool empty() const{
      return _index.size() == 1;
    }
    void clear(){
      _index.resize(1);
      _v.clear();
    }

    ///////Don't use these when dimensions are known, push_back is inefficient and not needed. Use resize instead
    void pushBackSubVector(){
      _index.push_back(_index.back());
    }
    void pushBackSubVector(size_type size){
      _index.push_back(_index.back() + size);
      _v.resize(_v.size() + size);
    }
    void pushBackInLastSubVector(const T& value){
      assert(!empty());
      _index.back()++;
      _v.push_back(value);
    }
    void popBackInLastSubVector(){
      assert(!empty());
      assert(subVectorSize(size() - 1) > 0);
      _index.back()--;
      _v.pop_back();
    }

    //I need member for resizing _index and _v given size as input.
    void resize(const long unsigned int row_size, const long unsigned int col_size){
      _v.resize(row_size * col_size);
      _index.resize(row_size + 1);
      for(size_type i = 0; i < row_size+1; i++){
        _index[i] = i*col_size;
      }
    }

    value_iterator vbegin(){
      return _v.begin();
    }
    value_iterator vend(){
      return _v.end();
    }
    
    index_iterator ibegin(){
      return _index.begin();
    }
    index_iterator iend(){
      return std::prev(_index.end(),1);
    }
    reverse_index_iterator irbegin(){
      return std::next(_index.rbegin(),1);
    }
    reverse_index_iterator irend(){
      return _index.rend();
    }

    value_iterator rowbegin(const index_iterator it){ 
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*it);
    }
    value_iterator rowend(const index_iterator it){
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*(it+1));
    }
    value_iterator rowbegin(const reverse_index_iterator it){ 
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*it);
    }
    value_iterator rowend(const reverse_index_iterator it){
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*(it-1));
    }
    value_iterator rowbegin(const unsigned int i){ 
      assert(i < _index.size() - 1);
      return std::next(_v.begin(),_index[i]);
    }
    value_iterator rowend(const unsigned int i){
      assert(i < _index.size() - 1);
      return std::next(_v.begin(),_index[i+1]);
    }
    const_value_iterator rowbegin(const index_iterator it) const{ 
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*it);
    }
    const_value_iterator rowend(const index_iterator it) const{
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*(it+1));
    }
    const_value_iterator rowbegin(const reverse_index_iterator it) const{ 
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*it);
    }
    const_value_iterator rowend(const reverse_index_iterator it) const{
      assert(*it < _index[_index.size() - 1]);
      return std::next(_v.begin(),*(it-1));
    }
    const_value_iterator rowbegin(const unsigned int i) const{ 
      assert(i < _index.size() - 1);
      return std::next(_v.begin(),_index[i]);
    }
    const_value_iterator rowend(const unsigned int i) const{
      assert(i < _index.size() - 1);
      return std::next(_v.begin(),_index[i+1]);
    }


    ///////////////////////

    size_type size() const{
      return _index.size() - 1;
    }
    const T* operator[](size_type i) const{
      return &_v[_index[i]];
    }
    T* operator[](size_type i){
      return &_v[_index[i]];
    }

    const T* back() const{
      assert(!empty());
      return &_v[_index[_index.size() - 2]];
    }
    T* back(){
      assert(!empty());
      return &_v[_index[_index.size() - 2]];
    }

    size_type subVectorSize(size_type i) const{
      return _index[i + 1] - _index[i];
    }

    //for general purpose
    void DumpToFile(FILE* pFile){ 
   
      assert(pFile != NULL);
      size_type isize          = this -> size();
      size_type isubVectorSize = this -> subVectorSize(0);
      fwrite(&isize, sizeof(size_type), 1, pFile);
      fwrite(&isubVectorSize, sizeof(size_type), 1, pFile);
      fwrite(&_v[0], sizeof(T), _v.size(), pFile);
    
    }

    void ReadFromFile(FILE* pFile){ 

      assert(pFile != NULL);
      size_type isize, isubVectorSize;
      //read into _v
      fread(&isize, sizeof(size_type), 1, pFile);
      fread(&isubVectorSize, sizeof(size_type), 1, pFile);
      resize(isize, isubVectorSize);
      fread(&_v[0], sizeof(T), isize * isubVectorSize, pFile);
     
    }

    //For stepping stones
    void DumpToFile(FILE* fp, int i, const std::vector<int>& boundarySNP, const std::vector<T>& logscales){ 

      //question is whether I can compress _v[_index[i]]
      
      size_type isubVectorSize = this -> subVectorSize(0);
      size_type isize = 1; 

      
      std::vector<int> times(isubVectorSize, 1);

      unique_v.resize(isubVectorSize);
      T current_val = _v[_index[i]];
      int k = 0;
      unique_v[k] = current_val;
      for(int j = 1; j < isubVectorSize; j++){
        if( std::fabs( current_val - _v[_index[i]+j] ) < 1e-3 * std::min(current_val,_v[_index[i]+j]) ){
          times[k]++;
        }else{
          current_val = _v[_index[i]+j];
          k++;
          unique_v[k] = current_val;
        }
      }
      k++;

      fwrite(&isize, sizeof(size_type), 1, fp);
      fwrite(&isubVectorSize, sizeof(size_type), 1, fp);

      fwrite(&boundarySNP[i], sizeof(int), 1, fp); 
      fwrite(&logscales[i], sizeof(T), 1, fp); 

      fwrite(&k, sizeof(int), 1, fp);
      fwrite(&unique_v[0], sizeof(T), k, fp);
      fwrite(&times[0], sizeof(int), k, fp);
     
      //fwrite(&_v[_index[i]], sizeof(T), isubVectorSize, fp);

    }

    //for stepping stone
    void ReadFromFile(FILE* pFile, int& boundarySNP, T& logscale){ 

      assert(pFile != NULL);
      size_type isize, isubVectorSize;
      //read into _v
      fread(&isize, sizeof(size_type), 1, pFile);
      fread(&isubVectorSize, sizeof(size_type), 1, pFile);
      fread(&boundarySNP, sizeof(int), isize, pFile);
      fread(&logscale, sizeof(T), isize, pFile);
      resize(isize, isubVectorSize);

      int k;
      fread(&k, sizeof(int), 1, pFile);
      std::vector<int> index(k, 0), times(k, 1);
      unique_v.resize(k);

      fread(&unique_v[0], sizeof(T), k, pFile);
      fread(&times[0], sizeof(int), k, pFile);

      int i = 0;
      for(int j = 0; j < k; j++){
        for(int tmp = 0; tmp < times[j]; tmp++){
          _v[i] = unique_v[j];
          i++;
        }
      }
      //fread(&_v[0], sizeof(T), isize * isubVectorSize, pFile);

    }




};

#endif //COLLAPSED_MATRIX_HPP
