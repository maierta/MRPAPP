//-*-C++-*-

#ifndef PSIMAG_STLUtil_H_
#define PSIMAG_STLUtil_H_


#include<vector>
#include<map>
#include<iostream>



//============================================================================
//==========     Compare the key values of two maps                 ==========
//============================================================================

/** Returns true if maps have identical key spaces.
 *  \author Gregory Brown
 */
template<class K, class T, class C, class A>
bool SameKeys(const std::map<K,T,C,A>& obj1, const std::map<K,T,C,A>& obj2)
{
  if( obj1.size() != obj2.size() ) return false;
  typename std::map<K,T,C,A>::const_iterator itr1;
  typename std::map<K,T,C,A>::const_iterator itr2;
  itr1=obj1.begin();
  itr2=obj2.begin();
  while( itr1!=obj1.end() )
  {
    //if( obj1.key_compare()(itr1->first,itr2->first) ||
    //    obj1.key_compare()(itr2->first,itr1->first) ) return false;
    if( (itr1->first<itr2->first) || (itr2->first<itr1->first) ) 
      return false;
    ++itr1; ++itr2;
  }
  return true;
}

//============================================================================
//==========     Stream input and output of std::vector             ==========
//============================================================================

/** Stream output of maps.
 *  The format of the output is a line with the number of elements in 
 *  the map followed by an integer. If that integer is less than zero
 *  each key is read own its own, but if the integer is positive all the
 *  keys have the same length and some compaction has been done. This 
 *  compaction necessitates a specialized operator>> for reading.
 *
 *  \author Gregory Brown
 */
template<class T, class A>
std::ostream& operator<<(std::ostream& os, const std::vector<T,A>& obj)
{
  typedef std::vector<T,A> VectorType;
  os << obj.size() << "\t" << "-1" << std::endl;
  for(typename VectorType::const_iterator itr=obj.begin(); itr!=obj.end(); itr++) 
  {
    os << *itr << std::endl;
  }
  return os;
}


template<class T, class A>
std::istream& operator>>(std::istream& is, std::vector<T,A>& obj)
{
  typedef std::vector<T,A> VectorType;
  long int nel,dummy;
  is >> nel >> dummy;
  if( dummy>=0 )
  {
    std::ostringstream msg;
    msg << "operator>> error: Generic version called, but element length not negative";
    throw std::logic_error(msg.str());
  }
  obj.resize(nel);
  for(long int iel=0; iel<nel; iel++)
  {
    is >> obj[iel];
  }
  return is;
}

//============================================================================
//==========     Stream input and output of std::map                ==========
//============================================================================

/** Stream output of maps.
 *  The format of the output is a line with the number of elements in 
 *  the map followed by an integer. If that integer is less than zero
 *  each key is read own its own, but if the integer is positive all the
 *  keys have the same length and some compaction has been done. This 
 *  compaction necessitates a specialized operator>> for reading.
 *
 *  \author Gregory Brown
 */
template<class K, class T, class C, class A>
std::ostream& operator<<(std::ostream& os, const std::map<K,T,C,A>& obj)
{
  typedef std::map<K,T,C,A> MapType;
  os << obj.size() << "\t" << "-1" << std::endl;
  for(typename MapType::const_iterator itr=obj.begin(); itr!=obj.end(); itr++) 
  {
    os << itr->first << "\t" << itr->second << std::endl;
  }
  return os;
}


template<class TO, 
class CO, 
class AO, 
class KI, 
class TI, 
class XI, 
class AI>
std::ostream& operator<<(std::ostream& os,
		         const std::map<std::map<KI,TI,XI,AI>,TO,CO,AO>& obj) 
{
  typedef std::map<KI,TI,XI,AI>    Inner;
  typedef std::map<Inner,TO,CO,AO> Outer;
  // Write the length of the object
  long int nel = obj.size();
  os << nel << "\t";
  // Short circuit if map is empty;
  if( nel==0 )
  {
    os << "0" << std::endl;
    return os;
  }
  // Find if all maps have the same keys
  long int nkey = 0;
  bool keysame = true;
  {
    typename Outer::const_iterator itr1,itr2;
    itr1 = obj.begin();
    itr2 = itr1; ++itr1;
    while( keysame && itr1 != obj.end() )
    {
      nkey = itr1->first.size();
      keysame = SameKeys(itr1->first,itr2->first);
      itr2 = itr1; ++itr1;
    }
  }
  if( keysame )
  {
    os << nkey;
    const Inner f = obj.begin()->first;
    for(typename Inner::const_iterator itr=f.begin(); itr!=f.end(); ++itr)
      os << "\t" << itr->first;
    os << std::endl;
  }
  else
    os << "-1" << std::endl;
  // Write out the contents
  for(typename Outer::const_iterator itrO=obj.begin(); itrO!=obj.end(); ++itrO)
  {
    // Output of the map itr0->first
    if( !keysame ) os << itrO->first.size() << "\t";
    for(typename Inner::const_iterator itrI=itrO->first.begin(); 
        itrI!=itrO->first.end(); ++itrI)
    {
      if( keysame )
        os << itrI->second << "\t";
      else
        os << itrI->first << " " << itrI->second << "\t";
    }
    os << itrO->second << std::endl;
  }
  return os;
}       // std::ostream operator<< map-keyed-with-map



template<class K, class T, class C, class A>
std::istream& operator>>(std::istream& is, std::map<K,T,C,A>& obj)
{
  typedef std::map<K,T,C,A> MapType;
  long int nel,dummy;
  is >> nel >> dummy;
  if( dummy>=0 )
  {
    std::ostringstream msg;
    msg << "operator>> error: Generic version called, by keysize not negative";
    throw std::logic_error(msg.str());
  }
  for(long int iel=0; iel<nel; iel++)
  {
    K key;
    T value;
    is >> key >> value;
    obj[ key ] = value;
  }
  return is;
}



template<class TO, class CO, class AO, class KI, class TI, class XI, class AI>
std::istream& operator>>(std::istream& is,
		         std::map<std::map<KI,TI,XI,AI>,TO,CO,AO>& obj) 
{
  typedef std::map<KI,TI,XI,AI>    Inner;
  typedef std::map<Inner,TO,CO,AO> Outer;
  // Key object for creating insertions
  Inner inner;
  // Reset state of the object
  obj.clear();
  // Get the size of the outer object
  long int nel;
  is >> nel; 
  // Short circuit if map is empty;
  if( nel==0 )
  {
    long int dummy;
    is >> dummy;
    return is;
  }
  // Find out if all maps have the same keys
  long int nkey;
  is >> nkey;
  bool keysame = nkey>=0;
  if( keysame )
  {
    // get the keys
    KI rkey;
    TI dummy;
    inner.clear();
    for(long int ikey=0; ikey<nkey; ikey++)
    {
      is >> rkey;
      inner[ rkey ] = dummy;
    }
  }
  // Read in the contents
  for(long int iel=0; iel<nel; iel++)
  {
    if( !keysame ) 
    {
      // each map is different
      is >> nkey;
      KI rkey;
      TI value;
      inner.clear();
      for(long int ikey=0; ikey<nkey; ikey++)
      {
        is >> rkey >> value;
        inner[ rkey ] = value;
      }
    }
    else
    {
      typename Inner::iterator kitr = inner.begin();
      for(long int ikey=0; ikey<nkey; ikey++)
      {
        is >> kitr->second;
	kitr++;
      }
    }
    TO value;
    is >> value;
    obj.insert( make_pair(inner,value) );
  }
  return is;
}       // std::istream operator>> map-keyed-with-map


//============================================================================
//==========     static_casting of objects to std::vector           ==========
//============================================================================

namespace psimag {


/** Copy various types into a vector
 *  This is primarily used in generic algorithms to convert types of unknown, 
 *  but fixed, dimension to a container whose size is determined at run time.
 *  The default behavior for scalars, defined here, creates a std::vector 
 *  with one element. 
 *
 *  \author Gregory Brown
 */
template<class T, class U>
void Copy(T i, std::vector<U>& result)
{ result.resize(1,i); }


/** Copy, potentially between vectors of different value_type.
 *  This is provided mostly so vectors can be created in generic
 *  situations, where the source could easily be a vector.
 *
 *  \author Gregory Brown
 */
template<class T, class U>
void Copy(const std::vector<T>& source, std::vector<U>& result)
{ result.assign(source.begin(),source.end()); }


/** Copy mapped_type of a map into a vector.
 *  The elements will be in the order defined by the map.
 *
 *  \author Gregory Brown
 */
template<class K, class M, class U>
void Copy(const std::map<K,M>& source, std::vector<U>& result)
{ 
  result.resize(source.size());
  typename std::vector<U>::iterator ritr=result.begin();
  for(typename std::map<K,M>::const_iterator itr=source.begin(); itr!=source.end(); itr++)
    *(ritr++) = itr->second;
}


}       // namespace psimag


#endif  // PSIMAG_STLUtil_H_
