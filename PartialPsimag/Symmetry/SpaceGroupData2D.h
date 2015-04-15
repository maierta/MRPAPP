//-*-C++-*-

#ifndef  Psimag_SpaceGroupData2D
#define  Psimag_SpaceGroupData2D

/** \ingroup symmetryConcepts */
/*@{*/

/** \file SpaceGroup2D.h
 *
 *  1999, Grosse-Kunstleve, "Algorithms for deriving crystallographic space-group information"
 *
 *  some of this will be implemented in the cctbx where I may be able to reuse it.
 */  

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <bitset>

#include "SymmetryOperation.h"
#include "SymmetryElement.h"
#include "SymmetryElements2D.h"
#include "CellPosition.h"

namespace psimag {
  

  /** \ingroup symmetryConcepts 
   *
   * \brief The SpaceGroup class
   */
  template<typename Field, size_t DIM, size_t NUM, size_t ALTERNATE, typename Algorithms> 
  class SpaceGroupData
  {};


  template<typename Field, size_t DIM, typename Algorithms> class SymmetryElements;

  //====================================================================== 1
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 1 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,1,1,Algorithms>
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=1, NumGenerators=1, NumOperations=1};

    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("1 Oblique p1");
    }

    static const std::string origin() {
      return std::string("arbitrary");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","1"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
    }
    
  };

  //====================================================================== 2
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 2 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,2,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=2, NumGenerators=2, NumOperations=5};

    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("2 Oblique p2");
    }

    static const std::string origin() {
      return std::string("at twoFold");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","0"));
      return result;
    }

     static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));

      
    }
  };

  //====================================================================== 3
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 3 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,3,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=3, NumGenerators=2, NumOperations=3};


    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 


    static const std::string name() {
      return std::string("3 Rectangular Pm");
    }

    static const std::string origin() {
      return std::string("on b mirror line");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror(1/2,0)"));
      
    }

  };

  //====================================================================== 4
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 4 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,4,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=4, NumGenerators=2, NumOperations=3};


    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("4 Rectangular pg");
    }

    static const std::string origin() {
      return std::string("on b glide line");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(1/2,0)"));
    }

  };
  //====================================================================== 5
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 5 class
   *
   * This will need more work later! &*&*&* depends on the reduced lattice closely.
   *
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,5,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=true, SpaceGroupNumber=5, NumGenerators=3, NumOperations=5};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("5 Rectangular cm");
    }

    static const std::string origin() {
      return std::string("on b mirror line");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(1/4,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(3/4,0)"));
    }

  };
  //====================================================================== 6
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData g class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,6,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=6, NumGenerators=3, NumOperations=8 };

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("6 Rectangular p2mm");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror(1/2,0)"));
    }

  };

  //====================================================================== 7
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 7 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,7,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=7 , NumGenerators=3, NumOperations=4};

    typedef SymmetryOperation<Field,DIM,Algorithms>  SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("7 Rectangular p2mg");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/4"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","1/4"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror(1/4,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror(3/4,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide(0,1/2)"));
    }

  };

  //====================================================================== 8
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 8 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,8,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=8, NumGenerators=3, NumOperations=4};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType;
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("8 Rectangular p2gg");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide(0,1/4)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide(0,3/4)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(1/4,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(3/4,0)"));
    }

  };
  //====================================================================== 9
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 9 class
   *
   * Will require some more work! &*&*&*
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,9,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=true, SpaceGroupNumber=9 , NumGenerators=4, NumOperations=11};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType;

    static const std::string name() {
      return std::string("9 Rectangular c2mm");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/4"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/4"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("2a -b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 2a -b glide"));
      //      ids.push_back(SymmetryElementsType::getIndexFor("2a -b glide(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(1/2,0)"));
      
    }

  };

  //====================================================================== 10
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 10 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,10,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=10, NumGenerators=3, NumOperations=4};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("10  Square p4");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("fourFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFold(1/2,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFoldN(1/2,1/2)"));

      
    }

  };

  //====================================================================== 11
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 11 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,11,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=11, NumGenerators=4, NumOperations=8};

    typedef SymmetryOperation<Field,DIM,Algorithms>  SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms>   SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("11 Square p4mm");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));

      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("fourFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFold(1/2,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFoldN(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("a mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("a mirror(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("a b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("a -b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a b glide(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a -b glide(1/2,0)"));
    }

  };
  //====================================================================== 12
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 12 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,12,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=12, NumGenerators=4, NumOperations=8};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("12 Square p4gm");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));

      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("fourFold"));     
      ids.push_back(SymmetryElementsType::getIndexFor("fourFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFold(1/2,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("fourFoldN(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide(0,1/4)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide(0,3/4)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(1/4,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(3/4,0)"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a b glide"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a -b glide"));

      ids.push_back(SymmetryElementsType::getIndexFor("a b mirror(1/2,0)"));

      ids.push_back(SymmetryElementsType::getIndexFor("a -b mirror(1/2,0)"));
      
    }

  };

  //====================================================================== 13
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 13 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,13,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=13 , NumGenerators=2, NumOperations=3};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType;
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("13 Hexagonal p3");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("2/3","1/3"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/3","2/3"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","1/2"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));

      ids.push_back(SymmetryElementsType::getIndexFor("threeFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(2/3,2/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(2/3,2/3)"));

      
    }

  };
  //====================================================================== 14
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 14 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,14,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=14 , NumGenerators=3, NumOperations=6};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("14 Hexagonal p3m1");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("2/3","1/3"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/3","2/3"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));

      ids.push_back(SymmetryElementsType::getIndexFor("threeFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(2/3,2/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(2/3,2/3)"));

      ids.push_back(SymmetryElementsType::getIndexFor("a b mirror"));

      ids.push_back(SymmetryElementsType::getIndexFor("2a -b mirror"));

      ids.push_back(SymmetryElementsType::getIndexFor("a -2b mirror(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("a -2b mirror"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a -2b glide(1/4,0)"));
      //      ids.push_back(SymmetryElementsType::getIndexFor("-a 2b glide(3/4,0)"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 2a -b glide(1/2,0)"));
      
    }

  };


  //====================================================================== 15
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 15 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,15,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=15, NumGenerators=3, NumOperations=6};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("15 Hexagonal p31m");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("2/3","1/3"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/2"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));

      ids.push_back(SymmetryElementsType::getIndexFor("threeFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(2/3,2/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(2/3,2/3)"));

      ids.push_back(SymmetryElementsType::getIndexFor("a mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("a -b mirror"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(1/2,0)"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a -b glide(1/2,0)"));
      
    }

  };

  //====================================================================== 16
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 16 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,16,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=16, NumGenerators=3, NumOperations=6};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements<Field,DIM,Algorithms> SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("16 Hexagonal p6");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("2/3","1/3"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","1/2"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));

      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("threeFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(2/3,2/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(2/3,2/3)"));

      ids.push_back(SymmetryElementsType::getIndexFor("sixFold"));

      
    }

  };

  //====================================================================== 17
  //
  /** \ingroup symmetryConcepts 
   *
   * \brief The 2D SpaceGroupData 17 class
   */
  template<typename Field, typename Algorithms>
  class SpaceGroupData<Field,2,17,1,Algorithms>   
  {
  public: 

    enum { DIM=2, Centered=false, SpaceGroupNumber=17, NumGenerators=4, NumOperations=12};

    typedef SymmetryOperation<Field,DIM,Algorithms>    SymmetryOperationType;
    typedef SymmetryElements <Field,DIM,Algorithms>    SymmetryElementsType; 
    typedef CellPosition<Field,DIM,Algorithms>       CellPositionType; 

    static const std::string name() {
      return std::string("17 Hexagonal p6mm");
    }

    static const std::vector<CellPositionType> verticies() {
      std::vector<CellPositionType> result;
      result.push_back(cellPosition<Field,std::string,Algorithms>("0","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("1/2","0"));
      result.push_back(cellPosition<Field,std::string,Algorithms>("2/3","1/3"));
      return result;
    }

    static void operations(std::vector<int>& ids) 
    {
      
      ids.push_back(SymmetryElementsType::getIndexFor("identity"));

      ids.push_back(SymmetryElementsType::getIndexFor("twoFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("twoFold(1/2,1/2)"));

      ids.push_back(SymmetryElementsType::getIndexFor("threeFold"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(1/3,1/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFold(2/3,2/3)"));
      ids.push_back(SymmetryElementsType::getIndexFor("threeFoldN(2/3,2/3)"));

      ids.push_back(SymmetryElementsType::getIndexFor("sixFold"));

      ids.push_back(SymmetryElementsType::getIndexFor("a mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a glide(0,1/2)"));
      ids.push_back(SymmetryElementsType::getIndexFor("1/2 b glide(1/2,0)"));

      ids.push_back(SymmetryElementsType::getIndexFor("b mirror"));
      ids.push_back(SymmetryElementsType::getIndexFor("a -b mirror"));

      ids.push_back(SymmetryElementsType::getIndexFor("a -2b mirror(1/2,0)"));
      ids.push_back(SymmetryElementsType::getIndexFor("a -2b mirror"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 a -2b glide(1/4,0)"));
      //      ids.push_back(SymmetryElementsType::getIndexFor("-a 2b glide(3/4,0)"));
      //ids.push_back(SymmetryElementsType::getIndexFor("-a 2b glide(5/4,0)"));

      ids.push_back(SymmetryElementsType::getIndexFor("2a -b mirror"));

      ids.push_back(SymmetryElementsType::getIndexFor("1/2 2a -b glide(1/2,0)"));

    }

  };

} /** namespace spimag **/

#endif // Psimag_SpaceGroupData2D

/*@}*/
