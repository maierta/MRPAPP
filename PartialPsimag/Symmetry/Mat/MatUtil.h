#ifndef PSIMAG_MatUtil_H_
#define PSIMAG_MatUtil_H_

#include "Mat.h"
#include "Vec.h"
#include<stdexcept>



namespace psimag {


//=====================================================================//
//====================== IMPLEMENT SOLVE ==============================//
//=====================================================================//

/** Solve a matrix equation a*x=v for the unknown vector x.
 *  Note that the arguments to the function are organized as
 *  input objects followed by output objects.
 *
 *  \author Gregory Brown
 */
template<typename T, std::size_t NROW, std::size_t NCOL>
void Solve(const Mat<T,NROW,NCOL>& a, 
           const Vec<T,NCOL>& v,
                 Vec<T,NCOL>& x) 
{
  // Still need to wrap a BLAS routine here for general cases
  ;
}


// Partial specialization used to make small solutions faster.
template<typename T>
void Solve(const Mat<T,3,3>& m,
           const Vec<T,3>& v,
                 Vec<T,3>& x)
{
  // Find reciprocal of determinant
  T iDet = Det(m);
  if(iDet==0)
  {
    std::ostringstream msg;
    msg << "MatUtil.h::Solve(Mat<T,3,3>): attempt to solve system with zero determinant";
    throw std::logic_error(msg.str());
    return;
  }
  iDet = static_cast<T>(1)/iDet;
  // This is just inverse matrix times vector
  x[0] = ( ( m[4]*m[8] - m[7]*m[5] ) * v[0] 
          +( m[6]*m[5] - m[3]*m[8] ) * v[1] 
          +( m[3]*m[7] - m[6]*m[4] ) * v[2] ) * iDet;
  //
  x[1] = ( ( m[7]*m[2] - m[1]*m[8] ) * v[0]
          +( m[0]*m[8] - m[6]*m[2] ) * v[1]
          +( m[6]*m[1] - m[0]*m[7] ) * v[2] ) * iDet;
  //
  x[2] = ( ( m[1]*m[5] - m[4]*m[2] ) * v[0]
          +( m[3]*m[2] - m[0]*m[5] ) * v[1]
          +( m[0]*m[4] - m[3]*m[1] ) * v[2] ) * iDet;
  return;
}


// Partial specialization used to make small solutions faster.
template<typename T>
void Solve(const Mat<T,2,2>& m,
           const Vec<T,2>& v,
                 Vec<T,2>& x)
{  
  // Find reciprocal of determinant
  T iDet = Det(m);
  if(iDet==0)
  {
    std::ostringstream msg;
    msg << "MatUtil.h::Solve(Mat<T,2,2>): attempt to solve system with zero determinant";
    throw std::logic_error(msg.str());
    return;
  }
  iDet = static_cast<T>(1)/iDet;
  // using array indexing
  x[0] = iDet * ( m[3]*v[0] - m[2]*v[1] );
  x[1] = iDet * ( m[0]*v[1] - m[1]*v[0] );
  //
  return;
}



}       // namespace psimag

#endif  // PSIMAG_MatUtil_H_
