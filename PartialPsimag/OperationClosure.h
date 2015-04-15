//-*- mode: C++; -*-

/** \ingroup PSIMAG */
/*@{*/

/*! \file OperationClosure.h  
 *
 */

#ifndef PSIMAG_Operation_Closures_H
#define PSIMAG_Operation_Closures_H


namespace psimag {
  
  class OP {
  public:
    typedef enum{FourierTransform,PLUS,MINUS,TIMES,DIVIDE,INV} Type;
  };

  //====================================================================== 
  
  template<typename Operand1Type,
	   int      Operator,
	   typename Operand2Type>
  class OperationClosure {
  public:
    
    const Operand1Type& lhs;
    const Operand2Type& rhs;
    
    OperationClosure(const Operand1Type& l,
		     const Operand2Type& r):
      lhs(l),
      rhs(r)
    {}
  };      

  //====================================================================== 
  
  template<int      Operator,
	   typename OperandType>
  class UnaryOperationClosure {
  public:
    
    const OperandType& operand;
    
    UnaryOperationClosure(const OperandType& theOperand):
      operand(theOperand)
    {}
  };      

} // end namespace rpa


/*@}*/
#endif
