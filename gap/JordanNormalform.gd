DeclareGlobalFunction("ConvertVecToRowMat");

DeclareGlobalFunction("GenerateRandomVector");

DeclareGlobalFunction("SpinUntil");

DeclareGlobalFunction("FindVectorNotInSubspaceNC");

DeclareGlobalFunction("FindCyclicVectorNC");

DeclareGlobalFunction("RemoveZeroRows");

DeclareGlobalFunction("MinPolVec");

DeclareGlobalFunction("PolyEvalFromSpan");

#! @Arguments A
#! @Description
#!  Returns a base change matrix <M>B</M> such that <M>B</M><A>A</A><M>B^{-1}</M> is a
#!  primary form of <A>A</A>.
#!  This function uses a modified version of Steel's algorithm.
#! 
#! @BeginExampleSession
#! gap> A := [ [ Z(5)^2, Z(5)^2, Z(5)^2, Z(5)^3, Z(5)^0, Z(5)^2 ], 
#! [ Z(5)^0, Z(5)^2, Z(5), Z(5)^0, Z(5)^0, Z(5) ], 
#!  [ Z(5)^2, Z(5)^2, Z(5)^0, 0*Z(5), Z(5)^2, Z(5)^0 ], 
#!  [ Z(5), Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5), Z(5) ], 
#!  [ Z(5)^3, 0*Z(5), Z(5)^0, Z(5)^0, Z(5)^3, Z(5)^0 ], 
#!  [ 0*Z(5), Z(5)^2, Z(5), Z(5), Z(5)^2, Z(5)^0 ] ]
#! gap> B := PrimaryDecomp(A);
#! < mutable compressed matrix 6x6 over GF(5) >
#! gap> Display(A^Inverse(B));
#!  . 1 . . . .
#!  . . 1 . . .
#!  3 . 4 . . .
#!  . . . . 1 .
#!  . . . . . 1
#!  . . . 3 3 3
#! @EndExampleSessiongap> B := PrimaryDecomp(A);
#!< mutable compressed matrix 6x6 over GF(5) >
#!gap> Display(A^Inverse(B));
#! . 1 . . . .
#! . . 1 . . .
#! 3 . 4 . . .
#! . . . . 1 .
#! . . . . . 1
#! . . . 3 3 3

DeclareGlobalFunction("PrimaryDecomp");

DeclareGlobalFunction("PrimaryDecompositionforJNF");

DeclareGlobalFunction("PrimaryDecompositionforJNFCyclic");

DeclareGlobalFunction("GetMinPolPowerWithVec");

DeclareGlobalFunction("FindLinearDependenceNC");

DeclareGlobalFunction("CyclicDecompositionOfPrimarySubspace");

DeclareGlobalFunction("JordanBlock");

DeclareGlobalFunction("JordanNormalformIrred");

#! @Arguments A
#! @Description
#!  Returns a base change matrix <M>B</M> such that <M>B</M><A>A</A><M>B^{-1}</M> is the Jordan 
#!  normal form of <A>A</A>. The algorithm first computes a primary decomposition
#!  of <A>A</A> following a modified version of Steel's algorithm and then 
#!  computes a cyclic decomposition of the primary components. Finally it computes 
#!  Jordan block form for each of the cyclic components. It works for matrices 
#!  over finite fields. 
#! 
#! @BeginExampleSession
#! gap> A := [ [ 0*Z(5), 0*Z(5), Z(5)^3, Z(5)^3, Z(5)^3, Z(5)^0 ], 
#! [ 0*Z(5), Z(5)^2, Z(5)^2, Z(5)^0, Z(5)^3, Z(5)^3 ], 
#! [ Z(5)^0, Z(5)^0, Z(5)^3, Z(5)^2, Z(5)^0, Z(5) ], 
#! [ 0*Z(5), Z(5)^3, Z(5), Z(5), 0*Z(5), Z(5)^2 ], 
#! [ Z(5)^2, Z(5)^0, Z(5)^0, 0*Z(5), Z(5), Z(5) ], 
#! [ 0*Z(5), Z(5)^0, Z(5)^2, Z(5), Z(5), Z(5) ] ];;
#! gap> B := JordanNormalform(A);
#! < mutable compressed matrix 6x6 over GF(5) >
#! gap> Display(A^Inverse(B));
#! 3 . . . . .
#! . 1 . . . .
#! . . . 1 . .
#! . . 2 . . .
#! . . . . . 1
#! . . . . 3 4
#! @EndExampleSession

DeclareGlobalFunction("JordanNormalform");
