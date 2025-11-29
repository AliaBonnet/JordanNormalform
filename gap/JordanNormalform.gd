DeclareGlobalFunction("ConvertVecToRowMat");

DeclareGlobalFunction("GenerateRandomVector");

DeclareGlobalFunction("SpinningCheck");

DeclareGlobalFunction("SpinUntil");

DeclareGlobalFunction("FindVectorNotInSubspaceNC");

DeclareGlobalFunction("FindCyclicVectorNC");

DeclareGlobalFunction("RemoveZeroRows");

DeclareGlobalFunction("MinPolVec");

DeclareGlobalFunction("PolyEvalFromSpan");

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
#!  Returns a matrix <M>B</M> such that <A>A</A>^(<M>B^-1</M>) is the Jordan 
#!  normal form of <A>A</A>. The algorithm first computes a primary decomposition
#!  of <A>A</A> following a modified version of Steel's algorithm and then 
#!  computes a cyclic decomposition of the primary components. Finally it computes 
#!  Jordan block form for each of the cyclic components. It works for matrices 
#!  over any field that is available in GAP.
#! 
#! @BeginExampleSession
#! gap> A:=[ [  2,  2,  0,  1,  0,  2,  1 ],
#! >         [  0,  4,  0,  0,  0,  1,  0 ],
#! >         [  0,  1,  1,  0,  0,  1,  1 ],
#! >         [  0, -1,  0,  1,  0, -1,  0 ],
#! >         [  0, -7,  0,  0,  1, -5,  0 ],
#! >         [  0, -2,  0,  0,  0,  1,  0 ],
#! >         [  0, -1,  0,  0,  0, -1,  1 ] ];
#! gap> f:=FrobeniusNormalForm(A);
#! [ [ x_1^4-7*x_1^3+17*x_1^2-17*x_1+6, x_1^2-3*x_1+2, x_1-1 ], 
#!                                  # f[1] = List of invariant factors
#!   [ [    1,   -2,    1,    1,    0,    0,    1 ],
#!     [    2,   -7,    1,    2,    0,   -1,    3 ],
#!     [    4,  -26,    1,    4,    0,   -8,    6 ],
#!     [    8,  -89,    1,    8,    0,  -35,   11 ],
#!     [ -1/2,   -2,    0,  1/2,    0,   -2, -3/2 ],
#!     [   -1,   -4,    0,    0,    0,   -4,   -2 ],
#!     [    0,  9/4,    0,   -3,    1,  5/4,  1/4 ] ],
#!                                  # f[2] = base change matrix P
#!   [ 1, 5, 7 ]  ]                 # f[3] = indices where the blocks begin
#! gap> PrintArray(f[2]*A*f[2]^-1);
#! [ [   0,   1,   0,   0,   0,   0,   0 ], 
#!   [   0,   0,   1,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   1,   0,   0,   0 ],
#!   [  -6,  17, -17,   7,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   1,   0 ],
#!   [   0,   0,   0,   0,  -2,   3,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   1 ] ]
#! (This is the Frobenius normal form; there are 3 diagonal blocks,
#!  one of size 4, one of size 2 and one of size 1.)
#! @EndExampleSession
#! 
#! You can also use  'CreateNormalForm( f[1] );' to produce the Frobenius
#! normal form. (This function just builds the block diagonal matrix with 
#! diagonal blocks given by the companion matrices corresponding to the 
#! various invariant factors of <A>A</A>.) 
DeclareGlobalFunction("JordanNormalform");
