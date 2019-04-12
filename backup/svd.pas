unit SVD;

{$mode objfpc}{$H+}

interface
uses
Classes, SysUtils, Math, MatrixOperations, uComplex, Dos;
type
DynamicMatrix = array of array of Real;
function SingularValues(Matrix : DynamicMatrix) : DynamicMatrix;
function TridiagonalHouseholder(Matrix : DynamicMatrix) : DynamicMatrix;
function QRDecomposition(Matrix : DynamicMatrix) : DynamicMatrix;
function Eigenvalue(Matrix : DynamicMatrix) : DynamicMatrix;
function NewtonRaphsonSecular(SubEigen, Z: DynamicMatrix) : DynamicMatrix;

implementation
function SingularValues(Matrix : DynamicMatrix) : DynamicMatrix;
var
Transpose : DynamicMatrix;
time_elapsed, chours, cminutes, cseconds, cmilliseconds, hours, minutes, seconds, milliseconds : word;
begin
GetTime(hours, minutes, seconds, milliseconds);
Transpose := MatrixTranspose(Matrix);
Matrix := MatrixMultiply(Matrix,Transpose);  //Will be symmetric
Matrix := TridiagonalHouseholder(Matrix);
Result := Eigenvalue(Matrix);
Result := MatrixSqrt(Result);
GetTime(chours, cminutes, cseconds, cmilliseconds);
time_elapsed := ((cmilliseconds-milliseconds) + (cseconds-seconds)*60 + (cminutes-minutes)*3600 + (chours-hours)*216000);
chours := time_elapsed div 216000;
time_elapsed := time_elapsed - chours*216000;
cminutes := time_elapsed div 3600;
time_elapsed := time_elapsed - cminutes*216000;
cseconds := time_elapsed div 60;
time_elapsed := time_elapsed - cseconds*216000;
cmilliseconds := time_elapsed;
writeln(cminutes,' minutes, ',cseconds, ' seconds, ',cmilliseconds, ' milliseconds');
end;
function TridiagonalHouseholder(Matrix : DynamicMatrix) : DynamicMatrix;
var
  I, J, K                           : Integer;
  S, R, Answer                      : Double;
  A, W, WT, V, C, Q, QT, Hold       : DynamicMatrix;
begin
     //Dimension Initialization
     SetLength(A,High(Matrix)+1,High(Matrix)+1);
     SetLength(Result,High(Matrix)+1,High(Matrix)+1);
     SetLength(Hold,High(Matrix)+1,High(Matrix)+1);
     SetLength(WT,1,High(Matrix)+1);
     SetLength(W,High(Matrix)+1,1);
     SetLength(V,High(Matrix),1);
     SetLength(QT,High(Matrix),1);
     SetLength(Q,High(Matrix),1);
     SetLength(C,1,1);

     for I := Low(Matrix) to High(Matrix)-2 do
     begin
     S := 0;
     R := 0;

     for J := Low(Matrix)+I+1 to High(Matrix) do
     begin
     S := S + power(Matrix[J,I],2);
     end;
     S := ((Matrix[Low(Matrix)+I+1,I])/abs((Matrix[Low(Matrix)+I+1,I])))*sqrt(S);
     R := sqrt(2*S*(Matrix[Low(Matrix)+I+1,I] + S));
     for J:= Low(Matrix) to Low(Matrix)+I do
     begin
     WT[0,J] := 0;
     end;
     WT[0,Low(Matrix)+I+1] := (S+Matrix[Low(Matrix)+I+1,I])/R;
     for J:= Low(Matrix)+I+2 to High(Matrix) do
         begin
         WT[0,J] := (Matrix[I,J])/R;
         end;
     W := MatrixTranspose(WT);
     V := MatrixMultiply(Matrix,W);
     C := MatrixMultiply(WT,V);
     Q := MatrixMultiply(W,C);
     Q := MatrixSubtract(V,Q);
     A := MatrixMultiply(Q,WT);
     A := MatrixMultiplyConst(A,2);
     QT := MatrixTranspose(Q);
     Hold := MatrixMultiply(W,QT);
     Hold := MatrixMultiplyConst(Hold,2);
     Hold := MatrixAdd(Hold,A);
     A := MatrixSubtract(Matrix,Hold);
     Matrix := A;
     end;
     Result := Matrix;
end;
function QRDecomposition(Matrix : DynamicMatrix) : DynamicMatrix;
var
  I, J, K                           : Integer;
  UNorm, Constant, Answer                      : Double;
  A, U, E, W, WT, V, C, Q, R, QT, Hold       : DynamicMatrix;
begin
    //Dimension Initialization
    SetLength(A,High(Matrix)+1,1);
    SetLength(U,High(Matrix)+1,1);
    SetLength(Hold,High(Matrix)+1,1);
    SetLength(E,1,High(Matrix)+1);
    SetLength(Q,High(Matrix)+1,High(Matrix)+1);
    SetLength(Result,High(Matrix)+1,High(Matrix)+1);
    SetLength(R,High(Matrix)+1,High(Matrix)+1);
    SetLength(C,1,1);

    for I := Low(Matrix) to High(Matrix) do
    begin
         for J:= Low(Matrix) to High(Matrix) do
         begin
         A[J,0]:= Matrix[J,I];
         end;
    J := I;
    Hold := MatrixSubtract(Hold,Hold);
         while J > 0 do
         begin
              for K := Low(Matrix) to High(Matrix) do
              begin
              E[0,K]:= Q[K,J-1];
              end;
         U := MatrixMultiply(E,A);
         U := MatrixMultiply(U,E);
         U := MatrixTranspose(U);
         Hold := MatrixAdd(Hold,U);
         J := J - 1;
         end;
    U := MatrixSubtract(A,Hold);
    UNorm := MatrixNorm(U);
    UNorm := 1/UNorm;
    E := MatrixMultiplyConst(U,UNorm);
    E := MatrixTranspose(E);
         for J:= Low(Matrix) to High(Matrix) do
         begin
         Q[J,I]:= E[0,J];
         end;
    end;
  Result := Q;
end;
function Eigenvalue(Matrix : DynamicMatrix) : DynamicMatrix;
var
  I, J, K, Count, Limit, EigenCount                                   : Integer;
  Answer, Tol, F, FPrime, Guess,NewGuess, OldGuess, Diff                           : Double;
  U : DynamicMatrix;
  Block,S,Q,QPrime,SubEigen,R,T,Hold,Z,Eigen      : DynamicMatrix;
begin
  if (High(Matrix)+1) mod 2 = 1  then
     begin
     Matrix := MatrixConcatenate(Matrix,Hold,0);
     end;

  WriteLn('In Eigenvalue Function:');
   MatrixPrint(Matrix);
  I := (High(Matrix)+1) div 2;
  SetLength(U,High(Matrix)+1,1);
  U[I,0] := sqrt(abs(Matrix[I-1,I]));
  U[I-1,0] := U[I,0];
  WriteLn('Matrix U');
   MatrixPrint(U);
  //S MATRIX
  S := MatrixTranspose(U);
  S := MatrixMultiply(U,S);
 WriteLn('Matrix S');
   MatrixPrint(S);
  //BLOCK MATRIX
  Block := MatrixSubtract(Matrix,S);
  WriteLn('Block Matrix');
   MatrixPrint(Block);
  for J:= 0 to 1 do
  begin
  T := MatrixCopy((J)*(I),(J+1)*(I)-1,(J)*(I),(J+1)*(I)-1,Block);
  //Tn MATRIX
  WriteLn('T',J);
   MatrixPrint(T);
  Q := MatrixIdentity(High(T[0])+1);
 WriteLn('Q',J);
   MatrixPrint(Q);
  Count := 0;
  Limit := 10000;
  wRITElN('Matrix T');
  MatrixPrint(T);
  while ((abs(T[0,1])) > Tol) and (Count < Limit) do
    begin
    write('.');
    Count := Count + 1;
    QPrime := QRDecomposition(T);
    T := MatrixDiv(T,QPrime);
    T := MatrixMultiply(T,QPrime);
    Q := MatrixMultiply(Q,QPrime);

    end;
  if J = 0 then
    begin
    SubEigen := T;
    R := Q;
    WriteLn('SubEigen Matrix');
    MatrixPrint(SubEigen);
    WriteLn('R Matrix');
    MatrixPrint(R);
    end;
  end;
  SubEigen := MatrixConcatenate(SubEigen,T,2);
  R := MatrixConcatenate(R,Q,2);
  //R MATRIX
  WriteLn('R Matrix');
  MatrixPrint(R);
  R := MatrixTranspose(R);
  //U MATRIX
  WriteLn('U Matrix');
  MatrixPrint(U);
  //Z MATRIX
  Z := MatrixMultiply(R,U);
  Eigen := NewtonRaphsonSecular(SubEigen, Z);
  Result := Eigen;
  end;
function NewtonRaphsonSecular(SubEigen, Z : DynamicMatrix) : DynamicMatrix;
var
I, J, K, Flag, Count, Limit, EigenCount,Skip              : Integer;
Answer, Tol, F, FPrime, FPPrime, Guess,NewGuess, OldGuess, Diff,Repeated, Scale : Double;
Block,S,Q,QPrime,R,U,T,Hold,Eigen      : DynamicMatrix;
begin
  SetLength(Eigen,High(SubEigen)+1,High(SubEigen)+1);
  //Starting Guess
  Skip := 0;
  NewGuess := 0;
  Flag := 0;
  J := 0;
  EigenCount := 0;
  Repeated := 0;
  Tol := 1e-12;
  WriteLn('Z MATRIX');
  MatrixPrint(Z);
  WriteLn('SUB EIGEN');
  MatrixPrint(SubEigen);
  //Eigenvalue loop
  while EigenCount < High(SubEigen)+1 do
    begin
    Count := 0;
    Diff := 1;
    Scale := Factor(NewGuess);
    //Newton-Raphson Loop
    while (abs(Diff) > Tol) and (Count < 10000) do
      begin
      OldGuess := NewGuess;
      F := Scale;
      FPrime := 0;
      FPPrime := 0;
      for I := 0 to High(SubEigen) do
        begin
        F := F + Scale*((power(Z[I,0],2))/((SubEigen[I,I])-OldGuess));
        FPrime := FPrime + Scale*(power(Z[I,0],2))/(power((SubEigen[I,I] - OldGuess),2));
        FPPrime := FPPrime + 2*(power(Z[I,0],2))/(power((SubEigen[I,I] - OldGuess),3));
        end;
      NewGuess := OldGuess - ((F)/(FPrime));
      //WriteLn('NewGuess := OldGuess - ((F)/(FPrime))', NewGuess,' := ',OldGuess,' - ((',F,',)/(',FPrime,')))');
      Count := Count + 1;
      Diff := OldGuess-NewGuess;
      if (NewGuess < 0) then
        begin
        Diff := 0;
        Flag := 1;
        end;
      end;

    //Conditions if max iterations haven't been reached
    if (Count < 10000) and (Flag <> 1) then
      begin
      WriteLn('EigenValue #',EigenCount,' = ',NewGuess:0:2,#9,'Number of iterations: ',Count);
      Eigen[EigenCount,EigenCount] := NewGuess;
      EigenCount := EigenCount +1;
      if (SubEigen[High(Subeigen)-EigenCount,High(Subeigen)-EigenCount]) < 0 then
        begin
        Skip := Skip+1;
        end;

      end;

      if (EigenCount) < High(SubEigen)+1 then
        begin
        Writeln('Skip',Skip);
        NewGuess := (SubEigen[High(Subeigen)-EigenCount-Skip,High(Subeigen)-EigenCount-Skip]);
        WriteLn('New Guess Reference',NewGuess:0:3,#9,'New Guess',NewGuess + Tol);
        NewGuess := NewGuess + Tol;
        end;
      Flag := 0;
      end;
  Result := Eigen;
  WriteLn('Tridiagonal EigenValues');
  MatrixPrint(Eigen);
end;
end.

