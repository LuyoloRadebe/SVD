unit MatrixOperations;

{$mode objfpc}{$H+}
interface
uses
Classes, SysUtils, Math;
type
DynamicMatrix = array of array of double;

procedure     MatrixPrint         (Matrix:DynamicMatrix);
function      MatrixMultiplyConst (Matrix:DynamicMatrix ; Constant : Double) : DynamicMatrix;
function      MatrixMultiply      (Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
function      MatrixDiv           (Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
function      MatrixAdd           (Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
function      MatrixSubtract      (Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
function      MatrixSqrt          (Matrix:DynamicMatrix) : DynamicMatrix;
function      MatrixDeterminant   (Matrix:DynamicMatrix):Double;
function      MatrixAdjugate      (Matrix:DynamicMatrix ; Index1, Index2:integer) : DynamicMatrix;
function      MatrixTranspose     (Matrix:DynamicMatrix):DynamicMatrix;
function      MatrixInverse       (Matrix:DynamicMatrix):DynamicMatrix;
function      MatrixNorm          (Matrix:DynamicMatrix):Double;
function      MatrixCopy          (MBegin,MEnd,NBegin,NEnd: Integer;Matrix:DynamicMatrix):DynamicMatrix;
function      MatrixPaste         (MBegin,NBegin: Integer;Mat1,Mat2:DynamicMatrix):DynamicMatrix;
function      MatrixConcatenate   (Mat1,Mat2: DynamicMatrix; Mode: Integer):DynamicMatrix;
function      MatrixIdentity      (Size: Integer):DynamicMatrix;
function      Factor              (Number: Double):Double;
implementation
procedure MatrixPrint(Matrix : DynamicMatrix);
var
   I,J: integer;
begin
  for I:= (low(Matrix)) to (high(Matrix)) do
      begin
           for J:= (low(Matrix[0])) to (high(Matrix[0])) do
               Write(Matrix[I,J]:1:5,0,#9:high(Matrix));
               WriteLn();
      end;
end;
function  Factor(Number: Double):Double;
var
I,J: Double;
begin
I := 10;
J := 1;
while Number > 0 do
  begin
  Number := Number / 10;
  J := J*10;
  end;
J := J/10;
Result := J;
end;
function MatrixSqrt(Matrix:DynamicMatrix) : DynamicMatrix;
var
   I,J: integer;
begin
for I:= Low(Matrix) to High(Matrix) do
    begin
    for J:= Low(Matrix[0]) to High(Matrix[0]) do
        begin
        Result[I,J] := power(Matrix[I,J],0.5);
        end;
    end;
end;
function MatrixMultiplyConst(Matrix:DynamicMatrix;Constant:Double):DynamicMatrix;
var
I,J: integer;
begin
for I:= Low(Matrix) to High(Matrix) do
    begin
    for J:= Low(Matrix[0]) to High(Matrix[0]) do
        begin
        Result[I,J] := Constant* Matrix[I,J];
        end;
    end;
end;
function MatrixMultiply(Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
var
   I,J,K,L: integer;
   Hold   : DynamicMatrix;
begin
SetLength(Hold, high(Mat1)+1,high(Mat2[0])+1);
SetLength(Result, high(Mat1)+1,high(Mat2[0])+1);

  if High(Mat1[0]) = High(Mat2) then
  begin
  for I:= Low(Mat1) to High(Mat1) do
      begin
      for J:= Low(Mat2[0]) to High(Mat2[0]) do
          begin
          Result[I,J] := 0;
          for K:= Low(Mat2) to High(Mat2) do
              Result[I,J] := Result[I,J] + (Mat1[I,K]*Mat2[K,J]);
          end;
      end;
  end;
  if High(Mat1[0]) <> High(Mat2) then
  WriteLn('Error: Matrix dimensions must agree');
end;
function MatrixDiv(Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
begin
SetLength(Result,(High(Mat1)+1),(High(Mat1[0])+1));
Mat2 := MatrixInverse(Mat2);
Result := MatrixMultiply(Mat2,Mat1);
end;
function MatrixDeterminant (Matrix:DynamicMatrix):Double;
var
   MatHold: DynamicMatrix;
   I,J,K,notbusy, count : integer;
   Sign, PosMult, NegMult, PosSum, NegSum, CofactorHold: Double;
   Det,DetHold : Double;
begin
  if High(Matrix[0]) <> High(Matrix[0]) then
     begin
     WriteLn('Error: Matrix dimensions must be square');
     end;
  if High(Matrix[0]) = 0 then
     begin
     Result := Matrix[0,0];
     end
  else
     begin
     Det:= 0;
     Sign := 1;
     for K:= Low(Matrix) to High(Matrix) do
          begin
          MatHold := MatrixAdjugate(Matrix,0,K);

          if High(Matrix) = 1 then
          CofactorHold := MatHold[0,0]
          else
          CofactorHold := MatrixDeterminant(MatHold);
          Det := Det + Sign*Matrix[0,K]*CofactorHold;
          Sign := Sign * (-1);
         end;
     Result := Det;
     end;
end;
function MatrixAdjugate(Matrix:DynamicMatrix ; Index1, Index2:integer) : DynamicMatrix;
var
   Size,I,J,K,Flag1, Flag2: integer;
begin
  SetLength(Result,high(Matrix),high(Matrix));
  Flag1 := 0;
  Flag2 := 0;
  for I:= Low(Matrix) to High(Matrix)-1 do
      begin
      if I = Index1 then
      Flag1 := 1;
      for J:= Low(Matrix) to High(Matrix)-1 do
          begin
          if J = Index2 then
          Flag2 := 1;
          if I+Flag2 < High(Matrix)+1 then
          Result[I,J] := Matrix[I+Flag1,J+Flag2];
          end;
      Flag2 := 0;
      end;
  Flag1 := 0;
end;
function MatrixTranspose(Matrix:DynamicMatrix) : DynamicMatrix;
var
I,J,K : integer;
begin
     SetLength(Result,(High(Matrix[0])+1),(High(Matrix)+1));
  for I:= Low(Matrix[0]) to High(Matrix[0]) do
      begin
      for J:= Low(Matrix) to High(Matrix) do
          begin
          Result[I,J] := Matrix[J,I];
          end;
      end;
end;
function MatrixInverse(Matrix :DynamicMatrix) : DynamicMatrix;
var
I,J,K : integer;
Det,Answer : double;
AnswerMatrix : DynamicMatrix;
begin
SetLength(Result,High(Matrix)+1,High(Matrix)+1);
SetLength(AnswerMatrix,High(Matrix),High(Matrix));
Det := MatrixDeterminant(Matrix);
for I:= Low(Matrix) to High(Matrix) do
    begin
    for J:= Low(Matrix) to High(Matrix) do
        begin
          AnswerMatrix := MatrixAdjugate(Matrix,I,J);
          Answer := MatrixDeterminant(AnswerMatrix);
          Result[I,J] := (power(-1,I)*power(-1,J)*Answer)/Det;
        end;
    end;
Result := MatrixTranspose(Result);
end;
function MatrixAdd(Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
var
I,J,K : integer;
begin
SetLength(Result,(High(Mat1)+1),(High(Mat1[0])+1));
for I:= Low(Mat1) to High(Mat1) do
    begin
    for J:= Low(Mat1[0]) to High(Mat1[0]) do
        begin
        Result[I,J] := Mat1[I,J] + Mat2[I,J];
        end;
    end;
end;
function MatrixSubtract(Mat1,Mat2:DynamicMatrix) : DynamicMatrix;
var
I,J,K : integer;
begin
SetLength(Result,(High(Mat1)+1),(High(Mat1[0])+1));
for I:= Low(Mat1) to High(Mat1) do
    begin
    for J:= Low(Mat1[0]) to High(Mat1[0]) do
        begin
        Result[I,J] := Mat1[I,J] - Mat2[I,J];
        end;
    end;
end;
function MatrixNorm (Matrix:DynamicMatrix):Double;
var
I,J,K : integer;
begin
Result := 0;
for I:= Low (Matrix) to High(Matrix) do
    begin
    for J := Low(Matrix[0]) to High(Matrix[0]) do
        begin
        Result := Result + power(Matrix[I,J],2);
        end;
    end;
Result := sqrt(Result);
end;
function MatrixCopy (MBegin,MEnd,NBegin,NEnd: Integer;Matrix:DynamicMatrix):DynamicMatrix;
var
I,J,K : integer;
begin
SetLength(Result, (MEnd-MBegin)+1,(NEnd-NBegin)+1);
for I:= 0 to (MEnd-MBegin) do
    begin
    for J := 0 to (NEnd-NBegin) do
        begin
        Result[I,J] := Matrix[MBegin+I,NBegin+J];
        end;
    end;
end;
function MatrixPaste (MBegin,NBegin: Integer;Mat1,Mat2:DynamicMatrix):DynamicMatrix;
var
I,J,K : integer;
begin
SetLength(Result, High(Mat1)+1,High(Mat1[0])+1);
for I:= 0 to High(Mat2) do
    begin
    for J := 0 to High(Mat2[0]) do
        begin
        Result[MBegin+I,NBegin+J] := Mat2[I,J];
        end;
    end;
end;
function MatrixConcatenate (Mat1,Mat2: DynamicMatrix; Mode: Integer):DynamicMatrix;
var
I,J,K : integer;
begin
case Mode of
0: //horizontal with mat1 left and mat2 right
begin
end;
1: //vertical with mat1 top and mat2 bot
begin
end;
2: //diagonal with mat1 top left and mat 2 bot right
begin
SetLength(Result, High(Mat1)+High(Mat2)+2,High(Mat1[0])+High(Mat2[0])+2);
Result := MatrixPaste(0,0,Result,Mat1);
Result := MatrixPaste(High(Mat1)+1,High(Mat1[0])+1,Result,Mat2);
end;
end;
end;
function MatrixIdentity (Size: Integer):DynamicMatrix;
var
I,J,K : integer;
Identity : DynamicMatrix;
begin
SetLength(Identity, Size,Size);
for I := 0 to Size-1 do
    begin
        Identity[I,I] := 1;
    end;
Result := Identity;
end;
end.
