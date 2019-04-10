{Main declaration}
program Demo;
{Link to external units similar to libraries in c}
uses
   SysUtils, MatrixOperations, Math, SVD;
{Declare any types if need be}
type
   DynamicMatrix = array of array of double;
{Main variables}
var
   Determinant: double;
   A, B, C: array [1..3, 1..3] of integer;
   AnswerMatrix, M, U, E, V: array of array of double;
   Rows, Columns, I, J, K, L: integer;
   AnswerDouble : double;
{Main begin}
begin
     Writeln('Matrix Operation test and demo');
     WriteLn();
     ReadLn();
     {Writeln(' Enter the number of ROWS of Matrix M');
     Readln(Rows);
     Writeln(' Enter the number of COLUMNS of Matrix M');
     Readln(Columns);}
        begin
        SetLength(M, 5, 5);
        SetLength(U, Rows, Rows);
        SetLength(E, 2, 2);
        SetLength(V, Columns, Columns);

        {Writeln(' Enter the elements of Matrix M');
             for I:= 0 to (Rows-1) do
             begin
                 for J:= 0 TO (Columns-1) do
                 begin
                 Readln(M[I,J]);
                 Writeln;
                 end;
             end;}


        end;


M[1,1] := 12;
M[1,2] := 10;
M[1,3] := 9;
M[1,0] := 2;
M[1,4] := 15;
M[2,1] := 0;
M[2,2] := 32;
M[2,3] := 12;
M[2,0] := 11;
M[2,4] := 2;
M[4,1] := 3;
M[4,2] := 17;
M[4,3] := 20;
M[4,4] := 16;
M[4,0] := 22;
M[3,1] := 6;
M[3,2] := 14;
M[3,3] := 81;
M[3,4] := 19;
M[3,0] := 13;
M[0,1] := 13;
M[0,2] := 4;
M[0,3] := 21;
M[0,4] := 6;
M[0,0] := 18;
{

M[0,0] := 3;
M[0,1] := 0;
M[0,2] := 2;
M[1,0] := 2;
M[1,1] := 0;
M[1,2] := -2;
M[2,0] := 0;
M[2,1] := 1;
M[2,2] := 1;
 }

E[0,0] := 4;
E[0,1] := 1;
E[1,0] := -2;
E[1,1] := 3;

// EXPECTED ANSWERS
// MATRIX M: 10
// MATRIX E: 14
//AnswerDouble :=  MatrixDeterminant(E);

AnswerDouble := MatrixDeterminant(E);
WriteLn('E Answer Double:',AnswerDouble);

MatrixPrint(M);
M := MatrixInverse(M);
MatrixPrint(M);

ReadLn;
{Display the elements of Matrix M}
{Multiply Matrix M with Matrix E}
{WriteLn('Input Matrix M:');
MatrixPrint(M);
WriteLn();
ReadLn;
U := MatrixInverse(M);
WriteLn('Matrix M Matrix of Minors:');
MatrixPrint(U);
WriteLn();
ReadLn;}
{ Transform the matrix into a bidiagonal, tridiagonal or hessenberg }
{ Perform a transpose and calculate the relevant eigenvalues }
{ Orthonormalize the resulting matrices (MTM and MMT) to find U and V*}
{ Perform a backwards operation to solve for Sigma E}
{ Return all relevant matrices}
end.

