PROGRAM SVD (Input, Output);

uses matrix;

type
  TestMatrix = array of array of integer;
var

function MultiplyMatrix(var arr: TestMatrix) : integer;
begin
        MatrixA[1,1] = MatrixB[1,1];
end;

VAR
   A, B, C: ARRAY[1..3, 1..3] of Integer;
   M, U, E, V: array of array of integer;
   Rows, Columns, I, J, K, SUM: Integer;
BEGIN
     Writeln('- - - Singular Value Decomposition - - -');
     Writeln('CHsARACTERISTIC EQUATION : M = UEV*');
     Writeln('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
     Writeln;
     Writeln(' Enter the number of ROWS of Matrix M');
     Readln(Rows);
     Writeln(' Enter the number of COLUMNS of Matrix M');
     Readln(Columns);
        BEGIN
             setLength(M, Rows, Columns);
             setLength(U, Rows, Rows);
             setLength(E, Rows, Columns);
             setLength(V, Columns, Columns);
             Writeln(' Enter the elements of Matrix M');
             For I:= 0 to (Rows-1) DO
                 BEGIN
                 For J:= 0 TO (Columns-1) DO
                     Readln(M[I,J]);
                     Writeln;
                 END;
{ Display the elements of Matrix M }
        Writeln(' MATRIX M');
        Writeln(' - - - - - - - -');
        For I:= 0 TO (Rows-1) DO
            BEGIN
                 For J:= 0 TO (Columns-1) DO
                     Write(M[I,J]:Columns-1);
                     WriteLn;
            END;
        WriteLn;
        Write(M[I,J]:Columns);
        M = MultiplyMatrix(M,M);

{ Transform the matrix into a bidiagonal, tridiagonal or hessenberg }
{ Perform a transpose and calculate the relevant eigenvalues }
{ Orthonormalize the resulting matrices (MTM and MMT) to find U and V*}
{ Perform a backwards operation to solve for Sigma E}
{ Return all relevant matrices}
        END;
        END.

