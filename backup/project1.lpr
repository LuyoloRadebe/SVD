
PROGRAM matrix (Input, Output);
VAR
   A, B, C: ARRAY[1..3, 1..3] OF Integer;
   Arow, AColumn, Brow, BColumn, I, J, K, SUM: Integer;
BEGIN
     Writeln('MATRIX MULTIPLICATION PROGRAM');
     Writeln('-----------------------------');
     Writeln;
     Writeln(' Enter the number of ROWS of Matrix A');
     Readln(Arow);
     Writeln(' Enter the number of COLUMNS of Matrix A');
     Readln(AColumn);
     Writeln(' Enter the number of ROWS of Matrix B');
     Readln(Brow);
     Writeln(' Enter the number of COLUMNS of Matrix B');
     Readln(BColumn);
     Writeln();
     IF (AColumn = Brow) THEN
        BEGIN
             Writeln(' Enter the elements of Matrix A');
             For I:= 1 to 3 DO
                 For J:= 1 TO 3 DO
                     Readln(A[I,J]);
             Writeln(' Enter the elements of Matrix B');
             For I:= 1 to 3 DO
                 For J:= 1 TO 3 DO<
                     Readln(B[I,J]);
             Writeln;
{ Display the elements of Matrix A }
        Writeln(' MATRIX A');
        Writeln(' --------');
        For I:= 1 TO 3 DO
            BEGIN
                 For J:= 1 TO 3 DO
                     Write(A[I,J]:3);
                     Writeln;
            END;
        Writeln;
{ Display the elements of Matrix B }
        Writeln(' MATRIX B');
        Writeln(' --------');
        For I:= 1 TO 3 DO
            BEGIN
                 For J:= 1 TO 3 DO
                 Write(B[I,J]:3);
                 Writeln;
            END;
        Writeln;
{ Compute the Product of Matrices A and B }
        For I:= 1 to 3 DO
            For J:= 1 TO 3 DO
                BEGIN
                     SUM:=0;
                     For K:= 1 TO 3 DO
                     SUM:= SUM + A[I,K] * B[K,J];
                     C[I,J]:= SUM;
                END;
        Writeln;
{ Display the Product of Matrices A and B }
        Writeln(' MATRIX AB');
        Writeln(' ---------');
        Writeln;
        For I:= 1 TO 3 DO
            BEGIN
                 For J:= 1 TO 3 DO
                 Write(C[I,J]:5);
                 Writeln;
            END;
        END
     ELSE
         Writeln(' Matrix A and B cannot be multiplied. Enter Matrices in the order m x n and n x p');
         Writeln();
         Writeln(' PRESS THE ENTER KEY TO EXIT');
         Readln;
END.

