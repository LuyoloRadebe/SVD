{System Eigenvalues:   { -3.645502505725471 ;  1.3542211142325042 ;  1.645619586029866 ;  6.645661805463101 }}
{System Eigenvalues QR Transform:   { -3.6457513110645894 ;  1.3542486889354095 ;  1.6457513110645905 ;  6.645751311064584 }}

{Real Eigenvalues:   { -4.249977875683766 ;  -2.7396032268707424 ;  2.9073897948309457 ;  12.082191307723573 }}

{Submatrix Eigenvalues: Real Eigenvalues:   { -2.5386983887223633 ;  5.376420728553985 }}
{Main declaration}
program Demo;
{Link to external units similar to libraries in c}
uses
   SysUtils, MatrixOperations, Math, SVD;
{Declare any types if need be}
type
   DynamicMatrix = array of array of double;
   //ComplexDynamicMatrix = array of array of Complex;
const
Tol = 1e-10;
{Main variables}
var
   Determinant: double;
   A, B, C: array [1..3, 1..3] of integer;
   Tridiagonal, Eigen, AnswerMatrix, Am, Bm,QPrime, Cm, M, T, U, E,G, V, Q, Q2, Q3,
      R, ZT, Z, ZFULL: array of array of double;
   Hold : array of array of double;
   Rows, Columns, Count, I, J, K, L: integer;
   AnswerDouble : double;
{Main begin}
begin
     Writeln('Matrix Operation test and demo v1.2');
     Writeln('Eigenvalue found, method of iteraton found.');
     Writeln('To do: Organize recursion and find eigenvector matrix');
//     ReadLn();


{SetLength(M, 14, 14);
M[0,0] := 4;
M[0,1] := 2;
M[0,2] := 0;
M[0,3] := 1;
M[0,4] := 8;
M[0,5] := 9;
M[0,6] := 9;
M[0,7] := 5;
M[0,8] := 7;
M[0,9] := -7;
M[0,10] := 0;
M[0,11] := 1;
M[0,12] := 4;
M[0,13] := -3;
M[1,0] := 2;
M[1,1] := 3;
M[1,2] := -1;
M[1,3] := 6;
M[1,4] := 3;
M[1,5] := 5;
M[1,6] := 9;
M[1,7] := 5;
M[1,8] := 7;
M[1,9] := -7;
M[1,10] := 0;
M[1,11] := 1;
M[1,12] := 10;
M[1,13] := 12;
M[2,0] := 0;
M[2,1] := -1;
M[2,2] := 9;
M[2,3] := 0;
M[2,4] := 11;
M[2,5] := 7;
M[2,6] := 9;
M[2,7] := 5;
M[2,8] := 7;
M[2,9] := -7;
M[2,10] := 0;
M[2,11] := 1;
M[2,12] := 10;
M[2,13] := 2.8;
M[3,0] := 1;
M[3,1] := 6;
M[3,2] := 0;
M[3,3] := 1;
M[3,4] := 12;
M[3,5] := -7;
M[3,6] := 9;
M[3,7] := 5;
M[3,8] := 7;
M[3,9] := -7;
M[3,10] := 0;
M[3,11] := 1;
M[3,12] := 10;
M[3,13] := 2.6;
M[4,0] := 8;
M[4,1] := 3;
M[4,2] := 11;
M[4,3] := 12;
M[4,4] := 6;
M[4,5] := 0;
M[4,6] := 9;
M[4,7] := 5;
M[4,8] := 7;
M[4,9] := -7;
M[4,10] := 0;
M[4,11] := 1;
M[4,12] := 1.2;
M[4,13] := 2.9;
M[5,0] := 9;
M[5,1] := 5;
M[5,2] := 7;
M[5,3] := -7;
M[5,4] := 0;
M[5,5] := 1;
M[5,6] := 9;
M[5,7] := 5;
M[5,8] := 7;
M[5,9] := -7;
M[5,10] := 0;
M[5,11] := 1;
M[5,12] := 5;
M[5,13] := 8;
M[6,0] := 8;
M[6,1] := 0;
M[6,2] := -7;
M[6,3] := -3;
M[6,4] := 12;
M[6,5] := 10;
M[6,6] := 4;
M[6,7] := 5;
M[6,8] := 5;
M[6,9] := -2;
M[6,10] := 0;
M[6,11] := 7;
M[6,12] := 0;
M[6,13] := 2;
M[7,0] := 1;
M[7,1] := 1;
M[7,2] := 4;
M[7,3] := -0.1;
M[7,4] := 12;
M[7,5] := 1;
M[7,6] := 2;
M[7,7] := 6;
M[7,8] := 0.5;
M[7,9] := -4;
M[7,10] := 1;
M[7,11] := 0.5;
M[7,12] := 0;
M[7,13] := 2;
M[8,0] := 2;
M[8,1] := 7;
M[8,2] := 0.3;
M[8,3] := -7;
M[8,4] := 5;
M[8,5] := 1;
M[8,6] := 9;
M[8,7] := 1;
M[8,8] := 2;
M[8,9] := -6;
M[8,10] := 4.7;
M[8,11] := -2;
M[8,12] := 0;
M[8,13] := 2;
M[9,0] := 11.3;
M[9,1] := 14;
M[9,2] := 6;
M[9,3] := -7;
M[9,4] := 0;
M[9,5] := 1;
M[9,6] := 9;
M[9,7] := 5;
M[9,8] := 7;
M[9,9] := -7;
M[9,10] := 0;
M[9,11] := 1;
M[9,12] := 0;
M[9,13] := 2;
M[10,0] := 2;
M[10,1] := 3;
M[10,2] := -1;
M[10,3] := 6;
M[10,4] := 3;
M[10,5] := 5;
M[10,6] := 4;
M[10,7] := 3;
M[10,8] := 7;
M[10,9] := -1;
M[10,10] := 6;
M[10,11] := -1;
M[10,12] := 0;
M[10,13] := 2;
M[11,0] := -2;
M[11,1] := -13;
M[11,2] := -3;
M[11,3] := 7;
M[11,4] := 1.5;
M[11,5] := 15;
M[11,6] := 3;
M[11,7] := 3;
M[11,8] := 1;
M[11,9] := -5;
M[11,10] := 13;
M[11,11] := 7;
M[11,12] := 0;
M[11,13] := 2;
M[12,0] := -2;
M[12,1] := 2.3;
M[12,2] := -0.3;
M[12,3] := 0.7;
M[12,4] := 1.5;
M[12,5] := 5;
M[12,6] := 12;
M[12,7] := 1;
M[12,8] := 8;
M[12,9] := -4;
M[12,10] := 0;
M[12,11] := 3;
M[12,12] := 0;
M[12,13] := 2;
M[13,0] := -2;
M[13,1] := -13;
M[13,2] := -3;
M[13,3] := 7;
M[13,4] := 15;
M[13,5] := 1.5;
M[13,6] := 2;
M[13,7] := 1.1;
M[13,8] := 43;
M[13,9] := 10;
M[13,10] := 0;
M[13,11] := 2;
M[13,12] := 0;
M[13,13] := 2;
}

SetLength(M, 6, 6);
M[0,0] := 4;
M[0,1] := 2;
M[0,2] := 0;
M[0,3] := 1;
M[0,4] := 8;
M[0,5] := 9;
M[1,0] := 2;
M[1,1] := 3;
M[1,2] := -1;
M[1,3] := 6;
M[1,4] := 3;
M[1,5] := 5;
M[2,0] := 0;
M[2,1] := -1;
M[2,2] := 9;
M[2,3] := 0;
M[2,4] := 11;
M[2,5] := 7;
M[3,0] := 1;
M[3,1] := 6;
M[3,2] := 0;
M[3,3] := 1;
M[3,4] := 12;
M[3,5] := -7;
M[4,0] := 8;
M[4,1] := 3;
M[4,2] := 11;
M[4,3] := 12;
M[4,4] := 6;
M[4,5] := 0;
M[5,0] := 9;
M[5,1] := 5;
M[5,2] := 7;
M[5,3] := -7;
M[5,4] := 0;
M[5,5] := 1;

WriteLn('InputMatrix');
//MatrixPrint(M);
ReadLn();

{Am := MatrixIdentity(High(M)+1);
Am := MatrixMultiplyConst(Am,8.964010707377847);
WriteLn('Eigenmatrix');
MatrixPrint(Am);

Am := MatrixSubtract(M,Am);
WriteLn('Input - Eigenmatrix');
MatrixPrint(Am);

WriteLn('LU Decomp');
M := MatrixLUDecomp(M);
MatrixPrint(M);
ReadLn();}
MatrixPrint(M);
M := SingularValues(M);
//WriteLn('Singular Values');
MatrixPrint(M);
ReadLn();


SetLength(M, 6, 6);
M[0,0] := 4;
M[0,1] := 2;
M[0,2] := 0;
M[0,3] := 1;
M[0,4] := 8;
M[0,5] := 9;
M[1,0] := 2;
M[1,1] := 3;
M[1,2] := -1;
M[1,3] := 6;
M[1,4] := 3;
M[1,5] := 5;
M[2,0] := 0;
M[2,1] := -1;
M[2,2] := 9;
M[2,3] := 0;
M[2,4] := 11;
M[2,5] := 7;
M[3,0] := 1;
M[3,1] := 6;
M[3,2] := 0;
M[3,3] := 1;
M[3,4] := 12;
M[3,5] := -7;
M[4,0] := 8;
M[4,1] := 3;
M[4,2] := 11;
M[4,3] := 12;
M[4,4] := 6;
M[4,5] := 0;
M[5,0] := 9;
M[5,1] := 5;
M[5,2] := 7;
M[5,3] := -7;
M[5,4] := 0;
M[5,5] := 1;

{Setlength(M,2,2);
M[0,0] := 5;
M[0,1] := -2;
M[1,0] := -2;
M[1,1] := 8;}
Setlength(Am, 2,1);
Am[0,0] := 1;
Am[1,0] := 1;

//Bm := MatrixIdentity(High(M)+1);
//Bm := MatrixMultiplyConst(Am,9);
//WriteLn('Eigenmatrix');
//MatrixPrint(Am);
//M := MatrixSubtract(M,Bm);

MatrixPrint(M);

Eigen := Singularvalues(M);

writeln('eigenvalues');
matrixprint(Eigen);

count := 0;
I := 0;
while (I < 20) do
begin
     Setlength(Am, 6,1);
Am[0,0] := 1;
Am[1,0] := 1;
Am[2,0] := 1;
Am[3,0] := 1;
Am[4,0] := 1;
Am[5,0] := 1;

     while ( I < 200) do
     begin
          Am := MatrixMultiply(M, Am);
          Bm := Am;
          AnswerDouble := 1/Am[high(Am),0];
          //writeln('answerdouble',answerdouble);
          Am := MatrixMultiplyConst(Am, AnswerDouble);
          count := count +1;
     end;
     writeln('scaled');
     MatrixPrint(Bm);
     writeln('multiple');
     MatrixPrint(Am);
     count := 0;
     AnswerDouble := 1/MatrixNorm(Bm);
     writeln('Norm',1/AnswerDouble);
     Bm := MatrixMultiplyConst(Bm,AnswerDouble);
     Am := MatrixTranspose(Am);
     Am := MatrixMultiply(Bm,Am);
     Am := MatrixMultiplyConst(Am,25.948283413374774);
     M := MatrixSubtract(M,Am);
     writeln('vector');
     matrixprint(bm);
     readln();
end;

Am := MatrixSubtract(M,Am);
WriteLn('Input - Eigenmatrix');
MatrixPrint(Am);

end.

