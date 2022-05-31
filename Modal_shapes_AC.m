function eq=Modal_shapes_AC(X)
global matrix
A=X(1);
C=X(2);

eq=[A*matrix(1,1)+C*matrix(1,2); A*matrix(2,1)+C*matrix(2,2)];