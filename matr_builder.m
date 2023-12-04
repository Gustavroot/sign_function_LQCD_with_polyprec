clc;clear;

M = csvread("matr_data.txt");
Dreal = zeros(3072);
Dimag = zeros(3072);

M(:,1) = M(:,1) + 1;
M(:,2) = M(:,2) + 1;

for i=1:size(M,1)

  Dreal(int32(M(i,1)),int32(M(i,2))) = M(i,3);
  Dimag(int32(M(i,1)),int32(M(i,2))) = M(i,4);

end

D = Dreal + 1i*Dimag;
