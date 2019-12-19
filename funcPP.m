function [val] = funcPP(z)
z1 = find((z>=50)&(z<=77));
z2 = find((z>77)&(z<=96));
z3 = find((z>96)&(z<=118));
z4 = find((z>118)&(z<=126));
z5 = find((z>126)&(z<=153));
z6 = find((z>153)&(z<=170));
z7 = find((z>170)&(z<=200));
A = fit([50 77]',[25.4542 46.86]','poly1');
val = zeros(length(z(:)),1);
val(z1) = A(z(z1));
B = fit([77 96]',[46.86 44.85]','poly1');
val(z2) = B(z(z2));
C = fit([96 118]',[44.85 21.45]','poly1');
val(z3) = C(z(z3));
D = fit([118 126]',[21.45 20.78]','poly1');
val(z4) = D(z(z4));
E = fit([126 153]',[20.78 51.27]','poly1');
val(z5) = E(z(z5));
F = fit([153 170]',[51.27 8]','poly1');
val(z6) = F(z(z6));
val(z7) = 8;
end

