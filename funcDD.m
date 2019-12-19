function [func2] = funcDD(z)
func2 = zeros(length(z(:)),1);
z1 = find((z>=-76)&(z<=-50));
z2 = find((z>=-94)&(z<-76));
z3 = find((z>=-117)&(z<-94));
z4 = find((z>=-126)&(z<-117));
z5 = find((z>=-156)&(z<-126));
z6 = find((z>=-191)&(z<-156));
z7 = find((z>=-200)&(z<-191));
func2(z1) = 47.81 - 0.8307 * (z(z1) + 76);
func2(z2) = 47.81;
func2(z3) = 47.81 + 1.12217 * (z(z3) + 94);
func2(z4) = 22;
func2(z5) = 22 - 1.66607 * (z(z5) + 126);
func2(z6) = 22 + 1.42343 * (z(z6) + 191);
func2(z7) = 22;
end

