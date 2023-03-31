namespace hfincpp::integral::faster_kernel::rys_quadrature {
double electron_repulsive_integral_kernel_0023(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;
  double g4=3*B01*g2+D00*g3;
  double g5=4*B01*g3+D00*g4;

  double result  = (QC*QC*QD*QD*QD)*g0;
         result += (3*QC*QC*QD*QD+2*QC*QD*QD*QD)*g1;
         result += (3*QC*QC*QD+6*QC*QD*QD+QD*QD*QD)*g2;
         result += (QC*QC+6*QC*QD+3*QD*QD)*g3;
         result += (2*QC+3*QD)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_0032(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;
  double g4=3*B01*g2+D00*g3;
  double g5=4*B01*g3+D00*g4;

  double result  = (QC*QC*QC*QD*QD)*g0;
         result += (2*QC*QC*QC*QD+3*QC*QC*QD*QD)*g1;
         result += (QC*QC*QC+6*QC*QC*QD+3*QC*QD*QD)*g2;
         result += (3*QC*QC+6*QC*QD+QD*QD)*g3;
         result += (3*QC+2*QD)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_0113(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;
  double g8=3*B01*g4+D00*g6;
  double g9=1*B00*g6+3*B01*g5+D00*g7;

  double result  = (PB*QC*QD*QD*QD)*g0;
         result += (QC*QD*QD*QD)*g1;
         result += (3*PB*QC*QD*QD+PB*QD*QD*QD)*g2;
         result += (3*QC*QD*QD+QD*QD*QD)*g3;
         result += (3*PB*QC*QD+3*PB*QD*QD)*g4;
         result += (3*QC*QD+3*QD*QD)*g5;
         result += (PB*QC+3*PB*QD)*g6;
         result += (QC+3*QD)*g7;
         result += (PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_0122(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;
  double g8=3*B01*g4+D00*g6;
  double g9=1*B00*g6+3*B01*g5+D00*g7;

  double result  = (PB*QC*QC*QD*QD)*g0;
         result += (QC*QC*QD*QD)*g1;
         result += (2*PB*QC*QC*QD+2*PB*QC*QD*QD)*g2;
         result += (2*QC*QC*QD+2*QC*QD*QD)*g3;
         result += (PB*QC*QC+4*PB*QC*QD+PB*QD*QD)*g4;
         result += (QC*QC+4*QC*QD+QD*QD)*g5;
         result += (2*PB*QC+2*PB*QD)*g6;
         result += (2*QC+2*QD)*g7;
         result += (PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_0131(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;
  double g8=3*B01*g4+D00*g6;
  double g9=1*B00*g6+3*B01*g5+D00*g7;

  double result  = (PB*QC*QC*QC*QD)*g0;
         result += (QC*QC*QC*QD)*g1;
         result += (PB*QC*QC*QC+3*PB*QC*QC*QD)*g2;
         result += (QC*QC*QC+3*QC*QC*QD)*g3;
         result += (3*PB*QC*QC+3*PB*QC*QD)*g4;
         result += (3*QC*QC+3*QC*QD)*g5;
         result += (3*PB*QC+PB*QD)*g6;
         result += (3*QC+QD)*g7;
         result += (PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_0203(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PB*PB*QD*QD*QD)*g0;
         result += (2*PB*QD*QD*QD)*g1;
         result += (QD*QD*QD)*g2;
         result += (3*PB*PB*QD*QD)*g3;
         result += (6*PB*QD*QD)*g4;
         result += (3*QD*QD)*g5;
         result += (3*PB*PB*QD)*g6;
         result += (6*PB*QD)*g7;
         result += (3*QD)*g8;
         result += (PB*PB)*g9;
         result += (2*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_0212(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PB*PB*QC*QD*QD)*g0;
         result += (2*PB*QC*QD*QD)*g1;
         result += (QC*QD*QD)*g2;
         result += (2*PB*PB*QC*QD+PB*PB*QD*QD)*g3;
         result += (4*PB*QC*QD+2*PB*QD*QD)*g4;
         result += (2*QC*QD+QD*QD)*g5;
         result += (PB*PB*QC+2*PB*PB*QD)*g6;
         result += (2*PB*QC+4*PB*QD)*g7;
         result += (QC+2*QD)*g8;
         result += (PB*PB)*g9;
         result += (2*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_0221(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PB*PB*QC*QC*QD)*g0;
         result += (2*PB*QC*QC*QD)*g1;
         result += (QC*QC*QD)*g2;
         result += (PB*PB*QC*QC+2*PB*PB*QC*QD)*g3;
         result += (2*PB*QC*QC+4*PB*QC*QD)*g4;
         result += (QC*QC+2*QC*QD)*g5;
         result += (2*PB*PB*QC+PB*PB*QD)*g6;
         result += (4*PB*QC+2*PB*QD)*g7;
         result += (2*QC+QD)*g8;
         result += (PB*PB)*g9;
         result += (2*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_0230(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PB*PB*QC*QC*QC)*g0;
         result += (2*PB*QC*QC*QC)*g1;
         result += (QC*QC*QC)*g2;
         result += (3*PB*PB*QC*QC)*g3;
         result += (6*PB*QC*QC)*g4;
         result += (3*QC*QC)*g5;
         result += (3*PB*PB*QC)*g6;
         result += (6*PB*QC)*g7;
         result += (3*QC)*g8;
         result += (PB*PB)*g9;
         result += (2*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_0302(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PB*PB*PB*QD*QD)*g0;
         result += (3*PB*PB*QD*QD)*g1;
         result += (3*PB*QD*QD)*g2;
         result += (QD*QD)*g3;
         result += (2*PB*PB*PB*QD)*g4;
         result += (6*PB*PB*QD)*g5;
         result += (6*PB*QD)*g6;
         result += (2*QD)*g7;
         result += (PB*PB*PB)*g8;
         result += (3*PB*PB)*g9;
         result += (3*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_0311(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PB*PB*PB*QC*QD)*g0;
         result += (3*PB*PB*QC*QD)*g1;
         result += (3*PB*QC*QD)*g2;
         result += (QC*QD)*g3;
         result += (PB*PB*PB*QC+PB*PB*PB*QD)*g4;
         result += (3*PB*PB*QC+3*PB*PB*QD)*g5;
         result += (3*PB*QC+3*PB*QD)*g6;
         result += (QC+QD)*g7;
         result += (PB*PB*PB)*g8;
         result += (3*PB*PB)*g9;
         result += (3*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_0320(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PB*PB*PB*QC*QC)*g0;
         result += (3*PB*PB*QC*QC)*g1;
         result += (3*PB*QC*QC)*g2;
         result += (QC*QC)*g3;
         result += (2*PB*PB*PB*QC)*g4;
         result += (6*PB*PB*QC)*g5;
         result += (6*PB*QC)*g6;
         result += (2*QC)*g7;
         result += (PB*PB*PB)*g8;
         result += (3*PB*PB)*g9;
         result += (3*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1013(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;
  double g8=3*B01*g4+D00*g6;
  double g9=1*B00*g6+3*B01*g5+D00*g7;

  double result  = (PA*QC*QD*QD*QD)*g0;
         result += (QC*QD*QD*QD)*g1;
         result += (3*PA*QC*QD*QD+PA*QD*QD*QD)*g2;
         result += (3*QC*QD*QD+QD*QD*QD)*g3;
         result += (3*PA*QC*QD+3*PA*QD*QD)*g4;
         result += (3*QC*QD+3*QD*QD)*g5;
         result += (PA*QC+3*PA*QD)*g6;
         result += (QC+3*QD)*g7;
         result += (PA)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_1022(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;
  double g8=3*B01*g4+D00*g6;
  double g9=1*B00*g6+3*B01*g5+D00*g7;

  double result  = (PA*QC*QC*QD*QD)*g0;
         result += (QC*QC*QD*QD)*g1;
         result += (2*PA*QC*QC*QD+2*PA*QC*QD*QD)*g2;
         result += (2*QC*QC*QD+2*QC*QD*QD)*g3;
         result += (PA*QC*QC+4*PA*QC*QD+PA*QD*QD)*g4;
         result += (QC*QC+4*QC*QD+QD*QD)*g5;
         result += (2*PA*QC+2*PA*QD)*g6;
         result += (2*QC+2*QD)*g7;
         result += (PA)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_1031(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;
  double g8=3*B01*g4+D00*g6;
  double g9=1*B00*g6+3*B01*g5+D00*g7;

  double result  = (PA*QC*QC*QC*QD)*g0;
         result += (QC*QC*QC*QD)*g1;
         result += (PA*QC*QC*QC+3*PA*QC*QC*QD)*g2;
         result += (QC*QC*QC+3*QC*QC*QD)*g3;
         result += (3*PA*QC*QC+3*PA*QC*QD)*g4;
         result += (3*QC*QC+3*QC*QD)*g5;
         result += (3*PA*QC+PA*QD)*g6;
         result += (3*QC+QD)*g7;
         result += (PA)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_1103(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PB*QD*QD*QD)*g0;
         result += (PA*QD*QD*QD+PB*QD*QD*QD)*g1;
         result += (QD*QD*QD)*g2;
         result += (3*PA*PB*QD*QD)*g3;
         result += (3*PA*QD*QD+3*PB*QD*QD)*g4;
         result += (3*QD*QD)*g5;
         result += (3*PA*PB*QD)*g6;
         result += (3*PA*QD+3*PB*QD)*g7;
         result += (3*QD)*g8;
         result += (PA*PB)*g9;
         result += (PA+PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1112(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PB*QC*QD*QD)*g0;
         result += (PA*QC*QD*QD+PB*QC*QD*QD)*g1;
         result += (QC*QD*QD)*g2;
         result += (2*PA*PB*QC*QD+PA*PB*QD*QD)*g3;
         result += (2*PA*QC*QD+2*PB*QC*QD+PA*QD*QD+PB*QD*QD)*g4;
         result += (2*QC*QD+QD*QD)*g5;
         result += (PA*PB*QC+2*PA*PB*QD)*g6;
         result += (PA*QC+PB*QC+2*PA*QD+2*PB*QD)*g7;
         result += (QC+2*QD)*g8;
         result += (PA*PB)*g9;
         result += (PA+PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1121(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PB*QC*QC*QD)*g0;
         result += (PA*QC*QC*QD+PB*QC*QC*QD)*g1;
         result += (QC*QC*QD)*g2;
         result += (PA*PB*QC*QC+2*PA*PB*QC*QD)*g3;
         result += (PA*QC*QC+PB*QC*QC+2*PA*QC*QD+2*PB*QC*QD)*g4;
         result += (QC*QC+2*QC*QD)*g5;
         result += (2*PA*PB*QC+PA*PB*QD)*g6;
         result += (2*PA*QC+2*PB*QC+PA*QD+PB*QD)*g7;
         result += (2*QC+QD)*g8;
         result += (PA*PB)*g9;
         result += (PA+PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1130(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PB*QC*QC*QC)*g0;
         result += (PA*QC*QC*QC+PB*QC*QC*QC)*g1;
         result += (QC*QC*QC)*g2;
         result += (3*PA*PB*QC*QC)*g3;
         result += (3*PA*QC*QC+3*PB*QC*QC)*g4;
         result += (3*QC*QC)*g5;
         result += (3*PA*PB*QC)*g6;
         result += (3*PA*QC+3*PB*QC)*g7;
         result += (3*QC)*g8;
         result += (PA*PB)*g9;
         result += (PA+PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1202(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PB*PB*QD*QD)*g0;
         result += (2*PA*PB*QD*QD+PB*PB*QD*QD)*g1;
         result += (PA*QD*QD+2*PB*QD*QD)*g2;
         result += (QD*QD)*g3;
         result += (2*PA*PB*PB*QD)*g4;
         result += (4*PA*PB*QD+2*PB*PB*QD)*g5;
         result += (2*PA*QD+4*PB*QD)*g6;
         result += (2*QD)*g7;
         result += (PA*PB*PB)*g8;
         result += (2*PA*PB+PB*PB)*g9;
         result += (PA+2*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1211(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PB*PB*QC*QD)*g0;
         result += (2*PA*PB*QC*QD+PB*PB*QC*QD)*g1;
         result += (PA*QC*QD+2*PB*QC*QD)*g2;
         result += (QC*QD)*g3;
         result += (PA*PB*PB*QC+PA*PB*PB*QD)*g4;
         result += (2*PA*PB*QC+PB*PB*QC+2*PA*PB*QD+PB*PB*QD)*g5;
         result += (PA*QC+2*PB*QC+PA*QD+2*PB*QD)*g6;
         result += (QC+QD)*g7;
         result += (PA*PB*PB)*g8;
         result += (2*PA*PB+PB*PB)*g9;
         result += (PA+2*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1220(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PB*PB*QC*QC)*g0;
         result += (2*PA*PB*QC*QC+PB*PB*QC*QC)*g1;
         result += (PA*QC*QC+2*PB*QC*QC)*g2;
         result += (QC*QC)*g3;
         result += (2*PA*PB*PB*QC)*g4;
         result += (4*PA*PB*QC+2*PB*PB*QC)*g5;
         result += (2*PA*QC+4*PB*QC)*g6;
         result += (2*QC)*g7;
         result += (PA*PB*PB)*g8;
         result += (2*PA*PB+PB*PB)*g9;
         result += (PA+2*PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_1301(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=D00*g0;
  double g6=1*B00*g0+D00*g1;
  double g7=2*B00*g1+D00*g2;
  double g8=3*B00*g2+D00*g3;
  double g9=4*B00*g3+D00*g4;

  double result  = (PA*PB*PB*PB*QD)*g0;
         result += (3*PA*PB*PB*QD+PB*PB*PB*QD)*g1;
         result += (3*PA*PB*QD+3*PB*PB*QD)*g2;
         result += (PA*QD+3*PB*QD)*g3;
         result += (QD)*g4;
         result += (PA*PB*PB*PB)*g5;
         result += (3*PA*PB*PB+PB*PB*PB)*g6;
         result += (3*PA*PB+3*PB*PB)*g7;
         result += (PA+3*PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_1310(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=D00*g0;
  double g6=1*B00*g0+D00*g1;
  double g7=2*B00*g1+D00*g2;
  double g8=3*B00*g2+D00*g3;
  double g9=4*B00*g3+D00*g4;

  double result  = (PA*PB*PB*PB*QC)*g0;
         result += (3*PA*PB*PB*QC+PB*PB*PB*QC)*g1;
         result += (3*PA*PB*QC+3*PB*PB*QC)*g2;
         result += (PA*QC+3*PB*QC)*g3;
         result += (QC)*g4;
         result += (PA*PB*PB*PB)*g5;
         result += (3*PA*PB*PB+PB*PB*PB)*g6;
         result += (3*PA*PB+3*PB*PB)*g7;
         result += (PA+3*PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_2003(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PA*QD*QD*QD)*g0;
         result += (2*PA*QD*QD*QD)*g1;
         result += (QD*QD*QD)*g2;
         result += (3*PA*PA*QD*QD)*g3;
         result += (6*PA*QD*QD)*g4;
         result += (3*QD*QD)*g5;
         result += (3*PA*PA*QD)*g6;
         result += (6*PA*QD)*g7;
         result += (3*QD)*g8;
         result += (PA*PA)*g9;
         result += (2*PA)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_2012(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PA*QC*QD*QD)*g0;
         result += (2*PA*QC*QD*QD)*g1;
         result += (QC*QD*QD)*g2;
         result += (2*PA*PA*QC*QD+PA*PA*QD*QD)*g3;
         result += (4*PA*QC*QD+2*PA*QD*QD)*g4;
         result += (2*QC*QD+QD*QD)*g5;
         result += (PA*PA*QC+2*PA*PA*QD)*g6;
         result += (2*PA*QC+4*PA*QD)*g7;
         result += (QC+2*QD)*g8;
         result += (PA*PA)*g9;
         result += (2*PA)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_2021(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PA*QC*QC*QD)*g0;
         result += (2*PA*QC*QC*QD)*g1;
         result += (QC*QC*QD)*g2;
         result += (PA*PA*QC*QC+2*PA*PA*QC*QD)*g3;
         result += (2*PA*QC*QC+4*PA*QC*QD)*g4;
         result += (QC*QC+2*QC*QD)*g5;
         result += (2*PA*PA*QC+PA*PA*QD)*g6;
         result += (4*PA*QC+2*PA*QD)*g7;
         result += (2*QC+QD)*g8;
         result += (PA*PA)*g9;
         result += (2*PA)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_2030(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;
  double g9=2*B01*g3+D00*g6;
  double g10=1*B00*g6+2*B01*g4+D00*g7;
  double g11=2*B00*g7+2*B01*g5+D00*g8;

  double result  = (PA*PA*QC*QC*QC)*g0;
         result += (2*PA*QC*QC*QC)*g1;
         result += (QC*QC*QC)*g2;
         result += (3*PA*PA*QC*QC)*g3;
         result += (6*PA*QC*QC)*g4;
         result += (3*QC*QC)*g5;
         result += (3*PA*PA*QC)*g6;
         result += (6*PA*QC)*g7;
         result += (3*QC)*g8;
         result += (PA*PA)*g9;
         result += (2*PA)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_2102(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PA*PB*QD*QD)*g0;
         result += (PA*PA*QD*QD+2*PA*PB*QD*QD)*g1;
         result += (2*PA*QD*QD+PB*QD*QD)*g2;
         result += (QD*QD)*g3;
         result += (2*PA*PA*PB*QD)*g4;
         result += (2*PA*PA*QD+4*PA*PB*QD)*g5;
         result += (4*PA*QD+2*PB*QD)*g6;
         result += (2*QD)*g7;
         result += (PA*PA*PB)*g8;
         result += (PA*PA+2*PA*PB)*g9;
         result += (2*PA+PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_2111(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PA*PB*QC*QD)*g0;
         result += (PA*PA*QC*QD+2*PA*PB*QC*QD)*g1;
         result += (2*PA*QC*QD+PB*QC*QD)*g2;
         result += (QC*QD)*g3;
         result += (PA*PA*PB*QC+PA*PA*PB*QD)*g4;
         result += (PA*PA*QC+2*PA*PB*QC+PA*PA*QD+2*PA*PB*QD)*g5;
         result += (2*PA*QC+PB*QC+2*PA*QD+PB*QD)*g6;
         result += (QC+QD)*g7;
         result += (PA*PA*PB)*g8;
         result += (PA*PA+2*PA*PB)*g9;
         result += (2*PA+PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_2120(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PA*PB*QC*QC)*g0;
         result += (PA*PA*QC*QC+2*PA*PB*QC*QC)*g1;
         result += (2*PA*QC*QC+PB*QC*QC)*g2;
         result += (QC*QC)*g3;
         result += (2*PA*PA*PB*QC)*g4;
         result += (2*PA*PA*QC+4*PA*PB*QC)*g5;
         result += (4*PA*QC+2*PB*QC)*g6;
         result += (2*QC)*g7;
         result += (PA*PA*PB)*g8;
         result += (PA*PA+2*PA*PB)*g9;
         result += (2*PA+PB)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_2201(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=D00*g0;
  double g6=1*B00*g0+D00*g1;
  double g7=2*B00*g1+D00*g2;
  double g8=3*B00*g2+D00*g3;
  double g9=4*B00*g3+D00*g4;

  double result  = (PA*PA*PB*PB*QD)*g0;
         result += (2*PA*PA*PB*QD+2*PA*PB*PB*QD)*g1;
         result += (PA*PA*QD+4*PA*PB*QD+PB*PB*QD)*g2;
         result += (2*PA*QD+2*PB*QD)*g3;
         result += (QD)*g4;
         result += (PA*PA*PB*PB)*g5;
         result += (2*PA*PA*PB+2*PA*PB*PB)*g6;
         result += (PA*PA+4*PA*PB+PB*PB)*g7;
         result += (2*PA+2*PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_2210(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=D00*g0;
  double g6=1*B00*g0+D00*g1;
  double g7=2*B00*g1+D00*g2;
  double g8=3*B00*g2+D00*g3;
  double g9=4*B00*g3+D00*g4;

  double result  = (PA*PA*PB*PB*QC)*g0;
         result += (2*PA*PA*PB*QC+2*PA*PB*PB*QC)*g1;
         result += (PA*PA*QC+4*PA*PB*QC+PB*PB*QC)*g2;
         result += (2*PA*QC+2*PB*QC)*g3;
         result += (QC)*g4;
         result += (PA*PA*PB*PB)*g5;
         result += (2*PA*PA*PB+2*PA*PB*PB)*g6;
         result += (PA*PA+4*PA*PB+PB*PB)*g7;
         result += (2*PA+2*PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_2300(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=4*B10*g3+C00*g4;

  double result  = (PA*PA*PB*PB*PB)*g0;
         result += (3*PA*PA*PB*PB+2*PA*PB*PB*PB)*g1;
         result += (3*PA*PA*PB+6*PA*PB*PB+PB*PB*PB)*g2;
         result += (PA*PA+6*PA*PB+3*PB*PB)*g3;
         result += (2*PA+3*PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_3002(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PA*PA*QD*QD)*g0;
         result += (3*PA*PA*QD*QD)*g1;
         result += (3*PA*QD*QD)*g2;
         result += (QD*QD)*g3;
         result += (2*PA*PA*PA*QD)*g4;
         result += (6*PA*PA*QD)*g5;
         result += (6*PA*QD)*g6;
         result += (2*QD)*g7;
         result += (PA*PA*PA)*g8;
         result += (3*PA*PA)*g9;
         result += (3*PA)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_3011(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PA*PA*QC*QD)*g0;
         result += (3*PA*PA*QC*QD)*g1;
         result += (3*PA*QC*QD)*g2;
         result += (QC*QD)*g3;
         result += (PA*PA*PA*QC+PA*PA*PA*QD)*g4;
         result += (3*PA*PA*QC+3*PA*PA*QD)*g5;
         result += (3*PA*QC+3*PA*QD)*g6;
         result += (QC+QD)*g7;
         result += (PA*PA*PA)*g8;
         result += (3*PA*PA)*g9;
         result += (3*PA)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_3020(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;
  double g8=1*B01*g0+D00*g4;
  double g9=1*B00*g4+1*B01*g1+D00*g5;
  double g10=2*B00*g5+1*B01*g2+D00*g6;
  double g11=3*B00*g6+1*B01*g3+D00*g7;

  double result  = (PA*PA*PA*QC*QC)*g0;
         result += (3*PA*PA*QC*QC)*g1;
         result += (3*PA*QC*QC)*g2;
         result += (QC*QC)*g3;
         result += (2*PA*PA*PA*QC)*g4;
         result += (6*PA*PA*QC)*g5;
         result += (6*PA*QC)*g6;
         result += (2*QC)*g7;
         result += (PA*PA*PA)*g8;
         result += (3*PA*PA)*g9;
         result += (3*PA)*g10;
         result += (1)*g11;
  return result;
}
double electron_repulsive_integral_kernel_3101(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=D00*g0;
  double g6=1*B00*g0+D00*g1;
  double g7=2*B00*g1+D00*g2;
  double g8=3*B00*g2+D00*g3;
  double g9=4*B00*g3+D00*g4;

  double result  = (PA*PA*PA*PB*QD)*g0;
         result += (PA*PA*PA*QD+3*PA*PA*PB*QD)*g1;
         result += (3*PA*PA*QD+3*PA*PB*QD)*g2;
         result += (3*PA*QD+PB*QD)*g3;
         result += (QD)*g4;
         result += (PA*PA*PA*PB)*g5;
         result += (PA*PA*PA+3*PA*PA*PB)*g6;
         result += (3*PA*PA+3*PA*PB)*g7;
         result += (3*PA+PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_3110(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=D00*g0;
  double g6=1*B00*g0+D00*g1;
  double g7=2*B00*g1+D00*g2;
  double g8=3*B00*g2+D00*g3;
  double g9=4*B00*g3+D00*g4;

  double result  = (PA*PA*PA*PB*QC)*g0;
         result += (PA*PA*PA*QC+3*PA*PA*PB*QC)*g1;
         result += (3*PA*PA*QC+3*PA*PB*QC)*g2;
         result += (3*PA*QC+PB*QC)*g3;
         result += (QC)*g4;
         result += (PA*PA*PA*PB)*g5;
         result += (PA*PA*PA+3*PA*PA*PB)*g6;
         result += (3*PA*PA+3*PA*PB)*g7;
         result += (3*PA+PB)*g8;
         result += (1)*g9;
  return result;
}
double electron_repulsive_integral_kernel_3200(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;
  double g5=4*B10*g3+C00*g4;

  double result  = (PA*PA*PA*PB*PB)*g0;
         result += (2*PA*PA*PA*PB+3*PA*PA*PB*PB)*g1;
         result += (PA*PA*PA+6*PA*PA*PB+3*PA*PB*PB)*g2;
         result += (3*PA*PA+6*PA*PB+PB*PB)*g3;
         result += (3*PA+2*PB)*g4;
         result += (1)*g5;
  return result;
}
}