namespace hfincpp::integral::faster_kernel::rys_quadrature {
double electron_repulsive_integral_kernel_0013(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;
  double g4=3*B01*g2+D00*g3;

  double result  = (QC*QD*QD*QD)*g0;
         result += (3*QC*QD*QD+QD*QD*QD)*g1;
         result += (3*QC*QD+3*QD*QD)*g2;
         result += (QC+3*QD)*g3;
         result += (1)*g4;
  return result;
}
double electron_repulsive_integral_kernel_0022(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;
  double g4=3*B01*g2+D00*g3;

  double result  = (QC*QC*QD*QD)*g0;
         result += (2*QC*QC*QD+2*QC*QD*QD)*g1;
         result += (QC*QC+4*QC*QD+QD*QD)*g2;
         result += (2*QC+2*QD)*g3;
         result += (1)*g4;
  return result;
}
double electron_repulsive_integral_kernel_0031(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;
  double g4=3*B01*g2+D00*g3;

  double result  = (QC*QC*QC*QD)*g0;
         result += (QC*QC*QC+3*QC*QC*QD)*g1;
         result += (3*QC*QC+3*QC*QD)*g2;
         result += (3*QC+QD)*g3;
         result += (1)*g4;
  return result;
}
double electron_repulsive_integral_kernel_0103(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PB*QD*QD*QD)*g0;
         result += (QD*QD*QD)*g1;
         result += (3*PB*QD*QD)*g2;
         result += (3*QD*QD)*g3;
         result += (3*PB*QD)*g4;
         result += (3*QD)*g5;
         result += (PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_0112(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PB*QC*QD*QD)*g0;
         result += (QC*QD*QD)*g1;
         result += (2*PB*QC*QD+PB*QD*QD)*g2;
         result += (2*QC*QD+QD*QD)*g3;
         result += (PB*QC+2*PB*QD)*g4;
         result += (QC+2*QD)*g5;
         result += (PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_0121(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PB*QC*QC*QD)*g0;
         result += (QC*QC*QD)*g1;
         result += (PB*QC*QC+2*PB*QC*QD)*g2;
         result += (QC*QC+2*QC*QD)*g3;
         result += (2*PB*QC+PB*QD)*g4;
         result += (2*QC+QD)*g5;
         result += (PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_0130(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PB*QC*QC*QC)*g0;
         result += (QC*QC*QC)*g1;
         result += (3*PB*QC*QC)*g2;
         result += (3*QC*QC)*g3;
         result += (3*PB*QC)*g4;
         result += (3*QC)*g5;
         result += (PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_0202(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PB*PB*QD*QD)*g0;
         result += (2*PB*QD*QD)*g1;
         result += (QD*QD)*g2;
         result += (2*PB*PB*QD)*g3;
         result += (4*PB*QD)*g4;
         result += (2*QD)*g5;
         result += (PB*PB)*g6;
         result += (2*PB)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_0211(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PB*PB*QC*QD)*g0;
         result += (2*PB*QC*QD)*g1;
         result += (QC*QD)*g2;
         result += (PB*PB*QC+PB*PB*QD)*g3;
         result += (2*PB*QC+2*PB*QD)*g4;
         result += (QC+QD)*g5;
         result += (PB*PB)*g6;
         result += (2*PB)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_0220(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PB*PB*QC*QC)*g0;
         result += (2*PB*QC*QC)*g1;
         result += (QC*QC)*g2;
         result += (2*PB*PB*QC)*g3;
         result += (4*PB*QC)*g4;
         result += (2*QC)*g5;
         result += (PB*PB)*g6;
         result += (2*PB)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_0301(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PB*PB*PB*QD)*g0;
         result += (3*PB*PB*QD)*g1;
         result += (3*PB*QD)*g2;
         result += (QD)*g3;
         result += (PB*PB*PB)*g4;
         result += (3*PB*PB)*g5;
         result += (3*PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_0310(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PB*PB*PB*QC)*g0;
         result += (3*PB*PB*QC)*g1;
         result += (3*PB*QC)*g2;
         result += (QC)*g3;
         result += (PB*PB*PB)*g4;
         result += (3*PB*PB)*g5;
         result += (3*PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_1003(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PA*QD*QD*QD)*g0;
         result += (QD*QD*QD)*g1;
         result += (3*PA*QD*QD)*g2;
         result += (3*QD*QD)*g3;
         result += (3*PA*QD)*g4;
         result += (3*QD)*g5;
         result += (PA)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_1012(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PA*QC*QD*QD)*g0;
         result += (QC*QD*QD)*g1;
         result += (2*PA*QC*QD+PA*QD*QD)*g2;
         result += (2*QC*QD+QD*QD)*g3;
         result += (PA*QC+2*PA*QD)*g4;
         result += (QC+2*QD)*g5;
         result += (PA)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_1021(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PA*QC*QC*QD)*g0;
         result += (QC*QC*QD)*g1;
         result += (PA*QC*QC+2*PA*QC*QD)*g2;
         result += (QC*QC+2*QC*QD)*g3;
         result += (2*PA*QC+PA*QD)*g4;
         result += (2*QC+QD)*g5;
         result += (PA)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_1030(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;
  double g6=2*B01*g2+D00*g4;
  double g7=1*B00*g4+2*B01*g3+D00*g5;

  double result  = (PA*QC*QC*QC)*g0;
         result += (QC*QC*QC)*g1;
         result += (3*PA*QC*QC)*g2;
         result += (3*QC*QC)*g3;
         result += (3*PA*QC)*g4;
         result += (3*QC)*g5;
         result += (PA)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_1102(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PA*PB*QD*QD)*g0;
         result += (PA*QD*QD+PB*QD*QD)*g1;
         result += (QD*QD)*g2;
         result += (2*PA*PB*QD)*g3;
         result += (2*PA*QD+2*PB*QD)*g4;
         result += (2*QD)*g5;
         result += (PA*PB)*g6;
         result += (PA+PB)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_1111(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PA*PB*QC*QD)*g0;
         result += (PA*QC*QD+PB*QC*QD)*g1;
         result += (QC*QD)*g2;
         result += (PA*PB*QC+PA*PB*QD)*g3;
         result += (PA*QC+PB*QC+PA*QD+PB*QD)*g4;
         result += (QC+QD)*g5;
         result += (PA*PB)*g6;
         result += (PA+PB)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_1120(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PA*PB*QC*QC)*g0;
         result += (PA*QC*QC+PB*QC*QC)*g1;
         result += (QC*QC)*g2;
         result += (2*PA*PB*QC)*g3;
         result += (2*PA*QC+2*PB*QC)*g4;
         result += (2*QC)*g5;
         result += (PA*PB)*g6;
         result += (PA+PB)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_1201(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PA*PB*PB*QD)*g0;
         result += (2*PA*PB*QD+PB*PB*QD)*g1;
         result += (PA*QD+2*PB*QD)*g2;
         result += (QD)*g3;
         result += (PA*PB*PB)*g4;
         result += (2*PA*PB+PB*PB)*g5;
         result += (PA+2*PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_1210(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PA*PB*PB*QC)*g0;
         result += (2*PA*PB*QC+PB*PB*QC)*g1;
         result += (PA*QC+2*PB*QC)*g2;
         result += (QC)*g3;
         result += (PA*PB*PB)*g4;
         result += (2*PA*PB+PB*PB)*g5;
         result += (PA+2*PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_1300(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;

  double result  = (PA*PB*PB*PB)*g0;
         result += (3*PA*PB*PB+PB*PB*PB)*g1;
         result += (3*PA*PB+3*PB*PB)*g2;
         result += (PA+3*PB)*g3;
         result += (1)*g4;
  return result;
}
double electron_repulsive_integral_kernel_2002(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PA*PA*QD*QD)*g0;
         result += (2*PA*QD*QD)*g1;
         result += (QD*QD)*g2;
         result += (2*PA*PA*QD)*g3;
         result += (4*PA*QD)*g4;
         result += (2*QD)*g5;
         result += (PA*PA)*g6;
         result += (2*PA)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_2011(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PA*PA*QC*QD)*g0;
         result += (2*PA*QC*QD)*g1;
         result += (QC*QD)*g2;
         result += (PA*PA*QC+PA*PA*QD)*g3;
         result += (2*PA*QC+2*PA*QD)*g4;
         result += (QC+QD)*g5;
         result += (PA*PA)*g6;
         result += (2*PA)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_2020(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;
  double g6=1*B01*g0+D00*g3;
  double g7=1*B00*g3+1*B01*g1+D00*g4;
  double g8=2*B00*g4+1*B01*g2+D00*g5;

  double result  = (PA*PA*QC*QC)*g0;
         result += (2*PA*QC*QC)*g1;
         result += (QC*QC)*g2;
         result += (2*PA*PA*QC)*g3;
         result += (4*PA*QC)*g4;
         result += (2*QC)*g5;
         result += (PA*PA)*g6;
         result += (2*PA)*g7;
         result += (1)*g8;
  return result;
}
double electron_repulsive_integral_kernel_2101(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PA*PA*PB*QD)*g0;
         result += (PA*PA*QD+2*PA*PB*QD)*g1;
         result += (2*PA*QD+PB*QD)*g2;
         result += (QD)*g3;
         result += (PA*PA*PB)*g4;
         result += (PA*PA+2*PA*PB)*g5;
         result += (2*PA+PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_2110(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PA*PA*PB*QC)*g0;
         result += (PA*PA*QC+2*PA*PB*QC)*g1;
         result += (2*PA*QC+PB*QC)*g2;
         result += (QC)*g3;
         result += (PA*PA*PB)*g4;
         result += (PA*PA+2*PA*PB)*g5;
         result += (2*PA+PB)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_2200(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;

  double result  = (PA*PA*PB*PB)*g0;
         result += (2*PA*PA*PB+2*PA*PB*PB)*g1;
         result += (PA*PA+4*PA*PB+PB*PB)*g2;
         result += (2*PA+2*PB)*g3;
         result += (1)*g4;
  return result;
}
double electron_repulsive_integral_kernel_3001(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PA*PA*PA*QD)*g0;
         result += (3*PA*PA*QD)*g1;
         result += (3*PA*QD)*g2;
         result += (QD)*g3;
         result += (PA*PA*PA)*g4;
         result += (3*PA*PA)*g5;
         result += (3*PA)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_3010(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=D00*g0;
  double g5=1*B00*g0+D00*g1;
  double g6=2*B00*g1+D00*g2;
  double g7=3*B00*g2+D00*g3;

  double result  = (PA*PA*PA*QC)*g0;
         result += (3*PA*PA*QC)*g1;
         result += (3*PA*QC)*g2;
         result += (QC)*g3;
         result += (PA*PA*PA)*g4;
         result += (3*PA*PA)*g5;
         result += (3*PA)*g6;
         result += (1)*g7;
  return result;
}
double electron_repulsive_integral_kernel_3100(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;
  double g4=3*B10*g2+C00*g3;

  double result  = (PA*PA*PA*PB)*g0;
         result += (PA*PA*PA+3*PA*PA*PB)*g1;
         result += (3*PA*PA+3*PA*PB)*g2;
         result += (3*PA+PB)*g3;
         result += (1)*g4;
  return result;
}
}