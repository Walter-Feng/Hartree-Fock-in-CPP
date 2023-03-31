namespace hfincpp::integral::faster_kernel::rys_quadrature {
double electron_repulsive_integral_kernel_0003(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;

  double result  = (QD*QD*QD)*g0;
         result += (3*QD*QD)*g1;
         result += (3*QD)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_0012(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;

  double result  = (QC*QD*QD)*g0;
         result += (2*QC*QD+QD*QD)*g1;
         result += (QC+2*QD)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_0021(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;

  double result  = (QC*QC*QD)*g0;
         result += (QC*QC+2*QC*QD)*g1;
         result += (2*QC+QD)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_0030(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;
  double g3=2*B01*g1+D00*g2;

  double result  = (QC*QC*QC)*g0;
         result += (3*QC*QC)*g1;
         result += (3*QC)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_0102(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;

  double result  = (PB*QD*QD)*g0;
         result += (QD*QD)*g1;
         result += (2*PB*QD)*g2;
         result += (2*QD)*g3;
         result += (PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_0111(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;

  double result  = (PB*QC*QD)*g0;
         result += (QC*QD)*g1;
         result += (PB*QC+PB*QD)*g2;
         result += (QC+QD)*g3;
         result += (PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_0120(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;

  double result  = (PB*QC*QC)*g0;
         result += (QC*QC)*g1;
         result += (2*PB*QC)*g2;
         result += (2*QC)*g3;
         result += (PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_0201(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;

  double result  = (PB*PB*QD)*g0;
         result += (2*PB*QD)*g1;
         result += (QD)*g2;
         result += (PB*PB)*g3;
         result += (2*PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_0210(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;

  double result  = (PB*PB*QC)*g0;
         result += (2*PB*QC)*g1;
         result += (QC)*g2;
         result += (PB*PB)*g3;
         result += (2*PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_0300(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;

  double result  = (PB*PB*PB)*g0;
         result += (3*PB*PB)*g1;
         result += (3*PB)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_1002(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;

  double result  = (PA*QD*QD)*g0;
         result += (QD*QD)*g1;
         result += (2*PA*QD)*g2;
         result += (2*QD)*g3;
         result += (PA)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_1011(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;

  double result  = (PA*QC*QD)*g0;
         result += (QC*QD)*g1;
         result += (PA*QC+PA*QD)*g2;
         result += (QC+QD)*g3;
         result += (PA)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_1020(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;
  double g4=1*B01*g0+D00*g2;
  double g5=1*B00*g2+1*B01*g1+D00*g3;

  double result  = (PA*QC*QC)*g0;
         result += (QC*QC)*g1;
         result += (2*PA*QC)*g2;
         result += (2*QC)*g3;
         result += (PA)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_1101(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;

  double result  = (PA*PB*QD)*g0;
         result += (PA*QD+PB*QD)*g1;
         result += (QD)*g2;
         result += (PA*PB)*g3;
         result += (PA+PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_1110(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;

  double result  = (PA*PB*QC)*g0;
         result += (PA*QC+PB*QC)*g1;
         result += (QC)*g2;
         result += (PA*PB)*g3;
         result += (PA+PB)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_1200(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;

  double result  = (PA*PB*PB)*g0;
         result += (2*PA*PB+PB*PB)*g1;
         result += (PA+2*PB)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_2001(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;

  double result  = (PA*PA*QD)*g0;
         result += (2*PA*QD)*g1;
         result += (QD)*g2;
         result += (PA*PA)*g3;
         result += (2*PA)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_2010(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=D00*g0;
  double g4=1*B00*g0+D00*g1;
  double g5=2*B00*g1+D00*g2;

  double result  = (PA*PA*QC)*g0;
         result += (2*PA*QC)*g1;
         result += (QC)*g2;
         result += (PA*PA)*g3;
         result += (2*PA)*g4;
         result += (1)*g5;
  return result;
}
double electron_repulsive_integral_kernel_2100(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;

  double result  = (PA*PA*PB)*g0;
         result += (PA*PA+2*PA*PB)*g1;
         result += (2*PA+PB)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_3000(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;
  double g3=2*B10*g1+C00*g2;

  double result  = (PA*PA*PA)*g0;
         result += (3*PA*PA)*g1;
         result += (3*PA)*g2;
         result += (1)*g3;
  return result;
}
}