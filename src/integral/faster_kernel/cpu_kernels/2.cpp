namespace hfincpp::integral::faster_kernel::rys_quadrature {
double electron_repulsive_integral_kernel_0002(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;

  double result  = (QD*QD)*g0;
         result += (2*QD)*g1;
         result += (1)*g2;
  return result;
}
double electron_repulsive_integral_kernel_0011(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;

  double result  = (QC*QD)*g0;
         result += (QC+QD)*g1;
         result += (1)*g2;
  return result;
}
double electron_repulsive_integral_kernel_0020(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;
  double g2=1*B01*g0+D00*g1;

  double result  = (QC*QC)*g0;
         result += (2*QC)*g1;
         result += (1)*g2;
  return result;
}
double electron_repulsive_integral_kernel_0101(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;

  double result  = (PB*QD)*g0;
         result += (QD)*g1;
         result += (PB)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_0110(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;

  double result  = (PB*QC)*g0;
         result += (QC)*g1;
         result += (PB)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_0200(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;

  double result  = (PB*PB)*g0;
         result += (2*PB)*g1;
         result += (1)*g2;
  return result;
}
double electron_repulsive_integral_kernel_1001(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;

  double result  = (PA*QD)*g0;
         result += (QD)*g1;
         result += (PA)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_1010(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=D00*g0;
  double g3=1*B00*g0+D00*g1;

  double result  = (PA*QC)*g0;
         result += (QC)*g1;
         result += (PA)*g2;
         result += (1)*g3;
  return result;
}
double electron_repulsive_integral_kernel_1100(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;

  double result  = (PA*PB)*g0;
         result += (PA+PB)*g1;
         result += (1)*g2;
  return result;
}
double electron_repulsive_integral_kernel_2000(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;
  double g2=1*B10*g0+C00*g1;

  double result  = (PA*PA)*g0;
         result += (2*PA)*g1;
         result += (1)*g2;
  return result;
}
}