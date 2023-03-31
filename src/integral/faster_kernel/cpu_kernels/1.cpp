namespace hfincpp::integral::faster_kernel::rys_quadrature {
double electron_repulsive_integral_kernel_0001(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;

  double result  = (QD)*g0;
         result += (1)*g1;
  return result;
}
double electron_repulsive_integral_kernel_0010(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=D00*g0;

  double result  = (QC)*g0;
         result += (1)*g1;
  return result;
}
double electron_repulsive_integral_kernel_0100(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;

  double result  = (PB)*g0;
         result += (1)*g1;
  return result;
}
double electron_repulsive_integral_kernel_1000(const double C00, const double D00, const double B01, const double B10, const double B00, const double PA, const double PB, const double QC, const double QD, const double g0){
  double g1=C00*g0;

  double result  = (PA)*g0;
         result += (1)*g1;
  return result;
}
}