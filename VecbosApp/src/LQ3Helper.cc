#include "LQ3Helper.hh"

FourVector::FourVector()
{
   P[0] = 0;
   P[1] = 1;
   P[2] = 2;
   P[3] = 3;
}

FourVector::FourVector(double p[4])
{
   P[0] = p[0];
   P[1] = p[1];
   P[2] = p[2];
   P[3] = p[3];
}

FourVector::FourVector(double e, double px, double py, double pz)
{
   P[0] = e;
   P[1] = px;
   P[2] = py;
   P[3] = pz;
}

FourVector::~FourVector()
{
}

void FourVector::SetPtEtaPhi(double pt, double eta, double phi)
{
   SetPtEtaPhiMass(pt, eta, phi, 0);
}

void FourVector::SetPtEtaPhiMass(double pt, double eta, double phi, double mass)
{
   P[1] = pt * cos(phi);
   P[2] = pt * sin(phi);
   P[3] = pt * sinh(eta);

   P[0] = sqrt(mass * mass + SpatialDot(*this));
}

void FourVector::SetSizeEtaPhi(double size, double eta, double phi)
{
   SetSizeEtaPhiMass(size, eta, phi, 0);
}

void FourVector::SetSizeEtaPhiMass(double size, double eta, double phi, double mass)
{
   P[0] = sqrt(size * size + mass * mass);
   P[1] = size / cosh(eta) * cos(phi);
   P[2] = size / cosh(eta) * sin(phi);
   P[3] = size * tanh(eta);
}

void FourVector::SetSizeEtaPhiEnergy(double size, double eta, double phi, double energy)
{
   P[0] = energy;
   P[1] = size / cosh(eta) * cos(phi);
   P[2] = size / cosh(eta) * sin(phi);
   P[3] = size * tanh(eta);
}

void FourVector::SetSizeThetaPhi(double size, double theta, double phi)
{
   SetSizeThetaPhiMass(size, theta, phi, 0);
}

void FourVector::SetSizeThetaPhiMass(double size, double theta, double phi, double mass)
{
   P[0] = sqrt(size * size + mass * mass);
   P[1] = size * sin(theta) * cos(phi);
   P[2] = size * sin(theta) * sin(phi);
   P[3] = size * cos(theta);
}

double &FourVector::operator [](int index)
{
   if(index >= 0 && index <= 3)
      return P[index];
   return P[0];
}

double FourVector::operator [](int index) const
{
   if(index >= 0 && index <= 3)
      return P[index];
   return 0;
}

FourVector &FourVector::operator =(const FourVector &Other)
{
   P[0] = Other.P[0];
   P[1] = Other.P[1];
   P[2] = Other.P[2];
   P[3] = Other.P[3];

   return *this;
}
   
FourVector FourVector::operator +(const FourVector &Other) const
{
   FourVector Out;
   Out.P[0] = P[0] + Other.P[0];
   Out.P[1] = P[1] + Other.P[1];
   Out.P[2] = P[2] + Other.P[2];
   Out.P[3] = P[3] + Other.P[3];
   return Out;
}

FourVector FourVector::operator -() const
{
   FourVector Out;
   Out.P[0] = -P[0];
   Out.P[1] = -P[1];
   Out.P[2] = -P[2];
   Out.P[3] = -P[3];
   return Out;
}

FourVector FourVector::operator -(const FourVector &Other) const
{
   FourVector Out;
   Out.P[0] = P[0] - Other.P[0];
   Out.P[1] = P[1] - Other.P[1];
   Out.P[2] = P[2] - Other.P[2];
   Out.P[3] = P[3] - Other.P[3];
   return Out;
}

FourVector FourVector::operator *(double Scale) const
{
   FourVector Out;
   Out.P[0] = P[0] * Scale;
   Out.P[1] = P[1] * Scale;
   Out.P[2] = P[2] * Scale;
   Out.P[3] = P[3] * Scale;
   return Out;
}

double FourVector::GetMass() const
{
   double Mass2 = GetMass2();
   
   if(Mass2 >= 0)
      return sqrt(Mass2);
   return 0;
}

double FourVector::GetMass2() const
{
   return MetricDot(*this);
}

double FourVector::GetP() const
{
   return sqrt(SpatialDot(*this));
}

double FourVector::GetPT() const
{
   return sqrt(P[1] * P[1] + P[2] * P[2]);
}

double FourVector::GetEta() const
{
   double Momentum = GetP();

   return 0.5 * log((Momentum + P[3]) / (Momentum - P[3]));
}

double FourVector::GetRapidity() const
{
   return 0.5 * log((P[0] + P[3]) / (P[0] - P[3]));
}

double FourVector::GetY() const
{
   return GetRapidity();
}

double FourVector::GetPhi() const
{
   double PT = GetPT();

   double Angle = acos(P[1] / PT);
   if(P[2] < 0)
      Angle = -Angle;

   return Angle;
}

double FourVector::GetTheta() const
{
   return acos(P[3] / GetP());
}

double FourVector::GetBeta() const
{
   double Gamma = GetGamma();
   return sqrt(1 - 1 / (Gamma * Gamma));
}

double FourVector::GetGamma() const
{
   return P[0] / GetMass();
}

FourVector FourVector::RotateX(double Angle) const
{
   FourVector Out;
   Out.P[0] = P[0];
   Out.P[1] = P[1];
   Out.P[2] = cos(Angle) * P[2] - sin(Angle) * P[3];
   Out.P[3] = sin(Angle) * P[2] + cos(Angle) * P[3];
   return Out;
}

FourVector FourVector::RotateY(double Angle) const
{
   FourVector Out;
   Out.P[0] = P[0];
   Out.P[1] = sin(Angle) * P[3] + cos(Angle) * P[1];
   Out.P[2] = P[2];
   Out.P[3] = cos(Angle) * P[3] - sin(Angle) * P[1];
   return Out;
}

FourVector FourVector::RotateZ(double Angle) const
{
   FourVector Out;
   Out.P[0] = P[0];
   Out.P[1] = cos(Angle) * P[1] - sin(Angle) * P[2];
   Out.P[2] = sin(Angle) * P[1] + cos(Angle) * P[2];
   Out.P[3] = P[3];
   return Out;
}

FourVector FourVector::Rotate(const FourVector Axis, double Angle) const
{
   // rotate "axis" and input to y-z plane, then rotate "axis" to z axis,
   //    rotate input with respect to z axis
   //    and then rotate back

   double Psi = PI / 2 - Axis.GetPhi();
   double Theta = acos(Axis[3] / Axis.GetP());

   return RotateZ(Psi).RotateX(Theta).RotateZ(Angle).RotateX(-Theta).RotateZ(-Psi);
}

FourVector FourVector::BoostX(double Beta) const
{
   double Gamma = BetaToGamma(Beta);

   FourVector Out;
   Out.P[0] = Gamma * P[0] - Beta * Gamma * P[1];
   Out.P[1] = -Beta * Gamma * P[0] + Gamma * P[1];
   Out.P[2] = P[2];
   Out.P[3] = P[3];
   return Out;
}

FourVector FourVector::BoostY(double Beta) const
{
   double Gamma = BetaToGamma(Beta);

   FourVector Out;
   Out.P[0] = Gamma * P[0] - Beta * Gamma * P[2];
   Out.P[1] = P[1];
   Out.P[2] = -Beta * Gamma * P[0] + Gamma * P[2];
   Out.P[3] = P[3];
   return Out;
}

FourVector FourVector::BoostZ(double Beta) const
{
   double Gamma = BetaToGamma(Beta);

   FourVector Out;
   Out.P[0] = Gamma * P[0] - Beta * Gamma * P[3];
   Out.P[1] = P[1];
   Out.P[2] = P[2];
   Out.P[3] = -Beta * Gamma * P[0] + Gamma * P[3];
   return Out;
}

FourVector FourVector::Boost(const FourVector Axis, double Beta) const
{
   if(Axis.GetPT() < 1e-8)   // axis along z direction
   {
      if(Axis[3] > 0)
         return BoostZ(Beta);
      else
         return BoostZ(-Beta);
   }

   double Psi = PI / 2 - Axis.GetPhi();
   double Theta = acos(Axis[3] / Axis.GetP());

   return RotateZ(Psi).RotateX(Theta).BoostZ(Beta).RotateX(-Theta).RotateZ(-Psi);
}
   
FourVector FourVector::SmearAngle(double Angle) const
{
   /*
   FourVector Reference(0, 1, 0, 0);
   if(fabs(P[2]) < 1e-6 && fabs(P[3]) < 1e-6)
      Reference[2] = 1;

   FourVector Axis = SpatialCross(Reference);   // so that axis is perpendicular to input momentum

   FourVector RealAxis;   // pick a random rotation axis perpendicular to input momentum
   double AxisRotation = DrawRandom(0, 2 * PI);
   RealAxis = Axis.Rotate(*this, AxisRotation);

   double SmearAngle = DrawGaussian(Angle);
   return Rotate(RealAxis, SmearAngle);
   */

   return *this;
}

FourVector FourVector::SmearMomentum(double Scale) const
{
   // double Factor = 1 + DrawGaussian(Scale);
   double Factor = 1;

   return (*this) * Factor;
}

FourVector FourVector::SpatialCross(const FourVector Other) const
{
   FourVector Out;
   Out.P[0] = P[0];
   Out.P[1] = P[2] * Other.P[3] - P[3] * Other.P[2];
   Out.P[2] = P[3] * Other.P[1] - P[1] * Other.P[3];
   Out.P[3] = P[1] * Other.P[2] - P[2] * Other.P[1];
}

double FourVector::SpatialDot(const FourVector &Other) const
{
   return P[1] * Other.P[1] + P[2] * Other.P[2] + P[3] * Other.P[3];
}

double FourVector::MetricDot(const FourVector &Other) const
{
   return P[0] * Other.P[0] - SpatialDot(Other);
}

double GetAngle(const FourVector P1, const FourVector P2)
{
   return acos(P1.SpatialDot(P2) / P1.GetP() / P2.GetP());
}

double GetDR(const FourVector P1, const FourVector P2)
{
   double DEta = P1.GetEta() - P2.GetEta();
   double DPhi = GetDPhi(P1, P2);

   return sqrt(DPhi * DPhi + DEta * DEta);
}

double GetDPhi(const FourVector P1, const FourVector P2)
{
   double DPhi = P1.GetPhi() - P2.GetPhi();

   if(DPhi > PI)
      DPhi = 2 * PI - DPhi;
   if(DPhi < -PI)
      DPhi = DPhi + 2 * PI;

   return DPhi;
}

double GetMT(const FourVector P1, const FourVector P2)
{
   double PT1 = P1.GetPT();
   double PT2 = P2.GetPT();

   return sqrt(2 * (PT1 * PT2 - P1[1] * P2[1] - P1[2] * P2[2]));
}

double GetMinRadius(const FourVector P1, const FourVector P2, const FourVector P3)   // in eta-phi space
{
   double Eta1 = P1.GetEta();
   double Phi1 = P1.GetPhi();
   double Eta2 = P2.GetEta();
   double Phi2 = P2.GetPhi();
   double Eta3 = P3.GetEta();
   double Phi3 = P3.GetPhi();

   double BestResidual2 = 9999999;
   double Best1 = 0;
   double Best2 = 0;
   double Best3 = 0;
   for(int i1 = 0; i1 <= 1; i1++)
   {
      for(int i2 = 0; i2 <= 1; i2++)
      {
         for(int i3 = 0; i3 <= 1; i3++)
         {
            double AveragePhi = (Phi1 + Phi2 + Phi3 + (i1 + i2 + i3) * 2 * PI) / 3;

            double Residual2 = (AveragePhi - Phi1 - i1 * 2 * PI) * (AveragePhi - Phi1 - i1 * 2 * PI)
               + (AveragePhi - Phi2 - i2 * 2 * PI) * (AveragePhi - Phi2 - i2 * 2 * PI)
               + (AveragePhi - Phi3 - i3 * 2 * PI) * (AveragePhi - Phi3 - i3 * 2 * PI);

            if(Residual2 < BestResidual2)
            {
               Best1 = i1;
               Best2 = i2;
               Best3 = i3;

               BestResidual2 = Residual2;
            }
         }
      }
   }

   return GetMinRadius(Eta1, Phi1 + Best1 * 2 * PI, Eta2, Phi2 + Best2 * 2 * PI, Eta3, Phi3 + Best3 * 2 * PI);
}

double GetMinRadius(const double X1, const double Y1, const double X2, const double Y2,
   const double X3, const double Y3)
{
   // compare two radii:
   //    - common circle radius
   //    - maximum of edge length (divided by two)
   // return the smaller of the two

   // calculate common circle radius
   double C1X0 = 2 * (X1 - X2);
   double C1Y0 = 2 * (Y1 - Y2);
   double C1 = X1 * X1 + Y1 * Y1 - X2 * X2 - Y2 * Y2;
   
   double C2X0 = 2 * (X1 - X3);
   double C2Y0 = 2 * (Y1 - Y3);
   double C2 = X1 * X1 + Y1 * Y1 - X3 * X3 - Y3 * Y3;

   double Distance2 = 99999999;

   if(fabs(C1X0 * C2Y0 - C2X0 * C1Y0) > 1e-8)   // Otherwise three points too close to a straight line
   {
      double M = C1X0 * C2Y0 - C2X0 * C1Y0;
      double MX = C1 * C2Y0 - C2 * C1Y0;
      double MY = C1X0 * C2 - C2X0 * C1;

      double X0 = MX / M;
      double Y0 = MY / M;

      Distance2 = (X1 - X0) * (X1 - X0) + (Y1 - Y0) * (Y1 - Y0);
   }

   // calculate max of edge
   double MaxEdge2 = 0;
   if(MaxEdge2 < (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2))
      MaxEdge2 = (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2);
   if(MaxEdge2 < (X2 - X3) * (X2 - X3) + (Y2 - Y3) * (Y2 - Y3))
      MaxEdge2 = (X2 - X3) * (X2 - X3) + (Y2 - Y3) * (Y2 - Y3);
   if(MaxEdge2 < (X3 - X1) * (X3 - X1) + (Y3 - Y1) * (Y3 - Y1))
      MaxEdge2 = (X3 - X1) * (X3 - X1) + (Y3 - Y1) * (Y3 - Y1);
   MaxEdge2 = MaxEdge2 / 4;

   // minimum of the two
   return sqrt(min(MaxEdge2, Distance2));
}

double GetMR(const FourVector P1, const FourVector P2)
{
   double Temp1 = P1[0] * P2[3] - P1[3] * P2[0];
   double Temp2 = P1[3] - P2[3];
   double Temp3 = P1[0] - P2[0];
   return 2 * sqrt(Temp1 * Temp1 / (Temp2 * Temp2 - Temp3 * Temp3));
}

double GetMRStar(const FourVector P1, const FourVector P2)
{
   double Temp1 = P1[0] + P2[0];
   double Temp2 = P1[3] + P2[3];
   double Temp3 = P1.GetPT() * P1.GetPT() - P2.GetPT() * P2.GetPT();
   double Temp4 = (P1 + P2).GetPT();
   return sqrt((Temp1 * Temp1) - (Temp2 * Temp2) - (Temp3 * Temp3) / (Temp4 * Temp4));
}

double GetMRT(const FourVector P1, const FourVector P2, const FourVector ME)
{
   double Temp1 = ME.GetPT() * (P1.GetPT() + P2.GetPT());
   double Temp2 = ME[1] * (P1[1] + P2[1]) + ME[2] * (P1[2] + P2[2]);
   return sqrt((Temp1 - Temp2) / 2);
}

double GetR(const FourVector P1, const FourVector P2, const FourVector ME)
{
   return GetMRT(P1, P2, ME) / GetMR(P1, P2);
}

double GetRStar(const FourVector P1, const FourVector P2, const FourVector ME)
{
   return GetMRT(P1, P2, ME) / GetMRStar(P1, P2) / GetGammaRStar(P1, P2);
}

double GetGammaRStar(const FourVector P1, const FourVector P2)
{
   double Temp1 = P1[0] + P2[0];
   double Temp2 = P1[3] + P2[3];
   double Temp3 = P1.GetPT() * P1.GetPT() - P2.GetPT() * P2.GetPT();
   double Temp4 = (P1 + P2).GetPT();

   double Upper = Temp1 * Temp1 - Temp2 * Temp2;
   double Lower = Temp1 * Temp1 - Temp2 * Temp2 - Temp3 * Temp3 / Temp4 / Temp4;
   return sqrt(Upper / Lower);
}

double BetaToGamma(double Beta)
{
   return 1 / sqrt(1 - Beta * Beta);
}

double GammaToBeta(double Gamma)
{
   return sqrt(1 - 1 / (Gamma * Gamma));
}

vector<FourVector> SplitIntoGroups(vector<FourVector> &Input, bool ZeroMass)
{
   vector<FourVector> Result;
   
   if(Input.size() == 0)
   {
      Result.push_back(FourVector());
      Result.push_back(FourVector());
      return Result;
   }
   if(Input.size() == 1)
   {
      Result.push_back(Input[0]);
      Result.push_back(FourVector());
      return Result;
   }
   if(Input.size() == 2)
      return Input;

   // let's start with easy (potentially slow) way: try out all possibilities
   // if speed becomes a problem then we'll come up with something else

   int InputSize = Input.size();

   vector<int> Groups(InputSize);
   for(int i = 0; i < InputSize; i++)
      Groups[i] = 0;

   FourVector Group1;
   FourVector Group2;
   double MinMass2 = -1;

   while(Groups[InputSize-1] == 0)   // last one is always in group "0"
   {
      FourVector Vector1Temp;
      FourVector Vector2Temp;

      for(int i = 0; i < InputSize; i++)
      {
         if(Groups[i] == 0)
            Vector1Temp = Vector1Temp + Input[i];
         else
            Vector2Temp = Vector2Temp + Input[i];
      }

      double MinMass2Temp = Vector1Temp.GetMass2() + Vector2Temp.GetMass2();
      if(MinMass2 < 0 || MinMass2Temp < MinMass2)
      {
         MinMass2 = MinMass2Temp;
         Group1 = Vector1Temp;
         Group2 = Vector2Temp;
      }

      Groups[0] = Groups[0] + 1;
      for(int i = 0; i < InputSize - 1; i++)
      {
         while(Groups[i] >= 2)   // just in case something went wrong....an if statement should be enough
         {
            Groups[i] = Groups[i] - 2;
            Groups[i+1] = Groups[i+1] + 1;
         }
      }
   }

   if(ZeroMass == true)
   {
      Group1[0] = Group1.GetP();
      Group2[0] = Group2.GetP();
   }

   Result.push_back(Group1);
   Result.push_back(Group2);

   return Result;
}



