
else if (model == "Tr.R.BFS.region.a5.k_froz_omega_froz") { // trained on a rectangular region , x/h=1:16 (streamwise), y/h=0:3 (wall_norm) conditioned on kfroz_min=3.0e-04 with tau = 1./( omega + osmall) 
   nonlinearStress_ = 2 * k_ * (- 1.0*T1*(2.0*I2 - 2.0*I1 + 1.14) - 1.0*T2*(2.0*I2 - 7.86*I1 + 1.355) - 1.0*T4*(7.0*I2 - 3.0) - 1.0*T3*(4.0*I2 + 1.0)*(I1 - 4.0)*(I2 - 1.0*I1 + 2.0) + T1); 
		   //which has a fitness of  (xs=250,ys=150)
        Info<< "\n Tr.R.BFS.region.a5.k_froz_omega_froz no-limiters \n" << endl;
}

else if (model == "Tr.aijx.RijHat.BFS.region.a5.newfrozenAppr") { // trained on a rectangular region , x/h=1:16 (streamwise), y/h=0:3 (wall_norm) conditioned on kfroz_min=1.0e-04 with tau = 1./( omega + osmall) 

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ *( T4*(6.0*I1 - 3.0*I1*I1 + 3.0*I1*I2 - 2.0) - 1.0*T2*(3.0*I2 + 2.0) - 1.0*T1*(I2 - 1.0*I1 + 0.15) + T3*(2.0*I1 - 2.0)*(I1 - 1.0)*(I2*I2 + 3.0*I2 + 3.0)) ; 
   //which has a fitness of 0.00866636869403 (xs=250,ys=150)
}
else {
   nonlinearStress_ =  2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * (T1*(I1 + I2 + (0.15*I2 + 0.15)*(2.0*I2 - 1.0*I1 + 2.302)) - 1.0*T3*(32.0*I2*(4.0*I1 - 8.0*I2 + 4.0) + 4.0) + T2*(2.0*I1 + 2.0*(2.0*I1 + 1.0)*(8.0*I1 - 40.0) + 4.0) + T4*(I1 - 2.0)*(I2 + 2.0)*(I2 - 4.15)) ; 
   //which has a fitness of 0.00392480083941 (xs=250,ys=150)
}
else {

  Rij_ = 2*k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.BFS.region.a5.newfrozenAppr \n" << endl;
}

else if (model == "Tr.aijx.squCyl.region_Ha5_newfrozenaijxRij_hat.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   nonlinearStress_ =  2 * k_ * (T8*(9.0*I2 - 4.0*I1 - 4.0*I4 + 5.0*I5 + 18.82) - 2.239*T4 + T2*(4.0*I1 - 6.0) + T1*(3.0*I1 + I5) + T9*(I5 + 3.0)); //which has a fitness of 0.0120780454326 (xs=160,ys=120,zs=200)

        Info<< "\n Tr.aijx.squCyl.region_Ha5_newfrozenaijxRij_hat.eve \n" << endl;
}

else if (model == "Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   Rij_ = 2 * k_ *(T1*(I2 - 1.0*I1 + I1*(I2 + 3.15)*(I1 - 2.0) + 2.58) + T9*(4.0*I1 + 8.0) + T8*(I2 - 2.43) + T2*(I3 - 0.205*I4 - 1.0*I1*(2.0*I1 - 1.0*I2) + 0.205) + T4*((3.0*I3 - 6.0)*(I2 - 4.0) + 1.0)*(I1 + I3 + 2.0))  ; //which has a fitness of  0.0275656458898 (xs=160,ys=120,zs=200)

        Info<< "\n Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.eve \n" << endl;
}

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr") { // 

dimensionedScalar one
     (
        "one",
        dimensionSet(0, 1, 0, 0, 0, 0 ,0),
        1.
     );
     volScalarField Xswtch0 = pos( this->mesh_.C().component(vector::X) + one*2.0);
     volScalarField Xswtch1 = 1.0 -pos( this->mesh_.C().component(vector::X) - one*12.0);

     volScalarField Yswtch3 = pos( this->mesh_.C().component(vector::Y) + one*6.0);

     volScalarField zswtch0 = pos( this->mesh_.C().component(vector::Z) - one*0.0);
     volScalarField zswtch1 = 1.0 -pos( this->mesh_.C().component(vector::Z) - one*10.0);

if (aijx_term == "2.0")
{
   nonlinearStress_ =  2 * k_ * ( T8*(9.0*I2 - 4.0*I1 - 4.0*I4 + 5.0*I5 + 18.82) - 2.239*T4 + T2*(4.0*I1 - 6.0) + T1*(3.0*I1 + I5) + T9*(I5 + 3.0) ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1 ; //which has a fitness of 0.0120780454326 (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( T1*(I2 - 1.0*I1 + I1*(I2 + 3.15)*(I1 - 2.0) + 2.58) + T9*(4.0*I1 + 8.0) + T8*(I2 - 2.43) + T2*(I3 - 0.205*I4 - 1.0*I1*(2.0*I1 - 1.0*I2) + 0.205) + T4*((3.0*I3 - 6.0)*(I2 - 4.0) + 1.0)*(I1 + I3 + 2.0)  ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1; //which has a fitness of  0.0275656458898 (xs=160,ys=120,zs=200)
}
else {

  Rij_ =  2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr \n" << endl;
}

else if (model == "Tr.aijx.squCyl.region_Ha5_newfrozenaijxRij_hat.2.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   nonlinearStress_ = 2 * k_ * (T1*(8.0*I1 - 0.399) - 1.0*T4*(2.0*I2 - 1.0*I5 + 2.0) + T2*(2.0*I1 - 3.43) - 1.0*T9*(I2 + I3 - 1.0*(I3 - 0.15)*(I5 - 0.15) - 1.178) + T8*(I1 - 1.0*I2 - 1.0*I3 + 6.15)) ; //which has a fitness of 0.00790156728824 (xs=160,ys=120,zs=200)

        Info<< "\n Tr.aijx.squCyl.region_Ha5_newfrozenaijxRij_hat.2.eve \n" << endl;
}

else if (model == "Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.2.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   Rij_ = 2 * k_ * (T1*(4.0*I2 - 2.0*I1 + 2.15) + T2*(I5 - 1.0*I1*(1.911*I3 - 1.0*I1 + I5) + 1.0) + T9*(6.0*I2 + I1*I5 + 18.089) - 1.0*T4*(2.0*I2 - 12.0) + I1*T8)  ; //which has a fitness of  0.0275141837708 (xs=160,ys=120,zs=200)

        Info<< "\n Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.2.eve \n" << endl;
}

else if (model == "Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.3.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   Rij_ = 2 * k_ * (T1*(3.0*I2 - 1.0*I1 + 2.3) + T2*(I2 - 2.0)*(I5 + (I2 - 2.0)*(I4 + 2.0)*(3.0*I1 - 3.0*I2 + I4 + 3.0) + 0.205))  ; //which has a fitness of  0.0276212674604  (xs=160,ys=120,zs=200)

        Info<< "\n Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.3.eve \n" << endl;
}

else if (model == "Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.4.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   Rij_ = 2 * k_ * (T4*(4.0*I1 - 1.0*I2 + 26.0) + T3*(3.0*I2 - 6.0) + I4*T2 + T1*(2.0*I2 - 3.0*I1 + 1.0)*(2.0*I2 - 1.0*I1 + 3.355)) ; //which has a fitness of  0.026800573799  (xs=160,ys=120,zs=200)

        Info<< "\n Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.4.eve \n" << endl;
}

else if (model == "Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.5.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   Rij_ = 2 * k_ * (2.0*T5 + T6*(16.0*I4 + 64.0) + T4*(4.0*I3 - 2.0*I4 + 2.0*I3*I5 + 5.0) + T3*(I5 - 2.0) + T2*(I2 - 2.0*I4*(I2*I5 - 1.0*I4 + 0.43)) - 1.0*T1*(I5 - 2.205)*(I5 - 2.0*I1 + 1.58)) ; //which has a fitness of   0.0260808271908 (xs=160,ys=120,zs=200)

        Info<< "\n Tr.Rij_hat.squCyl.region_Ha5_newfrozenaijxRij_hat.5.eve \n" << endl;
}

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.b") { // 

dimensionedScalar one
     (
        "one",
        dimensionSet(0, 1, 0, 0, 0, 0 ,0),
        1.
     );
     volScalarField Xswtch0 = pos( this->mesh_.C().component(vector::X) + one*2.0);
     volScalarField Xswtch1 = 1.0 -pos( this->mesh_.C().component(vector::X) - one*12.0);

     volScalarField Yswtch3 = pos( this->mesh_.C().component(vector::Y) + one*6.0);

     volScalarField zswtch0 = pos( this->mesh_.C().component(vector::Z) - one*0.0);
     volScalarField zswtch1 = 1.0 -pos( this->mesh_.C().component(vector::Z) - one*10.0);

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * ( T1*(8.0*I1 - 0.399) - 1.0*T4*(2.0*I2 - 1.0*I5 + 2.0) + T2*(2.0*I1 - 3.43) - 1.0*T9*(I2 + I3 - 1.0*(I3 - 0.15)*(I5 - 0.15) - 1.178) + T8*(I1 - 1.0*I2 - 1.0*I3 + 6.15) ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1 ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( T1*(4.0*I2 - 2.0*I1 + 2.15) + T2*(I5 - 1.0*I1*(1.911*I3 - 1.0*I1 + I5) + 1.0) + T9*(6.0*I2 + I1*I5 + 18.089) - 1.0*T4*(2.0*I2 - 12.0) + I1*T8  ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1; //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {

  Rij_ = 2 * k_ *  (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.b \n" << endl;
}

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.c") { // 

dimensionedScalar one
     (
        "one",
        dimensionSet(0, 1, 0, 0, 0, 0 ,0),
        1.
     );
     volScalarField Xswtch0 = pos( this->mesh_.C().component(vector::X) + one*2.0);
     volScalarField Xswtch1 = 1.0 -pos( this->mesh_.C().component(vector::X) - one*12.0);

     volScalarField Yswtch3 = pos( this->mesh_.C().component(vector::Y) + one*6.0);

     volScalarField zswtch0 = pos( this->mesh_.C().component(vector::Z) - one*0.0);
     volScalarField zswtch1 = 1.0 -pos( this->mesh_.C().component(vector::Z) - one*10.0);

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * ( T1*(8.0*I1 - 0.399) - 1.0*T4*(2.0*I2 - 1.0*I5 + 2.0) + T2*(2.0*I1 - 3.43) - 1.0*T9*(I2 + I3 - 1.0*(I3 - 0.15)*(I5 - 0.15) - 1.178) + T8*(I1 - 1.0*I2 - 1.0*I3 + 6.15) ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1 ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( T1*(3.0*I2 - 1.0*I1 + 2.3) + T2*(I2 - 2.0)*(I5 + (I2 - 2.0)*(I4 + 2.0)*(3.0*I1 - 3.0*I2 + I4 + 3.0) + 0.205)  ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1; //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {

  Rij_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.c \n" << endl;
}

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.d") { // 

dimensionedScalar one
     (
        "one",
        dimensionSet(0, 1, 0, 0, 0, 0 ,0),
        1.
     );
     volScalarField Xswtch0 = pos( this->mesh_.C().component(vector::X) + one*2.0);
     volScalarField Xswtch1 = 1.0 -pos( this->mesh_.C().component(vector::X) - one*12.0);
     
     volScalarField Yswtch3 = pos( this->mesh_.C().component(vector::Y) + one*6.0);

     volScalarField zswtch0 = pos( this->mesh_.C().component(vector::Z) - one*0.0);
     volScalarField zswtch1 = 1.0 -pos( this->mesh_.C().component(vector::Z) - one*10.0);

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * ( T1*(8.0*I1 - 0.399) - 1.0*T4*(2.0*I2 - 1.0*I5 + 2.0) + T2*(2.0*I1 - 3.43) - 1.0*T9*(I2 + I3 - 1.0*(I3 - 0.15)*(I5 - 0.15) - 1.178) + T8*(I1 - 1.0*I2 - 1.0*I3 + 6.15) ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1 ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ =  2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( T4*(4.0*I1 - 1.0*I2 + 26.0) + T3*(3.0*I2 - 6.0) + I4*T2 + T1*(2.0*I2 - 3.0*I1 + 1.0)*(2.0*I2 - 1.0*I1 + 3.355)   ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1; //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {

  Rij_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.d \n" << endl;
}

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.e") { // 

dimensionedScalar one
     (
        "one",
        dimensionSet(0, 1, 0, 0, 0, 0 ,0),
        1.
     );
     volScalarField Xswtch0 = pos( this->mesh_.C().component(vector::X) + one*2.0);
     volScalarField Xswtch1 = 1.0 -pos( this->mesh_.C().component(vector::X) - one*12.0);

     volScalarField Yswtch3 = pos( this->mesh_.C().component(vector::Y) + one*6.0);

     volScalarField zswtch0 = pos( this->mesh_.C().component(vector::Z) - one*0.0);
     volScalarField zswtch1 = 1.0 -pos( this->mesh_.C().component(vector::Z) - one*10.0);

if (aijx_term == "2.0")
{
   nonlinearStress_ =  2 * k_ * ( T1*(8.0*I1 - 0.399) - 1.0*T4*(2.0*I2 - 1.0*I5 + 2.0) + T2*(2.0*I1 - 3.43) - 1.0*T9*(I2 + I3 - 1.0*(I3 - 0.15)*(I5 - 0.15) - 1.178) + T8*(I1 - 1.0*I2 - 1.0*I3 + 6.15) ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1 ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( 2.0*T5 + T6*(16.0*I4 + 64.0) + T4*(4.0*I3 - 2.0*I4 + 2.0*I3*I5 + 5.0) + T3*(I5 - 2.0) + T2*(I2 - 2.0*I4*(I2*I5 - 1.0*I4 + 0.43)) - 1.0*T1*(I5 - 2.205)*(I5 - 2.0*I1 + 1.58)  ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1; //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {

  Rij_ =  2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.e \n" << endl;
}

else if (model == "Tr.aijx.squCyl.region_Ha5_newfrozenaijxRij_hat.3.eve") { // trained on a cubic region , x/d=-1:8 (streamwise), y/d=-6:0 (spanwise) and z/d=0:6 (wall_normal) conditioned on k_min=0.001 with tau = 1./( omega + osmall) 

   nonlinearStress_ =  2 * k_ * (T6 - 3.43*T4 - 1.0*T5*(2.0*I3 - 8.0) + T3*(2.0*I1 + I4 - 2.205) + T2*(1.097*I1 - 2.0*I2 + I2*(I5 - 4.0) - 3.247) + 5.464*I1*T1) ; //which has a fitness of 0.00803279443799 (xs=160,ys=120,zs=200)

        Info<< "\n Tr.aijx.squCyl.region_Ha5_newfrozenaijxRij_hat.3.eve \n" << endl;
}

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.f") { // 

dimensionedScalar one
     (
        "one",
        dimensionSet(0, 1, 0, 0, 0, 0 ,0),
        1.
     );
     volScalarField Xswtch0 = pos( this->mesh_.C().component(vector::X) + one*2.0);
     volScalarField Xswtch1 = 1.0 -pos( this->mesh_.C().component(vector::X) - one*12.0);

     volScalarField Yswtch3 = pos( this->mesh_.C().component(vector::Y) + one*6.0);

     volScalarField zswtch0 = pos( this->mesh_.C().component(vector::Z) - one*0.0);
     volScalarField zswtch1 = 1.0 -pos( this->mesh_.C().component(vector::Z) - one*10.0);

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * (  T6 - 3.43*T4 - 1.0*T5*(2.0*I3 - 8.0) + T3*(2.0*I1 + I4 - 2.205) + T2*(1.097*I1 - 2.0*I2 + I2*(I5 - 4.0) - 3.247) + 5.464*I1*T1   ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1 ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( T1*(4.0*I2 - 2.0*I1 + 2.15) + T2*(I5 - 1.0*I1*(1.911*I3 - 1.0*I1 + I5) + 1.0) + T9*(6.0*I2 + I1*I5 + 18.089) - 1.0*T4*(2.0*I2 - 12.0) + I1*T8  ) * Xswtch0* Xswtch1*Yswtch3*zswtch0*zswtch1; //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {

  Rij_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.f \n" << endl;
}
///////////////////////////////////////////////////////////////////////////////
else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.b.cube") { // 
	
if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * ( T1*(8.0*I1 - 0.399) - 1.0*T4*(2.0*I2 - 1.0*I5 + 2.0) + 
		        T2*(2.0*I1 - 3.43) - 1.0*T9*(I2 + I3 - 1.0*(I3 - 0.15)*(I5 - 0.15) - 1.178) + 
			T8*(I1 - 1.0*I2 - 1.0*I3 + 6.15) ) ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
        //Info<< " nonlinearStress_ " << nonlinearStress_ << endl;
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}
if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( T1*(4.0*I2 - 2.0*I1 + 2.15) + 
	    T2*(I5 - 1.0*I1*(1.911*I3 - 1.0*I1 + I5) + 1.0) + 
	    T9*(6.0*I2 + I1*I5 + 18.089) - 1.0*T4*(2.0*I2 - 12.0) + I1*T8 ) ; //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {
        //Info<< " Rij_ " << Rij_ << endl;
  Rij_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}
        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.b.cube \n" << endl;
}
///////////////////////////////////////////////////////////////////////////////

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.f.cube") { // 

if (aijx_term == "2.0")
{
   nonlinearStress_ =  2 * k_ * (  T6 - 3.43*T4 - 1.0*T5*(2.0*I3 - 8.0) + T3*(2.0*I1 + I4 - 2.205) + T2*(1.097*I1 - 2.0*I2 + I2*(I5 - 4.0) - 3.247) + 5.464*I1*T1   ) ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ =  2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ *  ( T1*(4.0*I2 - 2.0*I1 + 2.15) + T2*(I5 - 1.0*I1*(1.911*I3 - 1.0*I1 + I5) + 1.0) + T9*(6.0*I2 + I1*I5 + 18.089) - 1.0*T4*(2.0*I2 - 12.0) + I1*T8  ); //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {

  Rij_ = 2 * k_ *  (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.f.cube \n" << endl;
}

else if (model == "Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.b.No2.SFS2") { // 

dimensionedScalar one
     (
        "one",
        dimensionSet(0, 1, 0, 0, 0, 0 ,0),
        1.
     );

     volScalarField Xswtch2 = pos( this->mesh_.C().component(vector::X) - one*4.222);
     volScalarField Xswtch3 = 1.0 -pos( this->mesh_.C().component(vector::X) - one*12.0);
     
     volScalarField Yswtch3 = pos( this->mesh_.C().component(vector::Y) - one*0.0);
     volScalarField Yswtch4 = 1.0 - pos( this->mesh_.C().component(vector::Y) - one*3.0);

     volScalarField zswtch0 = pos( this->mesh_.C().component(vector::Z) - one*0.0);
     volScalarField zswtch1 = 1.0 -pos( this->mesh_.C().component(vector::Z) - one*3.0);

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * ( T1*(8.0*I1 - 0.399) - 1.0*T4*(2.0*I2 - 1.0*I5 + 2.0) + T2*(2.0*I1 - 3.43) - 1.0*T9*(I2 + I3 - 1.0*(I3 - 0.15)*(I5 - 0.15) - 1.178) + T8*(I1 - 1.0*I2 - 1.0*I3 + 6.15) ) * Xswtch2 *Xswtch3*Yswtch3*Yswtch4*zswtch0*zswtch1 ; //which has a fitness of  (xs=160,ys=120,zs=200)
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ =  2 * k_ * ( T1*(4.0*I2 - 2.0*I1 + 2.15) + T2*(I5 - 1.0*I1*(1.911*I3 - 1.0*I1 + I5) + 1.0) + T9*(6.0*I2 + I1*I5 + 18.089) - 1.0*T4*(2.0*I2 - 12.0) + I1*T8  ) * Xswtch2 *Xswtch3*Yswtch3*Yswtch4*zswtch0*zswtch1 ; //which has a fitness of   (xs=160,ys=120,zs=200)
}
else {

  Rij_ = 2 * k_ *(T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx.RijHat.squCyl.region_Ha5.newfrozenAppr.b.No2.SFS2 \n" << endl;
}

else if (model == "Tr.aijx_RijHat_PH_EVE-MO-52") { // trained on a rectangular region , x/h=1:16 (streamwise), y/h=0:3 (wall_norm) conditioned on kfroz_min=1.0e-04 with tau = 1./( omega + osmall) 

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * (T1*(-0.641+ 4.0*I1+4.0*I2) + T2*(1.166-30.305*I1+8.0*I1*I1)+ T3*(16.0+4.0*I1-2.0*I2)) ;
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ =  2 * k_ * (T1*(-1.049+6.0*I2 )+ T2*(-1.04-2.016*I1+2.0*I1*I1  ) + T3*( 6.653+2.0*I1+2.0*I2-2.0*I2*I2  )) ; //which has a fitness of }
}
else {

  Rij_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx_RijHat_PH_EVE-MO-52 \n" << endl;
}

else if (model == "Tr.aijx_RijHat_PH_EVE-MO-11") { // trained on a rectangular region , x/h=1:16 (streamwise), y/h=0:3 (wall_norm) conditioned on kfroz_min=1.0e-04 with tau = 1./( omega + osmall) 

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * (T1*(-0.641+ 4.0*I1+6.0*I2) + T2*(-26.18*I1-4*I2+8.0*I1*I1)+ T3*(8.0+6.0*I1+2.0*I2)) ;
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ *( T1*(-0.806+6.0*I2 )+ T2*(4.0+1.82*I2-2.0*I1*I1  ) + T3*( 5.064-2.0*I2  )) ; //which has a fitness of }
}
else {

  Rij_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx_RijHat_PH_EVE-MO-11 \n" << endl;
}

else if (model == "Tr.aijx_RijHat_PH_EVE-MO-57") { // trained on a rectangular region , x/h=1:16 (streamwise), y/h=0:3 (wall_norm) conditioned on kfroz_min=1.0e-04 with tau = 1./( omega + osmall) 

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * (T1*(-0.641+ 4.0*I1+4.0*I2) + T2*(4.0-30.0*I1+8.0*I1*I1)+ T3*(12.391+2.0*I1-2.0*I2+2.0*I1*I1)) ;
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * ( T1*(-1.049+6.0*I2 )+ T2*(-1.04*I1-1.04*I2+0.845*I1*I1+0.016*I1*I1*I1-0.016*I1*I2  ) + T3*( 6.591+2.0*I1-2.0*I2*I2  ) ); //which has a fitness of }
}
else {

  Rij_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aijx_RijHat_PH_EVE-MO-57 \n" << endl;
}

else if (model == "Tr.aij.sqcyl.SNH2020-Reg-B") { // trained on a rectangular region , x/h=1:16 (streamwise), y/h=0:3 (wall_norm) conditioned on kfroz_min=1.0e-04 with tau = 1./( omega + osmall) 

if (aijx_term == "2.0")
{
   nonlinearStress_ = 2 * k_ * (520.0*T9 + T4*(16.0*I2 - 46.0) - 1.0*T1*(5.0*I2 + 3.0) - 1.0*T2*(I4 + I5 + 11.3) + T8*(I4 + 3.0)*(1.15*I1 + 0.15*I2 - 4.15) + T1) ;
}
else {
   nonlinearStress_ = 2 * k_ * (T1-T1); 
        Info<< "\n +++++ aijx=0 +++++ \n" << endl;
}

if (Rterm_corr == "2.0")
{
   Rij_ = 2 * k_ * (T1-T1); ; //which has a fitness of }
}
else {

  Rij_ =2 * k_ * (T1-T1); 
        Info<< "\n +++++ Rij_hat_term=0 +++++ \n" << endl;
}

        Info<< "\n Tr.aij.sqcyl.SNH2020-Reg-B \n" << endl;
}
