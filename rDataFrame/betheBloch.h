//sternheimer parametrisation
auto densityEffect = [](long double beta, double gamma){

   double lar_C = 5.215, lar_x0 = 0.201, lar_x1 = 3, lar_a = 0.196, lar_k = 3;
   long double x = log10(beta * gamma);
   
   if( x >= lar_x1 ) return 2*log(10)*x - lar_C;

   else if ( lar_x0 <= x && x < lar_x1) return 2*log(10)*x - lar_C + lar_a * pow(( lar_x1 - x ) , lar_k );

   else return (long double) 0; //if x < lar_x0

};

//mean energy loss independent of pitch
auto betheBloch = [](double energy, double mass_particle){

   //energy needs to be in MeV! is kinetic energy!

   double K,rho,Z,A, charge, me, I, gamma,  momentum ,wmax, pitch;
   long double beta;
    K = 0.307;
    rho = 1.4;
    charge = 1;
    Z = 18;
    A = 39.948;
    I = pow(10,-6)*10.5*18; //MeV
    me = 0.51; //MeV me*c^2
    pitch = 1;
    
    //momentum = sqrt( pow(energy,2) - pow(mass_particle,2));
    //beta = momentum/sqrt(pow(mass_particle,2) + pow(momentum,2));
    //gamma =  1/sqrt(1 - pow(beta,2));
    
    gamma = (energy + mass_particle) / mass_particle;
    beta = sqrt( 1 - 1/pow(gamma,2));

    wmax = 2*me*pow(beta,2)*pow(gamma,2)/(1+2*gamma*me/mass_particle + pow(me,2)/pow(mass_particle,2));
    
    
    double dEdX;
    //multiply by rho to have dEdX MeV/cm in LAr

    dEdX = pitch*(rho*K*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );

   return dEdX;
};

auto betheBloch_mpv = [](double energy, double mass_particle){

   //energy needs to be in MeV!

   double K,rho,Z,A, charge, me, I, gamma,  momentum , j, pitch, xi;
  long double beta;
    K = 0.307;
    rho = 1.4;
    charge = 1;
    Z = 18;
    A = 39.948;
    I = pow(10,-6)*10.5*18; //MeV
    j = 0.2;
    pitch = 9*rho; //0.5*rho; // g*cm^-2
    me = 0.51; //MeV me*c^2
    
    //momentum = sqrt( pow(energy,2) - pow(mass_particle,2));
    //beta = momentum/sqrt(pow(mass_particle,2) + pow(momentum,2));
    //gamma =  1/sqrt(1 - pow(beta,2));
    
    gamma = (energy + mass_particle) / mass_particle;
    beta = sqrt( 1 - 1/pow(gamma,2));

    xi = ( K/2 )*( Z/A )* ( pitch / pow(beta,2));

    double eloss_mpv; //multiply by rho bc in PDG MPV is also given as MeV / g*cm^2

    eloss_mpv = xi*(log( 2*me*pow(gamma,2)*pow(beta,2) / I ) + log( xi / I ) + j - pow(beta,2) - densityEffect( beta, gamma ) );



   return eloss_mpv;
};

auto betaGamma = [](double energy){
   long double gamma, beta;
   double mass_particle = 139;
   gamma = (energy + mass_particle) / mass_particle;
   beta = sqrt( 1 - 1/pow(gamma,2));
   
   return beta*gamma;

};


