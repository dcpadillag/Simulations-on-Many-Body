#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto Nup = input.getInt("Nup",N/2); //number of particles with spin up, default is N/2 (half filling)
    auto Ndown = input.getInt("Ndown",N/2); //number of particles with spin down, default is N/2 (half filling)

    auto nsweeps = input.getInt("nsweeps");
    auto t1 = input.getReal("t1",1);
    auto U = input.getReal("U",0);
    auto Delta = input.getReal("D",0);

    //
    // Print the final energy reported by DMRG
    //

    println("Load sites from disk");
    Electron sites;
    readFromFile(format("sites_file_%d%e",Nup,Ndown),sites);
    println("Load psi from disk");
    MPS psi(sites);
    readFromFile(format("psi_file_%d%e",Nup,Ndown),psi);
    println("Load complete");

///*    //
    // Measure local operator: spin densities
    //
    Vector upd(N),dnd(N);
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Nup",j)*psi(j));
        dnd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Ndn",j)*psi(j));
        }

    println("Up Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,upd(j));
    println();

    println("Dn Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,dnd(j));
    println();

    println("Total Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,(upd(j)+dnd(j)));
    println();
//*/

    //
    // Measure two point correlation: bond order parameter
    //
    
    float BO = 0.0;

    for(int i=1; i<N; ++i)
	{
	int j = i+1;

	auto F_i = op(sites,"F",i);
	auto F_j = op(sites,"F",j);
	auto Adag_i_up = op(sites,"Adagup",i);
	auto Adag_i_dn = op(sites,"Adagdn",i);
	auto A_j_up = op(sites,"Aup",j);
	auto A_j_dn = op(sites,"Adn",j);


	psi.position(i);
	auto psidag = dag(psi);
	psidag.prime();

	auto li_1 = leftLinkIndex(psi,i);

	auto F_up = prime(psi(i),li_1)*F_i;
	F_up.noPrime("Site");
	auto AdagFi_Aj_up = F_up*Adag_i_up*psidag(i);
	auto Adagi_FAj_dn = prime(psi(i),li_1)*Adag_i_dn*psidag(i);


	auto lj = rightLinkIndex(psi,j);
	AdagFi_Aj_up *= prime(psi(j),lj);
	AdagFi_Aj_up *= A_j_up;
	AdagFi_Aj_up *= psidag(j);

	auto F_dn= F_j*psidag(j);
	F_dn.prime("Site");
	Adagi_FAj_dn *= prime(psi(j),lj);
	Adagi_FAj_dn *= A_j_dn;
	Adagi_FAj_dn *= F_dn;

	if(i%2==1)
		{
		BO += 2*elt(AdagFi_Aj_up)+2*elt(Adagi_FAj_dn);
		}
	else
		{
		BO += -2*elt(AdagFi_Aj_up)-2*elt(Adagi_FAj_dn);
		}
//	println("Cdagupi_Cupj=",elt(AdagFi_Aj_up));
//	println("Cdagdni_Cdnj=",elt(Adagi_FAj_dn));
//	println("Cdagupj_Cupi=",-elt(AFi_Adagj_up));
//	println("Cdagdnj_Cdni=",-elt(Ai_FAdagj_dn));
	}

    println("BO = ", abs(BO)/(N-1));
    return 0;
    }
