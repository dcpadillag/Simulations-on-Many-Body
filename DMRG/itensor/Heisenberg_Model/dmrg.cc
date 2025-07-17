#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto S = input.getReal("S");

    auto nsweeps = input.getInt("nsweeps");
    auto J = input.getReal("J",1);

    auto quiet = input.getYesNo("quiet",false);
    auto writedim = input.getInt("writedim",100);
    
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = CustomSpin(N,{"S=",S,"ConserveQNs=",false});

    //
    // Create the Hamiltonian using AutoMPO
    //
    auto ampo = AutoMPO(sites);

    for(int j = 1; j < N; ++j)
    	{
    	ampo += 0.5*J,"S+",j,"S-",j+1;
    	ampo += 0.5*J,"S-",j,"S+",j+1;
    	ampo +=     J,"Sz",j,"Sz",j+1;
    	}

    ampo += 0.5*J,"S+",N,"S-",1;
    ampo += 0.5*J,"S-",N,"S+",1;
    ampo +=     J,"Sz",N,"Sz",1;
    auto H = toMPO(ampo);

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
//    auto state = InitState(sites);
//    for(int i = 1; i <= N; ++i) 
//        {
//        if(i%2 == 1) state.set(i,"1");
//        else         state.set(i,"-1");
//        }
    
//    auto psi0 = MPS(state);
    auto psi0 = randomMPS(sites);
    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet=",quiet,"WriteDim=",writedim});
    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy per Site per Spin= %.10f", (2.*energy)/(J*N*S*S));
//    printfln("\nGround State Energy per Site = %.10f", energy/N);

    return 0;
    }
