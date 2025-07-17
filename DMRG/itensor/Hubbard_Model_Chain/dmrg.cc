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
    int Npart = Nup + Ndown; 

    auto nsweeps = input.getInt("nsweeps");
    auto t1 = input.getReal("t1",1);
    auto U = input.getReal("U",0);
    auto Delta = input.getReal("D",0);
    auto quiet = input.getYesNo("quiet",false);
    auto writedim = input.getInt("writedim",100);
    
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = Electron(N);

    //
    // Create the Hamiltonian using AutoMPO
    //
    auto ampo = AutoMPO(sites);

    for(int i = 1; i <= N; ++i) 
        {
        ampo += U,"Nupdn",i;
	if(i%2==1)
		{
		ampo += Delta,"Ntot",i;
		}
	else
		{
		ampo += -Delta,"Ntot",i;
		}
	}

    for(int b = 1; b < N; ++b)//Chain
	{
	ampo += -t1,"Cdagup",b,"Cup",b+1;
	ampo += -t1,"Cdagup",b+1,"Cup",b;
	ampo += -t1,"Cdagdn",b,"Cdn",b+1;
	ampo += -t1,"Cdagdn",b+1,"Cdn",b;
	}

    auto H = toMPO(ampo);

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    int p_up = Nup;
    int p_down = Ndown;
    int p = Npart;

    for(int i = N; i >= 1; --i)
       {
        if(p > i)
            {
//            println("Doubly occupying site ",i);
            state.set(i,"UpDn");
            p_up -= 1;
            p_down -= 1;
	    p -= 2;
            }
        else if(p > 0)
            {
            if (i%2 == 1 && p_up > 0)
               {
//               println("Singly occupying site Up ",i);            
               state.set(i,"Up");
               p_up -= 1;
	       p -= 1;
               }
            else if (i%2 == 1 && p_up==0)
               {
//               println("Singly occupying site Dn ",i);            
               state.set(i,"Dn");
               p_down -= 1;
               p -= 1;
               }
            if (i%2 == 0 && p_down > 0)
               {
//               println("Singly occupying site Dn ",i);            
               state.set(i,"Dn");
               p_down -= 1;
               p -= 1;
               }
            else if (i%2 == 0 && p_down==0)
               {
//               println("Singly occupying site Up ",i);            
               state.set(i,"Up");
               p_up -= 1;
               p -= 1;
               }
            }
         else
            {
//            println("Empty site Up ",i);            
            state.set(i,"Emp");
            }
       }

//    auto psi0 = MPS(state);
    auto psi0 = randomMPS(state);

    Print(totalQN(psi0));

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet=",quiet,"WriteDim=",writedim});

    //
    //Save sites and psi MPS to disk (correlations calculations)
    //
    if(Nup==Ndown)
	{
	println("Saving sites to disk");
//	writeToFile(format("Measurements/sites_file_%d%e",Nup,Ndown),sites);
	writeToFile(format("sites_file_%d%e",Nup,Ndown),sites);
	println("Saving psi to disk");
//	writeToFile(format("Measurements/psi_file_%d%e",Nup,Ndown),psi);
	writeToFile(format("psi_file_%d%e",Nup,Ndown),psi);
	println("Save complete");
	}

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);

    return 0;
    }
