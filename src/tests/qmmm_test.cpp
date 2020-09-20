#include <cmath>
#include <cassert>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include <numeric>

#include "io_tools/read_nrg.h"
#include "io_tools/write_nrg.h"

#include "bblock/system.h"

//#define NUMGRADS
//#define PRINT_GRADS
//#define PRINT_VIRIAL
namespace {

static std::vector<bblock::System> systems;
static std::vector<bblock::System> MM_system;
static std::vector<bblock::System> QM_system;
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.nrg> [mbx.json]" << std::endl;
        return 0;
    }

    try {
        std::ifstream ifs(argv[1]);

        if (!ifs) {
            throw std::runtime_error("could not open the NRG file");
        }

        tools::ReadNrg(argv[1], systems);
        ifs.close();

    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    std::vector<double> box;

    if (argc > 2) {
        systems[0].SetUpFromJson(argv[2]);
    } else {
        systems[0].SetUpFromJson();
    }


   // double qmmm_en = systems[0].QMMM_setup(true); 

    std::vector<int> qm_ind = systems[0].get_qm_indeces();

    size_t tnm=systems[0].GetNumMon();

    std::vector<int> mm_ind;

    std::vector<int> total_ind(tnm);

    int tot_num_mon = 0;

    for (size_t pp=0; pp<tnm; pp++) {
         tot_num_mon +=1;
    }

    std::iota(std::begin(total_ind), std::end(total_ind), 1);

    for (int pp=1; pp<=tot_num_mon; pp++) {
       if (std::find(qm_ind.begin(), qm_ind.end(),pp) != qm_ind.end()) {
          //std::cout<< "QM " << pp<< std::endl;
          ;
       } else {
           //std::cout<< "MM " << pp<< std::endl;
           mm_ind.push_back(pp);
       }

    }
    std::cout<< "mm system size " << mm_ind.size() <<" "<< mm_ind[0] << std::endl;

    tools::ReadNrg(argv[1], MM_system,&qm_ind); // MM system
    tools::ReadNrg(argv[1], QM_system,&mm_ind); // QM system 

    //std::cout << "potato" << std::endl;
    std::vector<double> MM_xyz = MM_system[0].GetRealXyz();
    std::vector<double> QM_xyz = QM_system[0].GetRealXyz();
    std::vector<double> MM_virt_xyz = MM_system[0].GetXyz();

    std::vector<double> MM_charges = MM_system[0].GetRealCharges();
    std::vector<double> QM_charges = QM_system[0].GetRealCharges();
    
    std::vector<double> MM_virt_charges = MM_system[0].GetCharges();

    std::vector<std::string> MM_atomnames = MM_system[0].GetRealAtomNames();
    std::vector<std::string> QM_atomnames = QM_system[0].GetRealAtomNames();

    std::vector<std::string> MM_virt_atomnames = MM_system[0].GetAtomNames();


    double QM_total_charge = std::accumulate(QM_charges.begin(), QM_charges.end(),
                               decltype(QM_charges)::value_type(0));

    //std::cout << "turnip" << std::endl;


    //std::cout << "MM SYSTEM " << std::endl;
    //for (size_t pp=0; pp < MM_xyz.size(); pp++) {
    //    std::cout << MM_xyz[pp] << std::endl;
    //}
    //std::cout << "QM SYSTEM " << std::endl;
    //for (size_t pp=0; pp < QM_xyz.size(); pp++) {
    //    std::cout << QM_xyz[pp] << std::endl;
    //}
    //std::cout << "TURNIP " << std::endl;
    
    std::vector<std::string> qm_auxparams = systems[0].get_qm_auxparams();
    if (systems[0].get_qm_code() == "qchem") {

            
        std::ofstream qin("qmmm_input.qin");

        qin << "$molecule\n";
        qin << systems[0].get_qm_charge() << " " << systems[0].get_qm_spin() << "\n";

        for (size_t pp=0; pp<QM_xyz.size()/3; pp++) {
             qin << QM_atomnames[pp] << " " << QM_xyz[3*pp] << " " << QM_xyz[3*pp+1] << " " << QM_xyz[3*pp+2] << "\n";
	}
        qin << "$end\n";


        qin << "\n";


        qin << "$rem\n";
        qin << "jobtype sp\n";
        qin << "method " << systems[0].get_qm_theory() << "\n";
        qin << "basis " << systems[0].get_qm_basis() << "\n";
        qin << "MBX true\n";
        qin << "MBX_molecule " << mm_ind.size() << "\n";
       
   

        for (size_t pp=0; pp<qm_auxparams.size(); pp++) {
             qin << qm_auxparams[pp] << "\n";
        }
        
        qin << "$end\n";
 
        qin << "\n"; 
      
        qin << "$MBX\n";
        for (size_t pp=0; pp<MM_virt_xyz.size()/3; pp++) {
             qin << MM_virt_atomnames[pp] << " " << MM_virt_xyz[3*pp] << " " << MM_virt_xyz[3*pp+1] << " " << MM_virt_xyz[3*pp+2] << " " << MM_virt_charges[pp] <<"\n";
        }
        qin << "$end\n";

        qin.close();

    } else if (systems[0].get_qm_code() == "psi4") {

        std::ofstream psi4in("qmmm_job.py");
        

        psi4in << "import psi4\n";

        psi4in << "psi4.set_memory(\'500 MB\')\n";

        psi4in << " molsys = psi4.geometry    

         std::cout << systems[0].get_qm_code() << " is not implemented yet. exiting. " << std::endl;
         std::exit(0);
         std::cout << " POTATO YOU SHOULD NOT BE SEEING THIS MESSAGE. SOMETHING IS WRONG. " << std::endl;

    } else if (systems[0].get_qm_code() == "chronus") {

         std::cout << systems[0].get_qm_code() << " is not implemented yet. exiting. " << std::endl;
         std::exit(0);
         std::cout << " POTATO YOU SHOULD NOT BE SEEING THIS MESSAGE. SOMETHING IS WRONG. " << std::endl;
    } else if (systems[0].get_qm_code() == "gaussian") {

         std::cout << systems[0].get_qm_code() << " is not implemented yet. exiting. " << std::endl;
         std::exit(0);
         std::cout << " POTATO YOU SHOULD NOT BE SEEING THIS MESSAGE. SOMETHING IS WRONG. " << std::endl;

    } else {

         std::cout << systems[0].get_qm_code() << " is not implemented yet. exiting. " << std::endl;
         std::exit(0);
         std::cout << " POTATO YOU SHOULD NOT BE SEEING THIS MESSAGE. SOMETHING IS WRONG. " << std::endl;
    }


} 
