#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <chrono>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <stdexcept>
//#include <mpi.h>
#include "../code/linterp.h"
#include "../code/Constants.hpp"
#include "../code/Model.hpp"
#include "../code/Simulation_Setups.hpp"

void Read_Setup(std::string filename, std::vector<double> &gene_lengths, std::vector<double> &TSSes, std::vector<int> &Gene_Directions, std::vector<double> &promoters)
{
	std::ifstream f;
	f.open(filename);
	double gene_length, TSS, promoter;
	int direction;
	std::string gene_name;

	while(f >> gene_name >> TSS >> gene_length >> direction >> promoter)
	{
		gene_lengths.push_back(gene_length);
		TSSes.push_back(TSS);
		Gene_Directions.push_back(direction);
		promoters.push_back(promoter);
	}

	f.close();

	return;
}

void Run_Simulation(int torque_flag, double segment, int segment_id)
{
	double force = 1.0;
	if(torque_flag == 0)
	{
		double sigma = -0.1, torque;
		std::ofstream sigma_file, torque_file;
		sigma_file.open("prokaryote/sigma_" + std::to_string(segment_id) + ".log");
		torque_file.open("prokaryote/torque_" + std::to_string(segment_id) + ".log");

		while(sigma < 0.1)
		{
			torque = Get_Prokaryotic_Torque(sigma, force, segment);
			sigma_file << sigma << "\n";
			torque_file << torque << "\n";

			sigma += 0.001;
		}

		sigma_file.close();
		torque_file.close();

		return;
	}

	InterpMultilinear<2, double> interp_ML = Setup_Interp();
	InterpMultilinear<1, double> interp_ML_Cutoff_0 = Setup_Interp_Cutoffs("../torque_interp/sigma_s.log");
	InterpMultilinear<1, double> interp_ML_Cutoff_1 = Setup_Interp_Cutoffs("../torque_interp/sigma_p.log");

	std::ofstream psi_file, sigma_file, torque_file, sigma_s_file, sigma_p_file;
	psi_file.open("eukaryote/psi_" + std::to_string(segment_id) + ".log");
	sigma_file.open("eukaryote/sigma_" + std::to_string(segment_id) + ".log");
	torque_file.open("eukaryote/torque_" + std::to_string(segment_id) + ".log");
	sigma_s_file.open("eukaryote/sigma_s_" + std::to_string(segment_id) + ".log");
	sigma_p_file.open("eukaryote/sigma_p_" + std::to_string(segment_id) + ".log");

	array<double, 1> args;
	array<double, 2> args_2;

	double psi = 0.0, sigma = -0.1;
	int flag = 0;
	while(psi < 0.995)
	{
		psi_file << psi << "\n";

		args = {psi};
		sigma_s_file << interp_ML_Cutoff_0.interp(args.begin()) << "\n";
		sigma_p_file << interp_ML_Cutoff_1.interp(args.begin()) << "\n";

		sigma = -0.1;
		while(sigma < 0.1)
		{
			if(flag == 0)
			{
				sigma_file << sigma << "\n";
			}
			args_2 = {psi, sigma};
			torque_file << Get_Eukaryotic_Torque(interp_ML, interp_ML_Cutoff_0, interp_ML_Cutoff_1, psi, sigma, force, segment) << "\t";
			sigma += 0.001;
		}
		torque_file << "\n";
		psi += 0.01;

		flag = 1;
	}

	psi_file.close();
	sigma_file.close();
	torque_file.close();
	sigma_s_file.close();
	sigma_p_file.close();

	return;
}

int main(int argc, char *argv[])
{
//	MPI_Init(NULL, NULL);
	int world_size;
//	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank = 0;
//	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int torque_flag = std::stoi(argv[1]);
	double segment = std::stod(argv[2]);
	segment = segment*0.34;
	int segment_id = std::stoi(argv[3]);

	try
	{
		Run_Simulation(torque_flag, segment, segment_id);
	}
	catch(const std::invalid_argument &e)
	{
		std::cout << "Exceptional circumstances" << "\n";
		return(0);
	}

//	MPI_Finalize();

	return(0);
}
