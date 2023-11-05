#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <chrono>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <stdexcept>
#include <mpi.h>
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

void Run_Simulation(std::string filename, int torque_flag, int config_id, int brute_force_flag, int world_rank, int restartflag)
{
	double force = 1.0;

	std::vector<double> gene_lengths, TSSes, promoters;
	std::vector<int> Gene_Directions;
	Read_Setup(filename, gene_lengths, TSSes, Gene_Directions, promoters);

	int clamp0_flag = 1, clamp1_flag = 1;
	double T = 6.0*3600.0;
	int finish_count_limit = -1;
	double topoisomerase = 1.0;
	double rna_degrad = 0.0;
	double barrier_on = 0.0, barrier_off = 0.0; // No barriers needed for this run
	double barrier = 10000.0;
	double nucl_par = 1.0;
	int fileflag = 0;

	int single_pol_II_flag = 0;

	int nucl_file_flag = 1;

	int probe_gene_bodies_flag = 1;

	std::string outputfolder = "outputfiles/RUN_" + std::to_string(config_id) + "_" + std::to_string(brute_force_flag);
	std::string inputfolder = "inputfiles/RUN_" + std::to_string(config_id) + "_" + std::to_string(brute_force_flag);

	Gillespie_Simulation_Multiple_Genes(force, restartflag, brute_force_flag, torque_flag, gene_lengths, TSSes, Gene_Directions, clamp0_flag, clamp1_flag, T, finish_count_limit, promoters, topoisomerase, rna_degrad, barrier_on, barrier_off, barrier, nucl_par, fileflag, outputfolder, inputfolder, world_rank, single_pol_II_flag, nucl_file_flag, probe_gene_bodies_flag);

	return;
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	std::string filename = argv[1];
	int restartflag = std::stoi(argv[2]);
	int config_id = std::stoi(argv[3]);
	int brute_force_flag = std::stoi(argv[4]);
	int torque_flag = 1;

	generator = std::mt19937(std::time(NULL) + config_id + brute_force_flag + world_rank);

	try
	{
		Run_Simulation(filename, torque_flag, config_id, brute_force_flag, world_rank, restartflag);
	}
	catch(const std::invalid_argument &e)
	{
		std::cout << "Exceptional circumstances" << "\n";
		return(0);
	}

	MPI_Finalize();

	return(0);
}
