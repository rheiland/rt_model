/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./radiot.h"

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom( parameters.ints("random_seed") ); 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; // true; 
	
	// set to no motility for cancer cells 
	cell_defaults.phenotype.motility.is_motile = false; 
	
	// use default proliferation and death 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[0] = 10; 
	cell_defaults.phenotype.secretion.saturation_densities[0] = 38; 

	// set the default cell type to o2-based proliferation with the effect of the 
	// on oncoprotein, and secretion of the immunostimulatory factor 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_therapy; 
	
	// add the extra bit of "attachment" mechanics 
	// cell_defaults.functions.custom_cell_rule = extra_elastic_attachment_mechanics; 
	
	// change the max cell-cell adhesion distance 
	// cell_defaults.phenotype.mechanics.set_relative_maximum_adhesion_distance(parameters.doubles("max_relative_cell_adhesion_distance") );
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	cell_defaults.custom_data.add_variable( "exposure" , "dimensionless", 0.0 ); 
	cell_defaults.custom_data[ "exposure" ] = 0.0;   // AUC  (TODO)

	// add custom data 
		
	Parameter<double> paramD; 
	
	// for therapy 
	
	// shall we introduce some param(s) for radiotherapy?

	// paramD = parameters.doubles["damage_rate"]; 
	// cell_defaults.custom_data.add_variable( "damage rate" , paramD.units, paramD.value ); 
	// paramD = parameters.doubles["repair_rate"]; 
	// cell_defaults.custom_data.add_variable( "repair rate" , paramD.units, paramD.value ); 
	// paramD = parameters.doubles["drug_death_rate"]; 	
	// cell_defaults.custom_data.add_variable( "drug death rate" , paramD.units, paramD.value ); 
	// cell_defaults.custom_data.add_variable( "damage" , "dimensionless", 0.0 ); 
	
	return; 
}

void setup_microenvironment( void )
{
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "WARNING: overriding from 3-D to 2-D" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// default_microenvironment_options.calculate_gradients = true; 
	default_microenvironment_options.calculate_gradients = false; 
	
	microenvironment.add_density( "therapeutic", "dimensionless" ); 
	microenvironment.diffusion_coefficients[1] = 1e3; 
	microenvironment.decay_rates[1] = 0.15625; 	
	
	// let BioFVM use oxygen as the default 
	
	default_microenvironment_options.use_oxygen_as_first_field = true; 

	// set Dirichlet conditions 
	
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	default_microenvironment_options.Dirichlet_condition_vector[0] = 38; // physioxic conditions 
	// default_microenvironment_options.Dirichlet_condition_vector[1] = 0; 
	// default_microenvironment_options.Dirichlet_condition_vector[2] = 0; 
	
	// default_microenvironment_options.Dirichlet_activation_vector[1] = false;  // no Dirichlet for the therapeutic  
	
	// set initial conditions 
	// default_microenvironment_options.initial_condition_vector = { 38.0 , 0.0, 0.0 }; 
	default_microenvironment_options.initial_condition_vector = { 38.0 , 0.0}; 
			
	initialize_microenvironment(); 	

	return; 
}	

void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles("tumor_radius"); // 200.0; 
	
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 
	
	int n = 0; 
	while( y < tumor_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell(); // tumor cell 
			pCell->assign_position( x , y , 0.0 );
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( x , -y , 0.0 );				
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( -x , y , 0.0 );
				
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(); // tumor cell 
					pCell->assign_position( -x , -y , 0.0 );
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	
	return; 
}

// done 
std::vector<std::string> radiot_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4, "black" ); 
	
	static int damage_i = pCell->custom_data.find_variable_index( "damage" ); 
	// static double max_damage = 1.0 * cell_defaults.custom_data["damage rate"] / (1e-16 + cell_defaults.custom_data[ "repair rate" ] );
	
	// cargo cell 
	// if( pCell->type == 1 )
	// {
	// 	output[0] = "blue";
	// 	output[1] = "blue";
	// 	output[2] = "blue"; 
	// 	output[3] = "none"; // no nuclear outline color 
	// 	return output;
	// }
	
	// // worker cell 
	// if( pCell->type == 2 )
	// {
	// 	output[0] = "red";
	// 	output[1] = "red";
	// 	output[2] = "red"; 
	// 	output[3] = "none"; // no nuclear outline color 
	// 	return output;
	// }
	
	// apoptotic tumor - cyan 
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - cyan
	{
		output[0] = "cyan";
		output[2] = "darkcyan"; 
		return output; 
	}	
	
	// Necrotic tumor - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
		return output; 
	}		
	
	// live tumor -- shade by level of damage 
	
	
	// if live: color by damage 
	if( pCell->phenotype.death.dead == false )
	{
		// TODO - fix for radiotherapy
		// int damage = (int) round( pCell->custom_data[damage_i] * 255.0 / max_damage ); 
		int damage = 128; 
		
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)" , damage , 255-damage , damage );
		output[0].assign( szTempString );
		output[1].assign( szTempString );
		sprintf( szTempString , "rgb(%u,%u,%u)" , damage/4 , (255-damage)/4 , damage/4 );
		output[2].assign( szTempString );
	}
	return output; 
}

/* TUMOR CELL RULES */ 
void tumor_cell_phenotype_with_therapy( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 

	// TODO - fix for radiotherapy
	// static int damage_i = pCell->custom_data.find_variable_index( "damage" ); 
	static int exposure_i = pCell->custom_data.find_variable_index( "exposure" ); 

	// static int damage_rate_i = pCell->custom_data.find_variable_index( "damage rate" ); 
	// static int repair_rate_i = pCell->custom_data.find_variable_index( "repair rate" ); 
	// static int death_rate_i = pCell->custom_data.find_variable_index( "drug death rate" ); 
	
	static int chemo_i = microenvironment.find_density_index( "therapeutic" ); 
	
	static int apoptosis_model_index = phenotype.death.find_death_model_index( "apoptosis" );	
	
	// static double max_damage = 1.0 * cell_defaults.custom_data["damage rate"] / (1e-16 + cell_defaults.custom_data[ "repair rate" ] );
	
	// if I'm dead, don't bother. disable my phenotype rule
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
		return; 
	}
	
	// first, vary the cell birth and death rates with oxygenation
	
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// then update the cell damage 
	
	// dD/dt = alpha*c - beta-D by implicit scheme 
	
	double temp = pCell->nearest_density_vector()[chemo_i];
	
	// reuse temp as much as possible to reduce memory allocations etc. 
	// temp *= dt; 
	// temp *= pCell->custom_data[damage_rate_i]; 
	
	// pCell->custom_data[damage_i] += temp; // d_prev + dt*chemo*damage_rate 
	
	// temp = pCell->custom_data[repair_rate_i];
	// temp *= dt; 
	// temp += 1.0; 
	// pCell->custom_data[damage_i] /= temp;  // (d_prev + dt*chemo*damage_rate)/(1 + dt*repair_rate)
	
	// then, see if the cell undergoes death from the therapy 
	
	temp = dt; 
	// temp *= pCell->custom_data[damage_i]; 
	// temp *= pCell->custom_data[death_rate_i]; 
	// temp /= max_damage; // dt*(damage/max_damage)*death_rate 


	if( UniformRandom() <= temp )
	{
		pCell->start_death( apoptosis_model_index );
		pCell->functions.update_phenotype = NULL; 		
		pCell->functions.custom_cell_rule = NULL; 
	}

	return; 
}

