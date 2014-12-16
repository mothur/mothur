//
//  getmimarkspackagecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/25/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "getmimarkspackagecommand.h"
#include "groupmap.h"


//**********************************************************************************************************************
vector<string> GetMIMarksPackageCommand::setParameters(){
	try {
        //files that have dependancies
        CommandParameter pgroup("group", "InputTypes", "", "", "groupOligos", "none", "none","",false,false); parameters.push_back(pgroup);
        CommandParameter pfile("file", "InputTypes", "", "", "groupOligos", "none", "none","",false,false); parameters.push_back(pfile);
        CommandParameter poligos("oligos", "InputTypes", "", "", "groupOligos", "none", "none","",false,false); parameters.push_back(poligos);
        CommandParameter ppackage("package", "Multiple", "air-host_associated-human_associated-human_gut-human_oral-human_skin-human_vaginal-microbial-miscellaneous-plant_associated-sediment-soil-wastewater-water", "miscellaneous", "", "", "","",false,false,true); parameters.push_back(ppackage);
        CommandParameter prequiredonly("requiredonly", "Boolean", "", "F", "", "", "","",false,false, true); parameters.push_back(prequiredonly);
  		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetMIMarksPackageCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The get.mimarkspackage command creates a mimarks package form with your groups. The required fields are flagged with * characters. \n";
        helpString += "Further documentation on the different packages and required formats can be found here, http://www.mothur.org/wiki/MIMarks_Data_Packages.\n";
		helpString += "The get.mimarkspackage command parameters are: oligos, group, package and requiredonly. oligos or group is required.\n";
		helpString += "The oligos parameter is used to provide your oligos file so mothur can extract your group names.\n";
        helpString += "The group parameter is used to provide your group file so mothur can extract your group names.\n";
        helpString += "The package parameter is used to select the mimarks package you would like to use. The choices are: air, host_associated, human_associated, human_gut, human_oral, human_skin, human_vaginal, microbial, miscellaneous, plant_associated, sediment, soil, wastewater or waterc. Default=miscellaneous.\n";
        helpString += "The requiredonly parameter is used to indicate you only want the required mimarks feilds printed. Default=F.\n";
		helpString += "The get.mimarkspackage command should be in the following format: get.mimarkspackage(oligos=yourOligosFile, package=yourPackage)\n";
		helpString += "get.mimarkspackage(oligos=GQY1XT001.oligos, package=human_gut)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetMIMarksPackageCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "tsv") {  pattern = "[filename],tsv"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetMIMarksPackageCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetMIMarksPackageCommand::GetMIMarksPackageCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["tsv"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "GetMIMarksPackageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetMIMarksPackageCommand::GetMIMarksPackageCommand(string option)  {
	try {
        
		abort = false; calledHelp = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
            vector<string> tempOutNames;
			outputTypes["tsv"] = tempOutNames;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter
			inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
				string path;
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
                }
				
            }
            
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") {  groupfile = "";  abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
            else {  m->setGroupFile(groupfile); inputfile = groupfile; }
            
            file = validParameter.validFile(parameters, "file", true);
			if (file == "not open") {  file = "";  abort = true; }
			else if (file == "not found") { file = ""; }
            else {  inputfile = file; }
            
            oligosfile = validParameter.validFile(parameters, "oligos", true);
			if (oligosfile == "not found")      {	oligosfile = "";	}
			else if(oligosfile == "not open")	{	abort = true;		}
			else {	m->setOligosFile(oligosfile); inputfile = oligosfile; }

            if ((groupfile != "") && (oligosfile != "") && (file != "")) {
                m->mothurOut("[ERROR]: You may not use a group file, file and an oligos file, only one."); m->mothurOutEndLine(); abort = true;
            }

            if ((groupfile == "") && (oligosfile == "") && (file == "")) {
                oligosfile = m->getOligosFile();
                if (oligosfile != "") { inputfile = oligosfile;  m->mothurOut("Using " + oligosfile + " as input file for the oligos parameter."); m->mothurOutEndLine(); }
                else {
                    groupfile = m->getGroupFile();
                    if (groupfile != "") { inputfile = groupfile;  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
                    else {
                        m->mothurOut("[ERROR]: You must provide file, groupfile or oligos file for the get.mimarkspackage command."); m->mothurOutEndLine(); abort = true;
                    }
                }
            }
            
            package = validParameter.validFile(parameters, "package", false);         if (package == "not found") { package = "miscellaneous"; }
            
            if ((package == "air") || (package == "host_associated") || (package == "human_associated") || (package == "human_gut") || (package == "human_oral") || (package == "human_skin") || (package == "human_vaginal") || (package == "microbial") || (package == "miscellaneous") || (package == "plant_associated") || (package == "sediment") || (package == "soil") || (package == "wastewater") || (package == "water")) {}
            else {
                m->mothurOut("[ERROR]: " + package + " is not a valid package selection. Choices are: air, host_associated, human_associated, human_gut, human_oral, human_skin, human_vaginal, microbial, miscellaneous, plant_associated, sediment, soil, wastewater or water. Aborting.\n."); abort = true;
            }
            
            string temp;
			temp = validParameter.validFile(parameters, "requiredonly", false);	if(temp == "not found"){	temp = "F";	}
			requiredonly = m->isTrue(temp);
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "GetMIMarksPackageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetMIMarksPackageCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        if (oligosfile != "") { Oligos oligos(oligosfile); createGroupNames(oligos);  } //createGroupNames fills in group names
        else if (file != "")  { readFile();     }
        else {  GroupMap groupmap(groupfile); groupmap.readMap();
            vector<string> tempGroups = groupmap.getNamesOfGroups();
            for (int i = 0; i < tempGroups.size(); i++) { Groups.insert(tempGroups[i]); }
        }
        
        if (outputDir == "") { outputDir += m->hasPath(inputfile); }
        map<string, string> variables;
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
		string outputFileName = getOutputFileName("tsv", variables);
		
        ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["tsv"].push_back(outputFileName);
        
        out << "#This is a tab-delimited file. Additional Documentation can be found at http://www.mothur.org/wiki/MIMarks_Data_Packages." << endl;
        out << "#Please fill all the required fields indicated with '*'" << endl;
        out << "#Unknown or inapplicable fields can be assigned NA value." << endl;
        out << "#You may add extra custom fields to this template. Make sure all the fields are separated by tabs." << endl;
        out << "#You may remove any fields not required (marked with '*'). Make sure all the fields are separated by tabs." << endl;
        out << "#You can edit this template using Microsoft Excel or any other editor. But while saving the file please make sure to save them as 'TAB-DELIMITED' TEXT FILE." << endl;
        
        if (package == "air") {
            out << "#Environmental:MIMARKS.specimen.air.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*altitude" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	rel_to_oxygen	samp_collect_device	samp_mat_process	*altitude	barometric_press	carb_dioxide	carb_monoxide	chem_administration	elev	humidity	methane	organism_count	oxygen	oxy_stat_samp	perturbation	pollutants	resp_part_matter	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	solar_irradiance	temp	ventilation_rate	ventilation_type	volatile_org_comp	wind_direction	wind_speed" << endl;
            }
        }else if (package == "host_associated") {
            out << "#Environmental:MIMARKS.specimen.host-associated.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host   " << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods    rel_to_oxygen	samp_collect_device	samp_mat_process	*host	age	altitude	blood_press_diast	blood_press_syst	body_habitat	body_product	tissue	chem_administration	depth	diet	disease_stat	dry_mass	elev	family_relationship	genotype	gravidity	height_or_length	host_body_temp	host_color	host_growth_cond	host_shape	host_subject_id	host_taxid	infra_specific_name	infra_specific_rank	last_meal	life_stage	organism_count	oxy_stat_samp	perturbation	phenotype	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	sex	substrate	temp	tot_mass" << endl;
            }
        }else if (package == "human_associated") {
            out << "#Environmental:MIMARKS.specimen.human-associated.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host	rel_to_oxygen	samp_collect_device	samp_mat_process	hiv_stat	ihmc_ethnicity	ihmc_medication_code	age	amniotic_fluid_color	fetal_health_stat	gestation_state	maternal_health_stat	blood_blood_disord	body_product	tissue	body_mass_index	chem_administration	diet	disease_stat	drug_usage	family_relationship	genotype	height	host_body_temp	host_subject_id	last_meal	nose_throat_disord	pulmonary_disord	diet_last_six_month	medic_hist_perform	occupation	organism_count	oxy_stat_samp	perturbation	phenotype	pet_farm_animal	pulse	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	sex	smoker	study_complt_stat	temp	tot_mass	travel_out_six_month	twin_sibling	urine_collect_meth	kidney_disord	urogenit_tract_disor	weight_loss_3_month" << endl;
            }
        }else if (package == "human_gut") {
            out << "#Environmental:MIMARKS.specimen.human-gut.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host	rel_to_oxygen	samp_collect_device	samp_mat_process	ihmc_ethnicity	ihmc_medication_code	age	body_product	tissue	body_mass_index	chem_administration	diet	disease_stat	family_relationship	gastrointest_disord	genotype	height	host_body_temp	host_subject_id	last_meal	liver_disord	medic_hist_perform	occupation	organism_count	oxy_stat_samp	perturbation	phenotype	pulse	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	sex	special_diet	temp	tot_mass" << endl;
            }
        }else if (package == "human_oral") {
            out << "#Environmental:MIMARKS.specimen.human-oral.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host	rel_to_oxygen	samp_collect_device	samp_mat_process	ihmc_ethnicity	ihmc_medication_code	age	body_product	tissue	body_mass_index	chem_administration	diet	disease_stat	family_relationship	genotype	height	host_body_temp	host_subject_id	last_meal	medic_hist_perform	nose_mouth_teeth_throat_disord	occupation	organism_count	oxy_stat_samp	perturbation	phenotype	pulse	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	sex	temp	time_last_toothbrush	tot_mass" << endl;
            }
        }else if (package == "human_skin") {
            out << "#Environmental:MIMARKS.specimen.human-skin.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host	rel_to_oxygen	samp_collect_device	samp_mat_process	ihmc_ethnicity	ihmc_medication_code	age	body_product	tissue	body_mass_index	chem_administration	dermatology_disord	diet	disease_stat	dominant_hand	family_relationship	genotype	height	host_body_temp	host_subject_id	last_meal	medic_hist_perform	occupation	organism_count	oxy_stat_samp	perturbation	phenotype	pulse	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	sex	temp	time_since_last_wash	tot_mass" << endl;
            }
        }else if (package == "human_vaginal") {
            out << "#Environmental:MIMARKS.specimen.human-vaginal.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host	rel_to_oxygen	samp_collect_device	samp_mat_process	hrt	ihmc_ethnicity	ihmc_medication_code	age	birth_control	body_product	tissue	body_mass_index	chem_administration	diet	disease_stat	douche	family_relationship	genotype	gynecologic_disord	height	host_body_temp	host_subject_id	hysterectomy	last_meal	medic_hist_perform	menarche	menopause	occupation	organism_count	oxy_stat_samp	perturbation	phenotype	pregnancy	pulse	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	sex	sexual_act	temp	tot_mass	urogenit_disord" << endl;
            }
        }else if (package == "microbial") {
            out << "#Environmental:MIMARKS.specimen.microbial.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth	*elev" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth	*elev	rel_to_oxygen	samp_collect_device	samp_mat_process	alkalinity	alkyl_diethers	altitude	aminopept_act	ammonium	bacteria_carb_prod	biomass	bishomohopanol	bromide	calcium	carb_nitro_ratio	chem_administration	chloride	chlorophyll	diether_lipids	diss_carb_dioxide	diss_hydrogen	diss_inorg_carb	diss_org_carb	diss_org_nitro	diss_oxygen	glucosidase_act	magnesium	mean_frict_vel	mean_peak_frict_vel	methane	n_alkanes	nitrate	nitrite	nitro	org_carb	org_matter	org_nitro	organism_count	oxy_stat_samp	ph	part_org_carb	perturbation	petroleum_hydrocarb	phaeopigments	phosphate	phosplipid_fatt_acid	potassium	pressure	redox_potential	salinity	samp_size	samp_store_dur	samp_store_loc	samp_store_temp	silicate	sodium	sulfate	sulfide	temp	tot_carb	tot_nitro	tot_org_carb	turbidity	water_content" << endl;
            }
        }else if (package == "miscellaneous") {
            out << "#Environmental:MIMARKS.specimen.miscellaneous.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name   *title  *seq_methods    *lat_lon" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	rel_to_oxygen	samp_collect_device	samp_mat_process	alkalinity	altitude	ammonium	biomass	bromide	calcium	chem_administration	chloride	chlorophyll	current	density	depth	diether_lipids	diss_carb_dioxide	diss_hydrogen	diss_inorg_carb	diss_org_nitro	diss_oxygen	elev	nitrate	nitrite	nitro	org_carb	org_matter	org_nitro	organism_count	oxy_stat_samp	ph	perturbation	phosphate	phosplipid_fatt_acid	potassium	pressure	salinity	samp_size	samp_store_dur	samp_store_loc	samp_store_temp	silicate	sodium	sulfate	sulfide	temp" << endl;
            }
        }else if (package == "plant_associated") {
            out << "#Environmental:MIMARKS.specimen.plant-associated.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*host	rel_to_oxygen	samp_collect_device	samp_mat_process	age	air_temp_regm	altitude	antibiotic_regm	body_product	chem_administration	chem_mutagen	climate_environment	depth	disease_stat	dry_mass	elev	fertilizer_regm	fungicide_regm	gaseous_environment	genotype	gravity	growth_hormone_regm	growth_med	height_or_length	herbicide_regm	host_taxid	humidity_regm	infra_specific_name	infra_specific_rank	life_stage	mechanical_damage	mineral_nutr_regm	non_mineral_nutr_regm	organism_count	oxy_stat_samp	ph_regm	perturbation	pesticide_regm	phenotype	tissue	plant_product	radiation_regm	rainfall_regm	salt_regm	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	season_environment	standing_water_regm	temp	tiss_cult_growth_med	tot_mass	water_temp_regm	watering_regm	wet_mass" << endl;
            }
        }else if (package == "sediment") {
            out << "#Environmental:MIMARKS.specimen.sediment.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth	*elev" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth	*elev	rel_to_oxygen	samp_collect_device	samp_mat_process	alkalinity	alkyl_diethers	aminopept_act	ammonium	bacteria_carb_prod	biomass	bishomohopanol	bromide	calcium	carb_nitro_ratio	chem_administration	chloride	chlorophyll	density	diether_lipids	diss_carb_dioxide	diss_hydrogen	diss_inorg_carb	diss_org_carb	diss_org_nitro	diss_oxygen	glucosidase_act	magnesium	mean_frict_vel	mean_peak_frict_vel	methane	n_alkanes	nitrate	nitrite	nitro	org_carb	org_matter	org_nitro	organism_count	oxy_stat_samp	ph	particle_class	part_org_carb	perturbation	petroleum_hydrocarb	phaeopigments	phosphate	phosplipid_fatt_acid	porosity	potassium	pressure	redox_potential	salinity	samp_size	samp_store_dur	samp_store_loc	samp_store_temp	sediment_type	silicate	sodium	sulfate	sulfide	temp	tidal_stage	tot_carb	tot_nitro	tot_org_carb	turbidity	water_content" << endl;
            }
        }else if (package == "soil") {
            out << "#Environmental:MIMARKS.specimen.soil.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth	*elev" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth	*elev	rel_to_oxygen	samp_collect_device	samp_mat_process	altitude	sieving	cur_land_use	cur_vegetation_meth	cur_vegetation	drainage_class	al_sat	al_sat_meth	heavy_metals_meth	heavy_metals	salinity_meth	extreme_salinity	fao_class	agrochem_addition	crop_rotation	extreme_event	fire	flooding	previous_land_use_meth	previous_land_use	tillage	horizon_meth	horizon	link_class_info	link_climate_info	link_addit_analys	annual_season_precpt	annual_season_temp	microbial_biomass_meth	microbial_biomass	other	ph_meth	ph	pool_dna_extracts	profile_position	samp_size	samp_weight_dna_ext	slope_aspect	slope_gradient	soil_type_meth	soil_type	local_class_meth	local_class	store_cond	texture_meth	texture	tot_n_meth	tot_n	tot_org_c_meth	tot_org_carb	water_content_soil_meth	water_content_soil" << endl;
            }
        }else if (package == "wastewater") {
            out << "#Environmental:MIMARKS.specimen.wastewater.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	rel_to_oxygen	samp_collect_device	samp_mat_process	alkalinity	biochem_oxygen_dem	chem_administration	chem_oxygen_dem	depth	efficiency_percent	emulsions	gaseous_substances	indust_eff_percent	inorg_particles	nitrate	org_particles	organism_count	oxy_stat_samp	ph	perturbation	phosphate	pre_treatment	primary_treatment	reactor_type	samp_size	samp_salinity	samp_store_dur	samp_store_loc	samp_store_temp	secondary_treatment	sewage_type	sludge_retent_time	sodium	soluble_inorg_mat	soluble_org_mat	suspend_solids	temp	tertiary_treatment	tot_nitro	tot_phosphate	wastewater_type" << endl;
            }
        }else if (package == "water") {
            out << "#Environmental:MIMARKS.specimen.water.4.0" << endl;
            if (requiredonly) {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth" << endl;
            }else {
                out << "*sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon    *title  *seq_methods	*depth	rel_to_oxygen	samp_collect_device	samp_mat_process	alkalinity	alkyl_diethers	aminopept_act	ammonium	atmospheric_data	bacteria_carb_prod	biomass	bishomohopanol	bromide	calcium	carb_nitro_ratio	chem_administration	chloride	chlorophyll	current	density	diether_lipids	diss_carb_dioxide	diss_hydrogen	diss_inorg_carb	diss_inorg_nitro	diss_inorg_phosp	diss_org_carb	diss_org_nitro	diss_oxygen	elev	glucosidase_act	light_intensity	magnesium	mean_frict_vel	mean_peak_frict_vel	n_alkanes	nitrate	nitrite	nitro	org_carb	org_matter	org_nitro	organism_count	oxy_stat_samp	ph	part_org_carb	part_org_nitro	perturbation	petroleum_hydrocarb	phaeopigments	phosphate	phosplipid_fatt_acid	photon_flux	potassium	pressure	primary_prod	redox_potential	salinity	samp_size	samp_store_dur	samp_store_loc	samp_store_temp	silicate	sodium	soluble_react_phosp	sulfate	sulfide	suspend_part_matter	temp	tidal_stage	tot_depth_water_col	tot_diss_nitro	tot_inorg_nitro	tot_nitro	tot_part_carb	tot_phosp" << endl;
            }
        }
        
        for (set<string>::iterator it = Groups.begin(); it != Groups.end(); it++) {  out << *it << '\t' << endl; }
        
        out.close();
        
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "execute");
		exit(1);
	}
}
//***************************************************************************************************************

// going to have to rework this to allow for other options --
/*
 file option 1
 
 sfffile1   oligosfile1
 sfffile2   oligosfile2
 ...
 
 file option 2
 
 fastqfile1 oligosfile1
 fastqfile2 oligosfile2
 ...
 
 file option 3
 
 fastqfile  fastqfile   group
 fastqfile  fastqfile   group
 fastqfile  fastqfile   group
 ...
 
 */

int GetMIMarksPackageCommand::readFile(){
	try {
        
        inputfile = file;
        int format = 2;
        
        ifstream in;
        m->openInputFile(file, in);
        
        while(!in.eof()) {
            
            Oligos oligos;
            
            if (m->control_pressed) { return 0; }
            
            string line = m->getline(in);  m->gobble(in);
            vector<string> pieces = m->splitWhiteSpace(line);
            
            string group = "";
            string thisFileName1, thisFileName2; thisFileName1 = ""; thisFileName2 = "";
            if (pieces.size() == 2) {
                thisFileName1 = pieces[0];
                thisFileName2 = pieces[1];
            }else if (pieces.size() == 3) {
                thisFileName1 = pieces[1];
                thisFileName2 = pieces[2];
                string group = pieces[0];
            }else {
                m->mothurOut("[ERROR]: file lines can be 2 or 3 columns. The 2 column files are sff file then oligos or fastqfile then oligos. You may have multiple lines in the file.  The 3 column files are for paired read libraries. The format is groupName, forwardFastqFile reverseFastqFile. \n"); m->control_pressed = true;
            }
            
            if (m->debug) { m->mothurOut("[DEBUG]: group = " + group + ", thisFileName1 = " + thisFileName1 + ", thisFileName2 = " + thisFileName2  + ".\n"); }
            
            if (inputDir != "") {
                string path = m->hasPath(thisFileName2);
                if (path == "") {  thisFileName2 = inputDir + thisFileName2;  }
                
                path = m->hasPath(thisFileName1);
                if (path == "") {  thisFileName1 = inputDir + thisFileName1;  }
            }
            
            //check to make sure both are able to be opened
            ifstream in2;
            int openForward = m->openInputFile(thisFileName1, in2, "noerror");
            
            //if you can't open it, try default location
            if (openForward == 1) {
                if (m->getDefaultPath() != "") { //default path is set
                    string tryPath = m->getDefaultPath() + m->getSimpleName(thisFileName1);
                    m->mothurOut("Unable to open " + thisFileName1 + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openForward = m->openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    thisFileName1 = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (openForward == 1) {
                if (m->getOutputDir() != "") { //default path is set
                    string tryPath = m->getOutputDir() + m->getSimpleName(thisFileName1);
                    m->mothurOut("Unable to open " + thisFileName1 + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openForward = m->openInputFile(tryPath, in4, "noerror");
                    thisFileName1 = tryPath;
                    in4.close();
                }
            }
            
            if (openForward == 1) { //can't find it
                m->mothurOut("[WARNING]: can't find " + thisFileName1 + ", ignoring.\n");
            }else{  in2.close();  }
            
            int openReverse = 1;
            
            ifstream in3;
            openReverse = m->openInputFile(thisFileName2, in3, "noerror");
            
            //if you can't open it, try default location
            if (openReverse == 1) {
                if (m->getDefaultPath() != "") { //default path is set
                    string tryPath = m->getDefaultPath() + m->getSimpleName(thisFileName2);
                    m->mothurOut("Unable to open " + thisFileName2 + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openReverse = m->openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    thisFileName2 = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (openReverse == 1) {
                if (m->getOutputDir() != "") { //default path is set
                    string tryPath = m->getOutputDir() + m->getSimpleName(thisFileName2);
                    m->mothurOut("Unable to open " + thisFileName2 + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openReverse = m->openInputFile(tryPath, in4, "noerror");
                    thisFileName2 = tryPath;
                    in4.close();
                }
            }
            
            if (openReverse == 1) { //can't find it
                m->mothurOut("[WARNING]: can't find " + thisFileName2 + ", ignoring pair.\n");
            }else{  in3.close();  }
            
            
            if ((pieces.size() == 2) && (openForward != 1) && (openReverse != 1)) { //good pair and sff or fastq and oligos
                    oligosfile = thisFileName2;
                    if (m->debug) { m->mothurOut("[DEBUG]: about to read oligos\n"); }
                    oligos.read(oligosfile);
                    createGroupNames(oligos); // adding in groupNames from this file
                format = 2;
            }else if((pieces.size() == 3) && (openForward != 1) && (openReverse != 1)) { //good pair and paired read
                Groups.insert(group);
                format = 3;
            }
        }
        in.close();
        
        inputfile = file;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "readFile");
		exit(1);
	}
}
//**********************************************************************************************************************

set<string> GetMIMarksPackageCommand::createGroupNames(Oligos& oligos) {
    try {
        bool pairedOligos = false;
        
        if (oligos.hasPairedPrimers() || oligos.hasPairedBarcodes()) {      pairedOligos = true;        }
        
        vector<string> groupNames = oligos.getGroupNames();
        if (groupNames.size() == 0) { return Groups;  }
        
        if (pairedOligos) {
            map<int, oligosPair> barcodes = oligos.getPairedBarcodes();
            map<int, oligosPair> primers = oligos.getPairedPrimers();
            for(map<int, oligosPair>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<int, oligosPair>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligos.getPrimerName(itPrimer->first);
                    string barcodeName = oligos.getBarcodeName(itBar->first);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string comboName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeName;
                        }else{
                            if(barcodeName == ""){
                                comboGroupName = primerName;
                            }
                            else{
                                comboGroupName = barcodeName + "." + primerName;
                            }
                        }
                        
                        if(((itPrimer->second).forward+(itPrimer->second).reverse) == ""){
                            comboName = ((itBar->second).forward+"."+(itBar->second).reverse);
                        }else{
                            if(((itBar->second).forward+(itBar->second).reverse) == ""){
                                comboName = (itPrimer->second).forward+"."+(itPrimer->second).reverse;
                            }
                            else{
                                comboName = ((itBar->second).forward+"."+(itBar->second).reverse) + "." + (itPrimer->second).forward+"."+(itPrimer->second).reverse;
                            }
                        }
                        
                        if (comboName != "") {  comboGroupName +=  "_" + comboName;  }
                        Groups.insert(comboGroupName);
                    }
                }
            }
        }else {
            map<string, int> barcodes = oligos.getBarcodes() ;
            map<string, int> primers = oligos.getPrimers();
            for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligos.getPrimerName(itPrimer->second);
                    string barcodeName = oligos.getBarcodeName(itBar->second);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string comboName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeName;
                        }else{
                            if(barcodeName == ""){
                                comboGroupName = primerName;
                            }
                            else{
                                comboGroupName = barcodeName + "." + primerName;
                            }
                        }
                        
                        if(itPrimer->first == ""){
                            comboName = itBar->first;
                        }else{
                            if(itBar->first == ""){
                                comboName = itPrimer->first;
                            }
                            else{
                                comboName = itBar->first + "." + itPrimer->first;
                            }
                        }
                        
                        if (comboName != "") {  comboGroupName +=  "_" + comboName;  }
                        Groups.insert(comboGroupName);
                    }
                }
            }
        }
        
        if (Groups.size() == 0) {
            m->mothurOut("[ERROR]: your oligos file does not contain any group names."); m->mothurOutEndLine(); m->control_pressed = true;
        }
        
        return Groups;
    }
    catch(exception& e) {
        m->errorOut(e, "GetMIMarksPackageCommand", "createGroupNames");
        exit(1);
    }
}
//**********************************************************************************************************************


