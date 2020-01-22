//
// Created by Bernardo on 11/23/2017.
//

#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>
#include <fstream>
#include <omp.h>
#ifdef  PARALLEL
#include <mpi.h>
#include "../../Borg/borgms.h"
#endif
#include "Triangle.h"
#include "../Controls/SeasonalMinEnvFlowControl.h"
#include "../Controls/Custom/FallsLakeMinEnvFlowControl.h"
#include "../Controls/StorageMinEnvFlowControl.h"
#include "../Controls/InflowMinEnvFlowControl.h"
#include "../Controls/FixedMinEnvFlowControl.h"
#include "../Controls/Custom/JordanLakeMinEnvFlowControl.h"
#include "../SystemComponents/WaterSources/AllocatedReservoir.h"
#include "../SystemComponents/WaterSources/Quarry.h"
#include "../SystemComponents/WaterSources/Relocation.h"
#include "../SystemComponents/WaterSources/ReservoirExpansion.h"
#include "../SystemComponents/WaterSources/WaterReuse.h"
//#include "../SystemComponents/WaterSources/SequentialJointTreatmentExpansion.h"
#include "../DroughtMitigationInstruments/Transfers.h"
#include "../DroughtMitigationInstruments/InsuranceStorageToROF.h"
#include "../Simulation/Simulation.h"
#include "../SystemComponents/Bonds/LevelDebtServiceBond.h"
#include "../SystemComponents/Bonds/BalloonPaymentBond.h"

#ifdef PARALLEL
void Triangle::setProblemDefinition(BORG_Problem &problem)
{
    // The parameter bounds are the same for all formulations
    BORG_Problem_set_bounds(problem, 0, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 1, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 2, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 3, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 4, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 5, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 6, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 7, 0.1, 0.47);
    BORG_Problem_set_bounds(problem, 8, 0.0, 0.37);
    BORG_Problem_set_bounds(problem, 9, 0.05, 0.42);
    BORG_Problem_set_bounds(problem, 10, 0.355, 0.725);
    BORG_Problem_set_bounds(problem, 11, 0.0, 0.1);
    BORG_Problem_set_bounds(problem, 12, 0.0, 0.1);
    BORG_Problem_set_bounds(problem, 13, 0.0, 0.1);
    BORG_Problem_set_bounds(problem, 14, 0.0, 0.1);
    BORG_Problem_set_bounds(problem, 15, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 16, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 17, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 18, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 19, 0.0, 0.02);
    BORG_Problem_set_bounds(problem, 20, 0.0, 0.02);
    BORG_Problem_set_bounds(problem, 21, 0.0, 0.02);
    BORG_Problem_set_bounds(problem, 22, 0.0, 0.02);
    BORG_Problem_set_bounds(problem, 23, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 24, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 25, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 26, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 27, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 28, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 29, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 30, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 31, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 32, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 33, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 34, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 35, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 36, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 37, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 38, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 39, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 40, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 41, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 42, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 43, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 44, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 45, 0.0, 1.0);
    BORG_Problem_set_bounds(problem, 46, 0.0, 100.0);
    BORG_Problem_set_bounds(problem, 47, 0.0, 100.0);
    BORG_Problem_set_bounds(problem, 48, 0.0, 100.0);
    BORG_Problem_set_bounds(problem, 49, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 50, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 51, 0.001, 1.0);
    BORG_Problem_set_bounds(problem, 52, 0.0, 10000.0);
    BORG_Problem_set_bounds(problem, 53, 0.0, 20.0);
    BORG_Problem_set_bounds(problem, 54, 0.0, 20.0);
    BORG_Problem_set_bounds(problem, 55, 0.0, 20.0);
    BORG_Problem_set_bounds(problem, 56, 0.0, 20.0);

    // Set epsilons for objectives
    BORG_Problem_set_epsilon(problem, 0, 0.002);
    BORG_Problem_set_epsilon(problem, 1, 0.02);
    BORG_Problem_set_epsilon(problem, 2, 10.);
    BORG_Problem_set_epsilon(problem, 3, 0.02);
    BORG_Problem_set_epsilon(problem, 4, 0.01);
    BORG_Problem_set_epsilon(problem, 5, 0.025);
}
#endif

Triangle::Triangle(unsigned long n_weeks, int n_realizations,
                                   int import_export_rof_table,
                                   string system_io,
                                   string &rof_tables_directory,
                                   int seed, unsigned long n_threads,
                                   string bootstrap_file,
                                   string &utilities_rdm_file,
                                   string &policies_rdm_file,
                                   string &water_sources_rdm_file,
                                   int n_sets, int n_bs_samples,
                                   string &solutions_file,
                                   vector<int> &solutions_to_run_range,
                                   bool plotting, bool print_obj_row)
        : HardCodedProblem(n_weeks, n_realizations, import_export_rof_table,
                           system_io, rof_tables_directory, seed, n_threads,
                           bootstrap_file, utilities_rdm_file,
                           policies_rdm_file,
                           water_sources_rdm_file, n_sets, n_bs_samples,
                           solutions_file, solutions_to_run_range, plotting,
                           print_obj_row) {
    readInputData();
}

Triangle::~Triangle() = default;


/**
 * Runs carolina problem.
 * @param vars
 * @param n_realizations
 * @param n_weeks
 * @param sol_number
 * @param output_directory
 * @todo check for solutions in which a utility does not have an allocation
 * on Jordan Lake (or any generic lake) but still pay for joint treatment
 * infrastructure).
 */
int Triangle::functionEvaluation(double *vars, double *objs, double *consts) {

    // ===================== SET UP DECISION VARIABLES  =====================

    Simulation *s = nullptr;
//    try {
    //throw invalid_argument("Test error");
    double Durham_restriction_trigger = vars[0];
    double OWASA_restriction_trigger = vars[1];
    double raleigh_restriction_trigger = vars[2];
    double cary_restriction_trigger = vars[3];
    double durham_transfer_trigger = vars[4];
    double owasa_transfer_trigger = vars[5];
    double raleigh_transfer_trigger = vars[6];
    double OWASA_JLA = vars[7];
    double Durham_JLA = vars[8];
    double Cary_JLA = vars[9];
    double Raleigh_JLA = vars[10];
    double durham_annual_payment = vars[11]; // contingency fund
    double owasa_annual_payment = vars[12];
    double raleigh_annual_payment = vars[13];
    double cary_annual_payment = vars[14];
    double durham_insurance_use = vars[15]; // insurance st_rof
    double owasa_insurance_use = vars[16];
    double raleigh_insurance_use = vars[17];
    double cary_insurance_use = vars[18];
    double durham_insurance_payment = vars[19];
    double owasa_insurance_payment = vars[20];
    double raleigh_insurance_payment = vars[21];
    double cary_insurance_payment = vars[22];



    /// Normalize Jordan Lake Allocations in case they exceed 1.
    double sum_jla_allocations = OWASA_JLA + Durham_JLA + Cary_JLA +
                                 Raleigh_JLA;
    if (sum_jla_allocations == 0.)
        throw invalid_argument("JLA allocations cannot be all "
                               "zero.");
    if (sum_jla_allocations > 0.69) { // At the time this study was done, 31% of JL had been allocated to other utilities.
        OWASA_JLA /= sum_jla_allocations / 0.69;
        Durham_JLA /= sum_jla_allocations / 0.69;
        Cary_JLA /= sum_jla_allocations / 0.69;
        Raleigh_JLA /= sum_jla_allocations / 0.69;
    }


       // ==================== SET UP RDM FACTORS ============================

    if (utilities_rdm.empty()) {
        /// All matrices below have dimensions n_realizations x nr_rdm_factors
        utilities_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(8, 1.));
        water_sources_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(55, 1.));
        policies_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(15, 1.));
    }


    // ===================== SET UP PROBLEM COMPONENTS =====================

    //cout << "BEGINNING TRIANGLE TEST" << endl << endl;
    // cout << "Using " << omp_get_num_threads() << " cores." << endl;
    //    cout << getexepath() << endl;

    /// Read streamflows
    int streamflow_n_weeks = (int) streamflows_durham[0].size();

    /// In case a vector containing realizations numbers to be calculated is passed, set
    /// number of realizations to number of realizations in that vector.

    //    vector<double> sewageFractions = Utils::parse1DCsvFile(
    //            io_directory + "/TestFiles/sewageFractions.csv");

    EvaporationSeries evaporation_durham(evap_durham, streamflow_n_weeks);
    EvaporationSeries evaporation_jordan_lake(
            evap_jordan_lake,
            streamflow_n_weeks);
    EvaporationSeries evaporation_falls_lake(
            evap_falls_lake,
            streamflow_n_weeks);
    EvaporationSeries evaporation_owasa(evap_owasa, streamflow_n_weeks);
    EvaporationSeries evaporation_little_river(
            evap_little_river,
            streamflow_n_weeks);
    EvaporationSeries evaporation_wheeler_benson(
            evap_wheeler_benson,
            streamflow_n_weeks);

    /// Create catchments and corresponding vectors
    //  Durham (Upper Neuse River Basin)
    Catchment durham_inflows(streamflows_durham, streamflow_n_weeks);

    //  Raleigh (Lower Neuse River Basin)
    Catchment lower_flat_river(streamflows_flat, streamflow_n_weeks);
    Catchment swift_creek(streamflows_swift, streamflow_n_weeks);
    Catchment little_river_raleigh(streamflows_llr, streamflow_n_weeks);
    Catchment crabtree_creek(streamflows_crabtree, streamflow_n_weeks);

    // OWASA (Upper Cape Fear Basin)
    Catchment phils_reek(streamflows_phils, streamflow_n_weeks);
    Catchment cane_creek(streamflows_cane, streamflow_n_weeks);
    Catchment morgan_creek(streamflows_morgan, streamflow_n_weeks);

    // Cary (Lower Cape Fear Basin)
    Catchment lower_haw_river(streamflows_haw, streamflow_n_weeks);

    // Downstream Gages
    Catchment neuse_river_at_clayton(streamflows_clayton, streamflow_n_weeks);
    Catchment cape_fear_river_at_lillington(
            streamflows_lillington,
            streamflow_n_weeks);

    vector<Catchment *> catchment_durham;

    vector<Catchment *> catchment_flat;
    vector<Catchment *> catchment_swift;
    vector<Catchment *> catchment_little_river_raleigh;
    vector<Catchment *> catchment_crabtree;

    vector<Catchment *> catchment_phils;
    vector<Catchment *> catchment_cane;
    vector<Catchment *> catchment_morgan;

    vector<Catchment *> catchment_haw;

    vector<Catchment *> gage_clayton;
    vector<Catchment *> gage_lillington;

    catchment_durham.push_back(&durham_inflows);

    catchment_flat.push_back(&lower_flat_river);
    catchment_swift.push_back(&swift_creek);
    catchment_little_river_raleigh.push_back(&little_river_raleigh);
    catchment_crabtree.push_back(&crabtree_creek);

    catchment_phils.push_back(&phils_reek);
    catchment_cane.push_back(&cane_creek);
    catchment_morgan.push_back(&morgan_creek);

    catchment_haw.push_back(&lower_haw_river);

    gage_clayton.push_back(&neuse_river_at_clayton);
    gage_lillington.push_back(&cape_fear_river_at_lillington);

    /// Storage vs. area reservoir curves.
    vector<double> falls_lake_storage = {0, 23266 * table_gen_storage_multiplier, 34700 * table_gen_storage_multiplier};
    vector<double> falls_lake_area = {0.32 * 5734, 0.32 * 29000, 0.28 * 40434};
    vector<double> wheeler_benson_storage = {0, 2789.66 * table_gen_storage_multiplier};
    vector<double> wheeler_benson_area = {0, 0.3675 * 2789.66};
    vector<double> teer_storage = {0, 1315.0};
    vector<double> teer_area = {20, 50};
    vector<double> little_river_res_storage = {0, 3700};
    vector<double> little_river_res_area = {0, 0.3675 * 3700};

    DataSeries falls_lake_storage_area(falls_lake_storage, falls_lake_area);
    DataSeries wheeler_benson_storage_area(wheeler_benson_storage,
                                           wheeler_benson_area);
    DataSeries teer_storage_area(teer_storage,
                                 teer_area);
    DataSeries little_river_storage_area(little_river_res_storage,
                                         little_river_res_area);

    /// Minimum environmental flow rules (controls)
    vector<int> dlr_weeks = {0, 21, 47, 53};
    vector<double> dlr_releases = {3.877 * 7, 9.05, 3.877 * 7};
    vector<double> wb_storage = {0.3 * 2789.66, 0.6 * 2789.66, 2789.66};
    vector<double> wb_releases = {0.646 * 7, 1.29 * 7, 1.94 * 7};
    vector<double> ccr_inflows = {0.1422 * 7, 0.5 * 7, 1 * 7, 1.5 * 7,
                                  1.797 * 7};
    vector<double> ccr_releases = {0.1422 * 7, 0.5 * 7, 1 * 7, 1.5 * 7,
                                   1.797 * 7};
    vector<int> falls_controls_weeks = {13, 43};
    vector<double> falls_base_releases = {64.64 * 7, 38.78 * 7};
    vector<double> falls_min_gage = {180 * 7, 130 * 7};

    SeasonalMinEnvFlowControl durham_min_env_control(0, dlr_weeks,
                                                     dlr_releases);
    //    FixedMinEnvFlowControl falls_min_env_control(1, 38.78 * 7);
    FallsLakeMinEnvFlowControl falls_min_env_control(1,
                                                     10,
                                                     falls_controls_weeks,
                                                     falls_base_releases,
                                                     falls_min_gage,
                                                     crabtree_creek);

    StorageMinEnvFlowControl wheeler_benson_min_env_control(2,
                                                            vector<int>(1, 2),
                                                            wb_storage,
                                                            wb_releases);
    FixedMinEnvFlowControl sq_min_env_control(3, 0);
    InflowMinEnvFlowControl ccr_min_env_control(4,
                                                ccr_inflows,
                                                ccr_releases);
    FixedMinEnvFlowControl university_min_env_control(5, 0);
    //    FixedMinEnvFlowControl jordan_min_env_control(6,
    //                                                  25.8527 * 7);
    JordanLakeMinEnvFlowControl jordan_min_env_control(
            6, cape_fear_river_at_lillington, 64.63, 129.26, 25.85, 193.89,
            290.84, 387.79, 30825.0, 14924.0);


    //    vector<int> eno_weeks = {7, 16, 53};
    //    vector<double> eno_releases = {6.49 * 7, 19.48 * 7, 6.49 * 7};
    //    SeasonalMinEnvFlowControl eno_min_env_control(&eno_weeks, &eno_releases);

    vector<MinEnvFlowControl *> min_env_flow_controls;
    min_env_flow_controls.push_back(&durham_min_env_control);
    min_env_flow_controls.push_back(&falls_min_env_control);
    min_env_flow_controls.push_back(&wheeler_benson_min_env_control);
    min_env_flow_controls.push_back(&sq_min_env_control);
    min_env_flow_controls.push_back(&ccr_min_env_control);
    min_env_flow_controls.push_back(&university_min_env_control);
    min_env_flow_controls.push_back(&jordan_min_env_control);

    double discount_rate = 0.05;

    /// Jordan Lake parameters
    double jl_supply_capacity = 14924.0 * table_gen_storage_multiplier;
    double jl_wq_capacity = 30825.0 * table_gen_storage_multiplier;
    double jl_storage_capacity = jl_wq_capacity + jl_supply_capacity;
    vector<int> jl_allocations_ids = {0, 1, 2, 3, WATER_QUALITY_ALLOCATION};
    vector<double> jl_allocation_fractions = {
            OWASA_JLA * jl_supply_capacity / jl_storage_capacity,
            Durham_JLA * jl_supply_capacity / jl_storage_capacity,
            Cary_JLA * jl_supply_capacity / jl_storage_capacity,
            Raleigh_JLA * jl_supply_capacity / jl_storage_capacity,
            jl_wq_capacity / jl_storage_capacity};
    vector<double> jl_treatment_allocation_fractions = {0.0, 0.0, 1.0, 0.0};

    /// Jordan Lake parameters
    double fl_supply_capacity = 14700.0 * table_gen_storage_multiplier;
    double fl_wq_capacity = 20000.0 * table_gen_storage_multiplier;
    double fl_storage_capacity = fl_wq_capacity + fl_supply_capacity;
    vector<int> fl_allocations_ids = {3, WATER_QUALITY_ALLOCATION};
    vector<double> fl_allocation_fractions = {
            fl_supply_capacity / fl_storage_capacity,
            fl_wq_capacity / fl_storage_capacity};
//        vector<double> fl_treatment_allocation_fractions = {0.0, 0.0, 0.0, 1.0};
    vector<double> fl_treatment_allocation_fractions = {1.0};

    // Existing Sources
    Reservoir durham_reservoirs("Lake Michie & Little River Res. (Durham)",
                                0,
                                catchment_durham,
                                6349.0 * table_gen_storage_multiplier,
                                ILLIMITED_TREATMENT_CAPACITY,
                                evaporation_durham, 1069);
    //    Reservoir falls_lake("Falls Lake", 1, catchment_flat,
    //                         34700.0, 99999,
    //                         &evaporation_falls_lake, &falls_lake_storage_area);
    AllocatedReservoir falls_lake("Falls Lake",
                                  1,
                                  catchment_flat,
                                  fl_storage_capacity,
                                  ILLIMITED_TREATMENT_CAPACITY,
                                  evaporation_falls_lake,
                                  falls_lake_storage_area,
                                  fl_allocations_ids,
                                  fl_allocation_fractions,
                                  fl_treatment_allocation_fractions);

    Reservoir wheeler_benson_lakes("Wheeler-Benson Lakes", 2, catchment_swift,
                                   2789.66 * table_gen_storage_multiplier,
                                   ILLIMITED_TREATMENT_CAPACITY,
                                   evaporation_wheeler_benson,
                                   wheeler_benson_storage_area);
    Reservoir stone_quarry("Stone Quarry",
                           3,
                           catchment_phils,
                           200.0 * table_gen_storage_multiplier,
                           ILLIMITED_TREATMENT_CAPACITY,
                           evaporation_owasa,
                           10);
    Reservoir ccr("Cane Creek Reservoir",
                  4,
                  catchment_cane,
                  2909.0 * table_gen_storage_multiplier,
                  ILLIMITED_TREATMENT_CAPACITY,
                  evaporation_owasa,
                  500);
    Reservoir university_lake("University Lake", 5, catchment_morgan,
                              449.0 * table_gen_storage_multiplier,
                              ILLIMITED_TREATMENT_CAPACITY,
                              evaporation_owasa,
                              212);
    AllocatedReservoir jordan_lake("Jordan Lake",
                                   6,
                                   catchment_haw,
                                   jl_storage_capacity,
                                   448,
                                   evaporation_jordan_lake,
                                   13940,
                                   jl_allocations_ids,
                                   jl_allocation_fractions,
                                   jl_treatment_allocation_fractions);

    // other than Cary WTP for Jordan Lake, assume no WTP constraints - each
    // city can meet its daily demands with available treatment infrastructure

    vector<double> construction_time_interval = {3.0, 5.0};
    LevelDebtServiceBond dummy_bond(7, 1., 1, 1., vector<int>(1, 0));
    Reservoir dummy_endpoint("Dummy Node", 7, vector<Catchment *>(), 1., 0, evaporation_durham, 1, {},
                             construction_time_interval, 0, dummy_bond);

    vector<WaterSource *> water_sources;
    water_sources.push_back(&durham_reservoirs);
    water_sources.push_back(&falls_lake);
    water_sources.push_back(&wheeler_benson_lakes);
    water_sources.push_back(&stone_quarry);
    water_sources.push_back(&ccr);
    water_sources.push_back(&university_lake);
    water_sources.push_back(&jordan_lake);
    water_sources.push_back(&dummy_endpoint);


    /*
     * System connection diagram (water
     * flows from top to bottom)
     * Potential projects and expansions
     * of existing sources in parentheses
     *
     *        3          4       5                 0
     *         \         /      /                  |
     *          \       /      /                   |
     *           \     /      /                    |
     *           |    /      /                     |
     *           |   /      /                      |
     *           |   |     /                        \
     *           |   |    /                           1                 2
     *           |   |   /                             |                |
     *            \__/__/                              |                |
     *              6                                  |                |
     *              |                                   \               |
     *              |                                    \              |
     *        Lillington Gage                             \             |
     *              |                                      |            |
     *              |                                      |           /
     *              |                                      |          /
     *              |                                 Clayton Gage   /
     *              |                                      |        /
     *               \                                     |   -----
     *                \                                     \ /
     *                 \                                     |
     *                  \                                    |
     *                   \                                   |
     *                    \                                   \
     *                     -------                            /
     *                            \             --------------
     *                             \           /
     *                              \     -----
     *                               \   /
     *                                \ /
     *                                 7
     */

    Graph g(8);
    g.addEdge(0, 1);
    g.addEdge(1, 7);
    g.addEdge(2, 7);

    g.addEdge(3, 6);
    g.addEdge(4, 6);
    g.addEdge(5, 6);
    g.addEdge(6, 7);

    // changing for 5 year sim period
    auto demand_n_weeks = (int) round(5 * WEEKS_IN_YEAR);

    vector<int> cary_ws_return_id;
    vector<vector<double>> cary_discharge_fraction_series;
    WwtpDischargeRule wwtp_discharge_cary(
            cary_discharge_fraction_series,
            cary_ws_return_id);
    vector<int> owasa_ws_return_id = {6};
    WwtpDischargeRule wwtp_discharge_owasa(
            demand_to_wastewater_fraction_owasa_raleigh,
            owasa_ws_return_id);
    vector<int> raleigh_ws_return_id = {7};
    WwtpDischargeRule wwtp_discharge_raleigh(
            demand_to_wastewater_fraction_owasa_raleigh,
            raleigh_ws_return_id);
    vector<int> durham_ws_return_id = {1, 6};
    WwtpDischargeRule wwtp_discharge_durham(
            demand_to_wastewater_fraction_durham,
            durham_ws_return_id);

    Utility owasa("OWASA", 0, demand_owasa, demand_n_weeks, owasa_annual_payment, owasaDemandClassesFractions,
                  owasaUserClassesWaterPrices, wwtp_discharge_owasa, 1.0);

    Utility durham("Durham", 1, demand_durham, demand_n_weeks, durham_annual_payment, durhamDemandClassesFractions,
                   durhamUserClassesWaterPrices, wwtp_discharge_durham, 1.0);

    Utility cary("Cary", 2, demand_cary, demand_n_weeks, cary_annual_payment, caryDemandClassesFractions,
            caryUserClassesWaterPrices, wwtp_discharge_cary, 1.0);

    Utility raleigh("Raleigh", 3, demand_raleigh, demand_n_weeks, raleigh_annual_payment, raleighDemandClassesFractions,
            raleighUserClassesWaterPrices, wwtp_discharge_raleigh, 1.0);


    vector<Utility *> utilities;
    utilities.push_back(&cary);
    utilities.push_back(&durham);
    utilities.push_back(&owasa);
    utilities.push_back(&raleigh);

    /// Water-source-utility connectivity matrix (each row corresponds to a utility and numbers are water
    /// sources IDs.
    vector<vector<int>> reservoir_utility_connectivity_matrix = {
            {3, 4,  5, 6},  //OWASA
            {0, 6},         //Durham
            {6},            //Cary
            {1, 2, 6}       //Raleigh
    };

    auto table_storage_shift = vector<vector<double>>(4, vector<double>(25, 0.));
    //table_storage_shift[3][17] = 2000.;
    //table_storage_shift[3][8] = 5000.;
    //table_storage_shift[1][14] = 100.;
    //table_storage_shift[1][20] = 500.;
    //table_storage_shift[1][21] = 500.;
    //table_storage_shift[1][15] = 700.;
    //table_storage_shift[1][9] = 700.;

    vector<DroughtMitigationPolicy *> drought_mitigation_policies;
    /// Restriction policies
    vector<double> initial_restriction_triggers = {OWASA_restriction_trigger,
                                                   Durham_restriction_trigger,
                                                   cary_restriction_trigger,
                                                   raleigh_restriction_trigger};

    vector<double> restriction_stage_multipliers_cary = {0.9, 0.8, 0.7, 0.6};
    vector<double> restriction_stage_triggers_cary = {initial_restriction_triggers[0],
                                                      initial_restriction_triggers[0] + 0.15f,
                                                      initial_restriction_triggers[0] + 0.35f,
                                                      initial_restriction_triggers[0] + 0.6f};
    vector<double> restriction_stage_multipliers_durham = {0.9, 0.8, 0.7, 0.6};
    vector<double> restriction_stage_triggers_durham = {initial_restriction_triggers[1],
                                                        initial_restriction_triggers[1] + 0.15f,
                                                        initial_restriction_triggers[1] + 0.35f,
                                                        initial_restriction_triggers[1] + 0.6f};
    vector<double> restriction_stage_multipliers_owasa = {0.9, 0.8, 0.7};
    vector<double> restriction_stage_triggers_owasa = {initial_restriction_triggers[2],
                                                       initial_restriction_triggers[2] + 0.15f,
                                                       initial_restriction_triggers[2] + 0.35f};
    vector<double> restriction_stage_multipliers_raleigh = {0.9, 0.8, 0.7, 0.6};
    vector<double> restriction_stage_triggers_raleigh = {initial_restriction_triggers[3],
                                                         initial_restriction_triggers[3] + 0.15f,
                                                         initial_restriction_triggers[3] + 0.35f,
                                                         initial_restriction_triggers[3] + 0.6f};

    Restrictions restrictions_o(0,
                                restriction_stage_multipliers_owasa,
                                restriction_stage_triggers_owasa,
                                &owasaDemandClassesFractions,
                                &owasaUserClassesWaterPrices,
                                &owasaPriceSurcharges);

    Restrictions restrictions_d(1,
                                restriction_stage_multipliers_durham,
                                restriction_stage_triggers_durham);

    Restrictions restrictions_c(2,
                                restriction_stage_multipliers_cary,
                                restriction_stage_triggers_cary);

    Restrictions restrictions_r(3,
                                restriction_stage_multipliers_raleigh,
                                restriction_stage_triggers_raleigh);

    drought_mitigation_policies = {&restrictions_o,  &restrictions_d, &restrictions_c,
                                    &restrictions_r};

    /// Transfer policy
    /*
     *      2
     *     / \
     *  0 v   v 1
     *   /     \
     *  3---><--1--><--0
     *      2       3
     */

    vector<int> buyers_ids = {0, 1, 3};
    //FIXME: CHECK IF TRANSFER CAPACITIES MATCH IDS IN BUYERS_IDS.
    vector<double> buyers_transfers_capacities = {10.8 * 7, 10.0 * 7, 11.5 * 7,
                                                  7.0 * 7};
    vector<double> buyers_transfers_trigger = {owasa_transfer_trigger,
                                               durham_transfer_trigger,
                                               raleigh_transfer_trigger};

    Graph ug(4);
    ug.addEdge(2, 1);
    ug.addEdge(2, 3);
    ug.addEdge(1, 3);
    ug.addEdge(1, 0);

    Transfers t(4, 2, 6, 35,
                buyers_ids,
                buyers_transfers_capacities,
                buyers_transfers_trigger,
                ug,
                vector<double>(),
                vector<int>());
    drought_mitigation_policies.push_back(&t);

    vector<double> insurance_triggers = {owasa_insurance_use,
                                         durham_insurance_use, cary_insurance_use,
                                         raleigh_insurance_use}; //FIXME: Change per solution
    vector<double> fixed_payouts = {owasa_insurance_payment,
                                    durham_insurance_payment,
                                    cary_insurance_payment,
                                    raleigh_insurance_payment};
    vector<int> insured_utilities = {0, 1, 2, 3};
    //double insurance_premium = 1.2;
    //InsuranceStorageToROF in(5, water_sources, g, reservoir_utility_connectivity_matrix, utilities,
    //                         min_env_flow_controls, utilities_rdm, water_sources_rdm, insurance_triggers,
    //                         insurance_premium, fixed_payouts, n_weeks);

    //drought_mitigation_policies.push_back(&in);

    /// Creates simulation object depending on use (or lack thereof) ROF tables
    double start_time = omp_get_wtime();
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        s = new Simulation(water_sources,
                           g,
                           reservoir_utility_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           utilities_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run,
                           rof_tables_directory);
        //double realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    } else if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        s = new Simulation (water_sources,
                            g,
                            reservoir_utility_connectivity_matrix,
                            utilities,
                            drought_mitigation_policies,
                            min_env_flow_controls,
                            utilities_rdm,
                            water_sources_rdm,
                            policies_rdm,
                            n_weeks,
                            realizations_to_run,
                            rof_tables,
                            table_storage_shift,
                            rof_tables_directory);
       // double realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    } else {
        s = new Simulation(water_sources,
                           g,
                           reservoir_utility_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           utilities_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    }
    double end_time = omp_get_wtime();
	//printf("Function evaluation time: %f s\n", end_time - start_time);

    //double realization_end = omp_get_wtime();
    std::cout << std::endl;
    std::cout << "Simulation took  " << end_time - start_time << "s" << std::endl;

    /// Calculate objectives and store them in Borg decision variables array.
#ifdef  PARALLEL
    objectives = calculateAndPrintObjectives(false);

        int i = 0;
        objs[i] = min(min(objectives[i], objectives[4 + i]),
        		   min(objectives[8 + i], objectives[12 + i])) - 1.;
        for (i = 1; i < 4; ++i) {
            objs[i] = max(max(objectives[i], objectives[4 + i]),
      	                  max(objectives[8 + i], objectives[12 + i]));
        }

        objs[4] = OWASA_JLA + Durham_JLA + Cary_JLA + Raleigh_JLA;

        objectives.push_back(objs[4]);

        if (s != nullptr) {
            delete s;
	}
	s = nullptr;
#endif
//    } catch (const std::exception& e) {
//        simulationExceptionHander(e, s, objs, vars);
//	return 1;
//    }

    delete s;

    return 0;
}

int Triangle::simulationExceptionHander(const std::exception &e, Simulation *s,
                                        double *objs, const double *vars) {
    int num_dec_var = 56;
//        printf("Exception called during calculations. Decision variables are below:\n");
    ofstream sol;
    int world_rank;

#ifdef  PARALLEL
    // int mpi_initialized;
	// MPI_Initialized(&mpi_initialized);
	// if (mpi_initialized)
 //            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// else
	    world_rank = 0;
#else
    world_rank = 0;
#endif
    string error_file = "sol_error_rank_" + to_string(world_rank) + ".csv";
    sol.open(error_file.c_str());
    for (int i = 0; i < num_dec_var; ++i) {
        sol << vars[i] << ",";
    }
    sol << flush;
    sol.close();
    printf("Error. Decision variables printed in %s\n", error_file.c_str());

#ifdef PARALLEL
    objs[0] = 0.;
	objs[1] = 1.1;
	objs[2] = 1000;
	objs[3] = 5.;
	objs[4] = 5.;
	objs[5] = 1.1;
	if (s != nullptr) {
	    delete s;
	    s = nullptr;
	}
#else
    Utils::print_exception(e);
#endif

    return 1;
}




/*
Triangle::Triangle(unsigned long n_weeks, int import_export_rof_table)
        : Problem(n_weeks) {
    if (import_export_rof_table == EXPORT_ROF_TABLES) {
        table_gen_storage_multiplier = BASE_STORAGE_CAPACITY_MULTIPLIER;
    } else {
        table_gen_storage_multiplier = 1.;
    }
}
*/


void Triangle::readInputData() {
    cout << "Reading input data." << endl;
    string data_dir = DEFAULT_DATA_DIR + BAR;

    string rdm_tseries_dir;
    if (rdm_no != NON_INITIALIZED) {
        rdm_tseries_dir = "rdm_tseries/rdm_" + to_string(rdm_no) + BAR;
    } else {
        rdm_tseries_dir = DEFAULT_RDM_TSERIES_DIR;
    }

#pragma omp parallel num_threads(omp_get_thread_num())
    {
#pragma omp single
        streamflows_durham = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "durham_inflows.csv", n_realizations);
#pragma omp single
        streamflows_flat = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "falls_lake_inflows.csv", n_realizations);
#pragma omp single
        streamflows_swift = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "lake_wb_inflows.csv", n_realizations);
#pragma omp single
        streamflows_llr = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "little_river_raleigh_inflows.csv", n_realizations);
        // }
#pragma omp single
        streamflows_crabtree = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "crabtree_inflows.csv", n_realizations);
#pragma omp single
        streamflows_phils = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "stone_quarry_inflows.csv", n_realizations);
#pragma omp single
        streamflows_cane = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "cane_creek_inflows.csv", n_realizations);
#pragma omp single
        streamflows_morgan = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "university_lake_inflows.csv", n_realizations);
#pragma omp single
        streamflows_haw = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "jordan_lake_inflows.csv", n_realizations);
#pragma omp single
        streamflows_clayton = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "clayton_inflows.csv", n_realizations);
#pragma omp single
        streamflows_lillington = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "lillington_inflows.csv", n_realizations);
// };
        //cout << "Reading evaporations." << endl;
#pragma omp single
        evap_durham = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "durham_evap.csv", n_realizations);
#pragma omp single
        evap_falls_lake = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "falls_lake_evap.csv", n_realizations);
#pragma omp single
        evap_owasa = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "owasa_evap.csv", n_realizations);
#pragma omp single
        evap_little_river = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "little_river_raleigh_evap.csv", n_realizations);
#pragma omp single
        {
            evap_wheeler_benson = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                    BAR + "wb_evap.csv", n_realizations);
            evap_jordan_lake = evap_owasa;
        }

        //cout << "Reading demands." << endl;
#pragma omp single
        demand_cary = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" + evap_inflows_suffix +
                BAR + "cary_demand.csv", n_realizations);
#pragma omp single
        demand_durham = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" + evap_inflows_suffix +
                BAR + "durham_demand.csv", n_realizations);
#pragma omp single
        demand_raleigh = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" + evap_inflows_suffix +
                BAR + "raleigh_demand.csv", n_realizations);
#pragma omp single
        demand_owasa = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" + evap_inflows_suffix +
                BAR + "owasa_demand.csv", n_realizations);

        //cout << "Reading others." << endl;
#pragma omp single
        {
            demand_to_wastewater_fraction_owasa_raleigh = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "demand_to_wastewater_fraction_owasa_raleigh.csv");
            demand_to_wastewater_fraction_durham = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "demand_to_wastewater_fraction_durham.csv");

            caryDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "caryDemandClassesFractions.csv");
            durhamDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "durhamDemandClassesFractions.csv");
            raleighDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "raleighDemandClassesFractions.csv");
            owasaDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "owasaDemandClassesFractions.csv");

            caryUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "caryUserClassesWaterPrices.csv");
            durhamUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "durhamUserClassesWaterPrices.csv");
            raleighUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "raleighUserClassesWaterPrices.csv");
            owasaUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "owasaUserClassesWaterPrices.csv");

            owasaPriceSurcharges = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "owasaPriceRestMultipliers.csv");
        }
//    cout << "Done reading input data." << endl;
    }

}
