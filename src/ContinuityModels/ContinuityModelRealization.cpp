//
// Created by bernardo on 1/26/17.
//

#include <iostream>
#include <algorithm>
#include "ContinuityModelRealization.h"

ContinuityModelRealization::ContinuityModelRealization(
        vector<WaterSource *> &water_sources,
        const Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_utilities,
        vector<Utility *> &utilities,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_control,
        const vector<double>& utilities_rdm,
        const vector<double>& water_sources_rdm,
        const vector<double>& policy_rdm,
        const unsigned int realization_id)
        : ContinuityModel(water_sources, utilities, min_env_flow_control, water_sources_graph,
                          water_sources_to_utilities, utilities_rdm, water_sources_rdm,
                          realization_id),
          drought_mitigation_policies(drought_mitigation_policies) {

    // Pass corresponding utilities to drought mitigation instruments.
    for (DroughtMitigationPolicy *dmp : this->drought_mitigation_policies) {
        dmp->addSystemComponents(utilities, water_sources, min_env_flow_control);
        dmp->setRealization(realization_id, utilities_rdm, water_sources_rdm,
                            policy_rdm);
    }

    // offset the jordan lake allocations

    vector<double>  realization_JLA_allocations (n_utilities+1);

    double adjusted_fraction;
    double jl_storage_capacity = 14924.0 + 30825.0;
    //double base_allocated_fraction;
    double DV_allocation;
    double sum_adjustments = 0;
    double jl_supply_capacity = this->continuity_water_sources[6]->getSupplyCapacity();
    vector<double> adjusted_supply_allocations (n_utilities);


    // apply offsets
/*
    for (int JLA_id = 0; JLA_id < this->n_utilities; JLA_id++){
        DV_allocation = this->continuity_water_sources[6]->getSupplyAllocatedFraction(JLA_id);

        // convert to original decision variable
        //DV_allocation = (base_allocated_fraction*jl_storage_capacity) /
        //        (this->continuity_water_sources[6]->getSupplyCapacity());

        // add offset from rdm input
        DV_allocation += water_sources_rdm[JLA_id+51];


        if (DV_allocation < 0){
            DV_allocation = 0;
        }

        // convert back to fraction of JLA
        adjusted_fraction = (DV_allocation * jl_supply_capacity) / jl_storage_capacity;
        sum_adjustments += (water_sources_rdm[JLA_id+51]*this->continuity_water_sources[6]->getSupplyCapacity()) /
                           jl_storage_capacity;

        // add to vector of adjusted allocations
        realization_JLA_allocations[JLA_id] = adjusted_fraction;
        adjusted_supply_allocations[JLA_id] = DV_allocation;
    }

    // add/remove the total sum of adjustments from the wq fraction
    realization_JLA_allocations[n_utilities] = this->continuity_water_sources[6]->getWqFraction() - sum_adjustments;

    // adjust the allocations using the resetAllocations function
    this->continuity_water_sources[6]->resetAllocations(&realization_JLA_allocations);
    //this->continuity_water_sources[6]->setSupplyAllocations(&adjusted_supply_allocations);
    */




}

ContinuityModelRealization::~ContinuityModelRealization() {
    // Delete drought mitigation policies.
    for (auto dmp : drought_mitigation_policies) {
        delete dmp;
    }
}

void ContinuityModelRealization::setShortTermROFs(const vector<double> &risks_of_failure) {
    for (unsigned long i = 0; i < continuity_utilities.size(); ++i) {
        continuity_utilities.at(i)->setRisk_of_failure(risks_of_failure.at(i));
    }
}

void ContinuityModelRealization::setLongTermROFs(const vector<double> &risks_of_failure, const int week) {
    vector<int> new_infra_triggered;
    int nit; // new infrastruction triggered - id.

    // Loop over utilities to see if any of them will build new infrastructure.
    for (unsigned long u = 0; u < continuity_utilities.size(); ++u) {
        // Runs utility's infrastructure construction handler and get the id
        // of new source built, if any.
        nit = continuity_utilities[u]->
                infrastructureConstructionHandler(risks_of_failure[u], week);
        // If new source was built, check add it to the list of sources
        // built by all utilities.
        if (nit != NON_INITIALIZED)
            new_infra_triggered.push_back(nit);
    }

    // Look for and remove duplicates, in the unlikely case two utilities
    // build the same source at the same time. This will prevent the source
    // from being erased from a utility which will later try to build it.
    sort(new_infra_triggered.begin(),
         new_infra_triggered.end());
    new_infra_triggered.erase(unique(new_infra_triggered.begin(),
                                     new_infra_triggered.end()),
                              new_infra_triggered.end());

    // If infrastructure was built, force utilities to build their share of
    // that infrastructure option (which will only happen it the listed
    // option is in the list of sources to be built for other utilities.
    if (!new_infra_triggered.empty())
        for (Utility *u : continuity_utilities) {
            u->forceInfrastructureConstruction(week, new_infra_triggered);
        }
}

void ContinuityModelRealization::applyDroughtMitigationPolicies(int week) {
    for (DroughtMitigationPolicy* dmp : drought_mitigation_policies) {
        dmp->applyPolicy(week);
    }
}

const vector<DroughtMitigationPolicy *> ContinuityModelRealization::getDrought_mitigation_policies() const {
    return drought_mitigation_policies;
}
