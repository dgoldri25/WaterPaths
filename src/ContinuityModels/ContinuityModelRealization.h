//
// Created by bernardo on 1/26/17.
//

#ifndef TRIANGLEMODEL_CONTINUITYMODELREALIZATION_H
#define TRIANGLEMODEL_CONTINUITYMODELREALIZATION_H


#include "Base/ContinuityModel.h"
#include "../DroughtMitigationInstruments/Base/DroughtMitigationPolicy.h"
#include "../Controls/AllocationModifier.h"

class ContinuityModelRealization : public ContinuityModel {
private:
    vector<DroughtMitigationPolicy *> drought_mitigation_policies;

public:
    ContinuityModelRealization(
            vector<WaterSource *> &water_sources,
            const Graph &water_sources_graph,
            const vector<vector<int>> &water_sources_to_utilities,
            vector<Utility *> &utilities,
            const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
            vector<MinEnvironFlowControl *> &min_env_flow_control,
            vector<vector<double>> *utilities_rdm,
            vector<vector<double>> *water_sources_rdm,
            vector<vector<double>> *policy_rdm,
            const unsigned int realization_index);

    ContinuityModelRealization(
            ContinuityModelRealization &continuity_model_realization);

    vector<WaterSource *> getWater_sources();

    void setShortTermROFs(const vector<double> &risks_of_failure);

    void applyDroughtMitigationPolicies(int week);

    const vector<DroughtMitigationPolicy *> getDrought_mitigation_policies() const;

    void setLongTermROFs(const vector<double> &risks_of_failure, const int week);

    virtual ~ContinuityModelRealization();
};


#endif //TRIANGLEMODEL_CONTINUITYMODELREALIZATION_H
