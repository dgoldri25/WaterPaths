//
// Created by Bernardo on 11/29/2019.
//

#ifndef TRIANGLEMODEL_DROUGHTINSURANCEPARSER_H
#define TRIANGLEMODEL_DROUGHTINSURANCEPARSER_H


#include "../Base/DroughtMitigationPolicyParser.h"

class DroughtInsuranceParser : public DroughtMitigationPolicyParser {
private:
    vector<double> insurance_triggers;
    double insurance_premium = NON_INITIALIZED;
    vector<double> fixed_payouts;
public:

    DroughtInsuranceParser();

    ~DroughtInsuranceParser() override;

    void parseVariables(vector<vector<string>> &block, int n_realizations,
                        int n_weeks, int line_no, Graph &utilities_graph,
                        Graph &ws_graph,
                        const map<string, int> &utility_name_to_id) override;

    void checkMissingOrExtraParams(int line_no,
                                   vector<vector<string>> &block) override;

    DroughtMitigationPolicy *
    generateSource(int id, vector<vector<string>> &block, int line_no,
                   int n_realizations, int n_weeks,
                   const map<string, int> &ws_name_to_id,
                   const map<string, int> &utility_name_to_id,
                   Graph &utilities_graph, Graph &ws_graph,
                   const vector<vector<int>> &water_sources_to_utilities,
                   vector<Utility *> &utilities,
                   vector<WaterSource *> &water_sources,
                   vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
                   vector<MinEnvFlowControl *> min_env_flow_controls,
                   vector<vector<double>> &utilities_rdm,
                   vector<vector<double>> &water_sources_rdm,
                   vector<vector<double>> &policy_rdm) override;

};


#endif //TRIANGLEMODEL_DROUGHTINSURANCEPARSER_H
