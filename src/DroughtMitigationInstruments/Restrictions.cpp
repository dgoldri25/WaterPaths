//
// Created by bernardo on 2/3/17.
//

#include "Restrictions.h"

Restrictions::Restrictions(const int id, const vector<double> &stage_multipliers,
                           const vector<double> &stage_triggers, const vector<int> utilities_ids)
        : DroughtMitigationPolicy(id, utilities_ids),
          stage_multipliers(stage_multipliers),
          stage_triggers(stage_triggers) {}

Restrictions::Restrictions(const Restrictions &restrictions) : DroughtMitigationPolicy(restrictions.id,
                                                                                       restrictions.utilities_ids),
                                                               stage_multipliers(restrictions.stage_multipliers),
                                                               stage_triggers(restrictions.stage_triggers),
                                                               utility(restrictions.utility) {}

Restrictions::~Restrictions() {}

void Restrictions::applyPolicy(int week) {

    current_multiplier = 1.0;
    for (int i = 0; i < stage_triggers.size(); ++i) {
        if (utility->getRisk_of_failure() > stage_triggers[i]) {
            current_multiplier = stage_multipliers[i];
        } else break;
    }

    utility->setDemand(week, current_multiplier);
}

double Restrictions::getCurrent_multiplier() const {
    return current_multiplier;
}

void Restrictions::addUtility(Utility *utility) {
    this->utility = utility;
}
