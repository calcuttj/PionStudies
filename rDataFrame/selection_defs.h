auto interaction_topology = [](int origin, int pdg, int ID, double endZ,
                               std::string process, int nPi0,
                               int beam_hit_true_origin,
                               int beam_hit_true_ID,
                               std::vector<int> & true_daughter_pdg,
                               std::vector<double> & true_daughter_startP) {
  int topology = -1;
  if (origin == 2) {
    topology = 7; 
  }
  else if (pdg == 211) {
    //if (!beam_hit_true_ID.size())
    if (endZ < 0.) {
      topology = 4;
    }
    else if (endZ > 226.) {
      topology = 7;
    }
    else if (process == "pi+Inelastic") {
      //daughters with & without thresholds
      bool has_pion_above_threshold = false;
      for (size_t i = 0; i < true_daughter_startP.size(); ++i) {
        if (abs(true_daughter_pdg[i]) == 211 &&
            true_daughter_startP[i] > .150) {
          has_pion_above_threshold = true;
          break;
        }
      }

      if (has_pion_above_threshold) {
        topology = 3;
      }
      else if (nPi0 == 0) {
        topology = 1;
      }
      else if (nPi0 > 0) {
        topology = 2;
      }
    }
    else {
      topology = 7;
    }
  }
  else if (pdg == -13) {
    if (beam_hit_true_origin == 2) {
      topology = 7; //Other for now
    }
    else if (beam_hit_true_ID == ID) {
      topology = 5;
    }
    else {
      topology = 7;
    }
  }
  else {
    topology = 7;
  }

  return topology;
};

auto new_interaction_topology = [](int pdg, double endZ,
                                   std::string process, int nPi0,
                                   std::vector<int> & true_daughter_pdg,
                                   std::vector<double> & true_daughter_startP,
                                   std::vector<double> & incidentEnergies) {

  int topology = -1;
  if (pdg == 211) {
    //Before FV
    //if (!incidentEnergies.size()) {
    if (endZ < -.49375) {
      topology = 4;
    }
    //After FV
    //else if (endZ > 225.94353) {
    else if (endZ > 222.10561) {
      topology = 6;
    }
    else if (process == "pi+Inelastic") {
      //daughters with & without thresholds
      bool has_pion_above_threshold = false;
      for (size_t i = 0; i < true_daughter_startP.size(); ++i) {
        if (abs(true_daughter_pdg[i]) == 211 &&
            true_daughter_startP[i] > .150) {
          has_pion_above_threshold = true;
          break;
        }
      }

      if (has_pion_above_threshold) {
        topology = 3;
      }
      else if (nPi0 == 0) {
        topology = 1;
      }
      else if (nPi0 > 0) {
        topology = 2;
      }
    }
    else {
      topology = 7;
    }
  }
  else if (pdg == -13) {
    topology = 5;
  }
  else {
    topology = 7;
  }

  return topology;

};

auto selection_ID = [](bool beam_is_track, bool ends_in_APA3,
                       bool no_pion_daughter,
                       bool beam_cuts, bool has_shower) {
  if (!beam_is_track) {
    return 6;
  }

  if (!beam_cuts) {
    return 5;
  }
  
  if (!ends_in_APA3) {
    return 4;
  }

  if (!no_pion_daughter) {
    return 3;
  }

  if (has_shower) {
    return 2;
  }
  else {
    return 1;
  }
};
