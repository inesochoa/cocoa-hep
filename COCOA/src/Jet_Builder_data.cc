#include "Jet_Builder_data.hh"

Jet_Builder_data::Jet_Builder_data(std::string prefix)
{
    Prefix = prefix;
}

Jet_Builder_data::~Jet_Builder_data()
{;}

void Jet_Builder_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch(TString(Prefix + "_jet_pt"), "vector<float>", &jet_pt);
    outTree->Branch(TString(Prefix + "_jet_eta"), "vector<float>", &jet_eta);
    outTree->Branch(TString(Prefix + "_jet_phi"), "vector<float>", &jet_phi);
    outTree->Branch(TString(Prefix + "_jet_m"), "vector<float>", &jet_m);
    outTree->Branch(TString(Prefix + "_jet_constituents_jetIndex"), "vector<int>", &jet_constituents_jetIndex);
}

void Jet_Builder_data::fill_cell_var()
{
    jet_constituents_jetIndex = std::vector<int>( n_constituents, -1 );
    int size_jets = jets.size();
    for (int ijet = 0; ijet < size_jets; ijet++)
    {
        jet_pt.push_back(jets.at(ijet).pt());
        jet_eta.push_back(jets.at(ijet).eta());
        float phi = jets.at(ijet).phi();
        if (phi < -M_PI) phi += 2*M_PI;
        if (phi > M_PI) phi -= 2*M_PI;
        jet_phi.push_back(phi);
        jet_m.push_back(jets.at(ijet).m());
        int size_constituents = jets.at(ijet).constituents().size();
        std::vector<float> constituents;
        for (int icon = 0; icon < size_constituents; icon++)
        {
	    int constituent_index                        = jets.at(ijet).constituents().at(icon).user_index();
	    jet_constituents_jetIndex[constituent_index] = ijet;
        }
    }
}

float Jet_Builder_data::R_distance_func(float phi_1, float eta_1, float phi_2, float eta_2)//, float sigma_phi, float sigma_eta)
{
    float d_phi = fabs(phi_1 - phi_2);
    if (d_phi > M_PI)
        d_phi -= 2 * M_PI;
    float d_eta = eta_1 - eta_2;
    return sqrtf(sqr(d_phi) + sqr(d_eta));
}


void Jet_Builder_data::track_match_to_jet(std::vector<Track_struct> &track_list, float radius)
{

    int size_jets = jets.size();
    int size_track_list = track_list.size();

    for (int itrack = 0; itrack < size_track_list; itrack++) {
        Track_struct &track = track_list.at(itrack);
        if ((abs(track.pdgcode)!=13) && track.Is_Track_Useable()) {
            for (int ijet = 0; ijet < size_jets; ijet++) {
                float track_eta = -log(tan(track.theta*0.5));
                float R_distance = R_distance_func(track.phiHelix, track_eta,
                                                    jets.at(ijet).phi(),
                                                    jets.at(ijet).eta());
                if (track.Rprime_to_closest_jet.empty())
                {
                    track.Rprime_to_closest_jet.push_back(R_distance);
                    track.index_of_closest_jet.push_back(ijet);
                }
                else if (R_distance < track.Rprime_to_closest_jet.front())
                {
                    track.Rprime_to_closest_jet.front() = R_distance;
                    track.index_of_closest_jet.front() = ijet;
                }
            }
        }
        if (track.Rprime_to_closest_jet.empty())
        {
            track.Rprime_to_closest_jet.push_back(1000);
            track.index_of_closest_jet.push_back(-1);
        }  

        if (track.Rprime_to_closest_jet.at(0) < radius) {
            track.index_of_matched_jet.push_back(track.index_of_closest_jet.at(0));
        }
        else {
            track.index_of_matched_jet.push_back(-1);
        }
    }
}

/*
void Jet_Builder_data::flavour_labelling(std::vector<Track_struct> &track_list, )
{

    int size_jets = jets.size();
    int size_track_list = track_list.size();

    for (int ijet = 0; ijet < size_jets; ijet++) {
        for (int itrack = 0; itrack < size_track_list; itrack++) {
            Track_struct &track = track_list.at(itrack);
            if ( track.index_of_matched_jet == ijet ) { // found track
                int size_trajectories_list = fAllTrajectoryInfo.size();
                for (int iparticle = 0; iparticle < size_trajectories_list; iparticle++)
                {
                    if ( tracks.nFinal_State_Particles == iparticle ) { // found particle
                        std::vector<HepMC::GenParticle *> nodes =;
                    }
                    if fAllTrajectoryInfo.at(iparticle)
                    float px    = fAllTrajectoryInfo.at(iparticle).fMomentum.x();
                    float py    = fAllTrajectoryInfo.at(iparticle).fMomentum.y();
                }
            }
        }
    }
}
*/

void Jet_Builder_data::clear()
{
    jet_pt.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_m.clear();
    jet_constituents_jetIndex.clear();
}
