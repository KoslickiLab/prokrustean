#ifndef APPLICATION_CDBG_COUNTING_HPP_
#define APPLICATION_CDBG_COUNTING_HPP_
#include <algorithm>
#include "../prokrustean.support.hpp"


int count_maximal_unitigs_single_k(int k, ProkrusteanExtension &prokrustean_ext, bool verbose=false){
    Prokrustean &prokrustean = prokrustean_ext.prokrustean;
    optional<int> turn_on=nullopt;
    if(turn_on.has_value()){
        cout << "filter turned on: investigating maximal unitig starting point by components" << endl;
    }

    int cnt=0;
    vector<int> stats(7);
    vector<Edge> spectrum;
    for(int i=0; i<prokrustean.sequence_count; i++){
        SequenceVertex seq = prokrustean.get_sequence(i);
        if(seq.size<k-1){
            continue;
        }
        // tip of sequence
        prokrustean.get_spectrum(seq, k-1, spectrum);
        if(spectrum[0].is_reflected){
            // maximal case 1
            if(!turn_on.has_value() || turn_on.value()==1){
                cnt++;
            }
            stats[0]++;
        }
    }
    for(int i=0; i<prokrustean.stratum_count; i++){
        StratumVertex stra = prokrustean.get_stratum(i);
        if(stra.size<k-1){
            continue;
        }
        prokrustean.get_spectrum(stra, k-1, spectrum);
        if(stra.size>k-1){
            if(spectrum[0].is_reflected){
                if(prokrustean_ext.stratum_left_ext_count[i]==0){
                    // tip
                    // maximal case 2
                    if(!turn_on.has_value() || turn_on.value()==2){
                        cnt++;
                    }
                    stats[1]++;
                } else if(prokrustean_ext.stratum_left_ext_count[i]>1){
                    // convergence
                    // maximal case 3
                    if(!turn_on.has_value() || turn_on.value()==3){
                        cnt++;
                    }
                    stats[2]++;
                }
            } 
            if(spectrum[spectrum.size()-1].is_reflected){
                    if(prokrustean_ext.stratum_right_ext_count[i]>1){
                    // divergence
                    // maximal case 4
                    if(!turn_on.has_value() || turn_on.value()==4){
                        cnt+=prokrustean_ext.stratum_right_ext_count[i];
                    }
                    stats[3]+=prokrustean_ext.stratum_right_ext_count[i];
                }
            }
        } else {
            // special case k-1
            if(prokrustean_ext.stratum_right_ext_count[i]>1){
                // divergence multiple -> convergence does not matter
                // maximal case 5
                if(!turn_on.has_value() || turn_on.value()==5){
                    cnt+=prokrustean_ext.stratum_right_ext_count[i];
                }
                stats[4]+=prokrustean_ext.stratum_right_ext_count[i];
            } else if(prokrustean_ext.stratum_right_ext_count[i]==1){
                
                if(prokrustean_ext.stratum_left_ext_count[i]==0){
                    // cout << "tip at stratum of k-1 " << i << endl; divergence single -> convergence can work
                    /*stratum tip at stratum of k-1 && right ext count 1 is impossible:
                    If right ext count is 1, that means at least one extension was terminal (to make stratum)
                    and being left side tip means all occurrences are terminal -> meaning the string is k-1
                    */
                    assert(false); // leave this as assert
                    // cnt++;
                    stats[5]++;
                } else if(prokrustean_ext.stratum_left_ext_count[i]>1){
                    // cout << "convergence at stratum of k-1 " << i << endl;  
                    // maximal case 6
                    if(!turn_on.has_value() || turn_on.value()==6){
                        cnt++;
                    }
                    stats[6]++;
                }
            }
        }
    }
    if(verbose){
        cout << "maximal unitigs of " << k << ": " << cnt << endl;
        for(int i=0; i<stats.size(); i++){
            if(i==0) cout << "tips at sequence " << stats[i] << endl; 
            if(i==1) cout << "tips at stratum " << stats[i] << endl; 
            if(i==2) cout << "convergences at stratum " << stats[i] << endl; 
            if(i==3) cout << "divergneces at stratum " << stats[i] << endl; 
            if(i==4) cout << "divergneces at stratum of k-1 " << stats[i] << endl; 
            if(i==5) cout << "tips at stratum of k-1 " << stats[i] << endl; 
            if(i==6) cout << "convergences at stratum of k-1 " << stats[i] << endl; 
        }
    }
    return cnt;
}

int _count_max_unitig_start_at_seq(int k, SeqId seq_id, ProkrusteanExtension &ext){
    if(ext.prokrustean.sequences__size[seq_id]>=k-1 
    && ext.seq__refracted_at_front_of_cover(seq_id, k-1)){
        return 1;
    } else {
        return 0;
    }
}

int _count_max_unitig_start_at_stratum(int k, StratumId stratum_id, ProkrusteanExtension &ext){
    int stratum_size=ext.prokrustean.stratums__size[stratum_id];
    int cnt=0;
    if(stratum_size<k-1){
        return 0;
    } else if(stratum_size>k-1){
        if(ext.stratum__refracted_at_front_of_cover(stratum_id, k-1)){
            // convergence or tip
            cnt += ext.stratum_left_ext_count[stratum_id]!=1? 1 : 0;
        }
        if(ext.stratum__refracted_at_back_of_cover(stratum_id, k-1)){
            // divergence
            cnt += ext.stratum_right_ext_count[stratum_id]>1? ext.stratum_right_ext_count[stratum_id] : 0;
        }
    } else { // k-1 case
        if(ext.stratum_right_ext_count[stratum_id]>1){
            // divergence implication
            cnt += ext.stratum_right_ext_count[stratum_id];    
        } else if(ext.stratum_right_ext_count[stratum_id]==1 &&  ext.stratum_left_ext_count[stratum_id]){
            // convergence implication
            cnt += 1;
        }
    }
    return cnt;
}
int _count_events_in_range(int k, ProkrusteanExtension &prokrustean_ext, bool verbose=false){
    Prokrustean &prokrustean = prokrustean_ext.prokrustean;
    optional<int> turn_on=nullopt;
    if(turn_on.has_value()){
        cout << "filter turned on: investigating maximal unitig starting point by components" << endl;
    }

    int cnt=0;
    vector<int> stats(7);
    vector<Edge> spectrum;
    for(int i=0; i<prokrustean.sequence_count; i++){
        SequenceVertex seq = prokrustean.get_sequence(i);
        if(seq.size<k-1){
            continue;
        }
        // tip of sequence
        prokrustean.get_spectrum(seq, k-1, spectrum);
        if(spectrum[0].is_reflected){
            // maximal case 1
            if(!turn_on.has_value() || turn_on.value()==1){
                cnt++;
            }
            stats[0]++;
        }
    }
    for(int i=0; i<prokrustean.stratum_count; i++){
        StratumVertex stra = prokrustean.get_stratum(i);
        if(stra.size<k-1){
            continue;
        }
        prokrustean.get_spectrum(stra, k-1, spectrum);
        if(stra.size>k-1){
            if(spectrum[0].is_reflected){
                if(prokrustean_ext.stratum_left_ext_count[i]==0){
                    // tip
                    // maximal case 2
                    if(!turn_on.has_value() || turn_on.value()==2){
                        cnt++;
                    }
                    stats[1]++;
                } else if(prokrustean_ext.stratum_left_ext_count[i]>1){
                    // convergence
                    // maximal case 3
                    if(!turn_on.has_value() || turn_on.value()==3){
                        cnt++;
                    }
                    stats[2]++;
                }
            } 
            if(spectrum[spectrum.size()-1].is_reflected){
                    if(prokrustean_ext.stratum_right_ext_count[i]>1){
                    // divergence
                    // maximal case 4
                    if(!turn_on.has_value() || turn_on.value()==4){
                        cnt+=prokrustean_ext.stratum_right_ext_count[i];
                    }
                    stats[3]+=prokrustean_ext.stratum_right_ext_count[i];
                }
            }
        } else {
            // special case k-1
            if(prokrustean_ext.stratum_right_ext_count[i]>1){
                // divergence multiple -> convergence does not matter
                // maximal case 5
                if(!turn_on.has_value() || turn_on.value()==5){
                    cnt+=prokrustean_ext.stratum_right_ext_count[i];
                }
                stats[4]+=prokrustean_ext.stratum_right_ext_count[i];
            } else if(prokrustean_ext.stratum_right_ext_count[i]==1){
                
                if(prokrustean_ext.stratum_left_ext_count[i]==0){
                    // cout << "tip at stratum of k-1 " << i << endl; divergence single -> convergence can work
                    /*stratum tip at stratum of k-1 && right ext count 1 is impossible:
                    If right ext count is 1, that means at least one extension was terminal (to make stratum)
                    and being left side tip means all occurrences are terminal -> meaning the string is k-1
                    */
                    assert(false); // leave this as assert
                    // cnt++;
                    stats[5]++;
                } else if(prokrustean_ext.stratum_left_ext_count[i]>1){
                    // cout << "convergence at stratum of k-1 " << i << endl;  
                    // maximal case 6
                    if(!turn_on.has_value() || turn_on.value()==6){
                        cnt++;
                    }
                    stats[6]++;
                }
            }
        }
    }
    if(verbose){
        cout << "maximal unitigs of " << k << ": " << cnt << endl;
        for(int i=0; i<stats.size(); i++){
            if(i==0) cout << "tips at sequence " << stats[i] << endl; 
            if(i==1) cout << "tips at stratum " << stats[i] << endl; 
            if(i==2) cout << "convergences at stratum " << stats[i] << endl; 
            if(i==3) cout << "divergneces at stratum " << stats[i] << endl; 
            if(i==4) cout << "divergneces at stratum of k-1 " << stats[i] << endl; 
            if(i==5) cout << "tips at stratum of k-1 " << stats[i] << endl; 
            if(i==6) cout << "convergences at stratum of k-1 " << stats[i] << endl; 
        }
    }
    return cnt;
}

void count_maximal_unitigs_range_of_k(uint64_t from, uint64_t to, ProkrusteanExtension &ext, vector<uint64_t> &output){
    // Definitions:
    assert(from>0 && from<=to);
    output.clear();
    output.resize(to+1, 0);

    vector<int64_t> partial_C(to+1, 0);
    vector<int64_t> partial_partial_C(to+1, 0);
    vector<int64_t> dOutput(to+1, 0);
    
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        output[from]+=_count_max_unitig_start_at_seq(from, i, ext);
        // _viable_kmer_decreases(sequence.size, from, sequence.s_edges, partial_partial_C);
        // _each_stra_rgn_range_decided_by_intersection(sequence.size, from, sequence.s_edges, partial_partial_C);
    }

    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        output[from]+=_count_max_unitig_start_at_stratum(from, i, ext);
    }
    for(int k=from+1; k<to+1; k++){
        dOutput[k]=dOutput[k-1] + partial_partial_C[k];
        output[k]=output[k-1] + dOutput[k];
    }
}

    
#endif