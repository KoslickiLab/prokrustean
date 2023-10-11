#ifndef APPLICATION_CDBG_COUNTING_HPP_
#define APPLICATION_CDBG_COUNTING_HPP_
#include <algorithm>
#include "../prokrustean.enhance.hpp"


int count_maximal_unitigs_single_k(int k, ProkrusteanEnhancement &prokrustean_ext, bool verbose=false){
    Prokrustean &prokrustean = prokrustean_ext.prokrustean;

    int cnt=0;
    vector<int> stats(7);
    vector<Region> spectrum;
    for(int i=0; i<prokrustean.sequence_count(); i++){
        Sequence seq = prokrustean.get_sequence(i);
        if(seq.size<k){
            continue;
        }
        // tip of sequence
        prokrustean.get_spectrum(seq, k-1, spectrum);
        if(spectrum[0].is_reflected){
            // cnt++;
            stats[0]++;
        }
    }
    for(int i=0; i<prokrustean.stratum_count(); i++){
        Stratum stra = prokrustean.get_stratum(i);
        if(stra.size<k-1){
            continue;
        }
        prokrustean.get_spectrum(stra, k-1, spectrum);
        if(stra.size>k-1){
            if(spectrum[0].is_reflected){
                if(prokrustean_ext.stratum_left_ext_count[i]==0){
                    // tip
                    // cout << "tip at stratum " << i << endl; 
                    // cnt++;
                    stats[1]++;
                } else if(prokrustean_ext.stratum_left_ext_count[i]>1){
                    // convergence
                    // cout << "convergence at stratum " << i << endl; 
                    // cnt++;
                    stats[2]++;
                }
            } 
            if(spectrum[spectrum.size()-1].is_reflected){
                    if(prokrustean_ext.stratum_right_ext_count[i]>1){
                    // divergence
                    // cout << "divergence at stratum " << i << " (" << (int)prokrustean_optional.stratum_right_ext_count[i] << ")" << endl;  
                    cnt+=prokrustean_ext.stratum_right_ext_count[i];
                    stats[3]+=prokrustean_ext.stratum_right_ext_count[i];
                }
            }
        } else {
            // special case
            if(prokrustean_ext.stratum_right_ext_count[i]>1){
                // divergence multiple -> convergence does not matter
                // cnt+=prokrustean_ext.stratum_right_ext_count[i];
                stats[4]+=prokrustean_ext.stratum_right_ext_count[i];
            } else if(prokrustean_ext.stratum_right_ext_count[i]==1){
                if(prokrustean_ext.stratum_left_ext_count[i]==0){
                    // divergence single -> convergence can work
                    // cout << "tip at stratum of k-1 " << i << endl;  
                    // cnt++;
                    stats[5]++;
                } else if(prokrustean_ext.stratum_left_ext_count[i]>1){
                    // cout << "convergence at stratum of k-1 " << i << endl;  
                    // cnt++;
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
    
#endif