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
                if(prokrustean_ext.prokrustean.get_left_cnt(i)==0){
                    // tip
                    // maximal case 2
                    if(!turn_on.has_value() || turn_on.value()==2){
                        cnt++;
                    }
                    stats[1]++;
                } else if(prokrustean_ext.prokrustean.get_left_cnt(i)>1){
                    // convergence
                    // maximal case 3
                    if(!turn_on.has_value() || turn_on.value()==3){
                        cnt++;
                    }
                    stats[2]++;
                }
            } 
            if(spectrum[spectrum.size()-1].is_reflected){
                    if(prokrustean_ext.prokrustean.get_right_cnt(i)>1){
                    // divergence
                    // maximal case 4
                    if(!turn_on.has_value() || turn_on.value()==4){
                        cnt+=prokrustean_ext.prokrustean.get_right_cnt(i);
                    }
                    stats[3]+=prokrustean_ext.prokrustean.get_right_cnt(i);
                }
            }
        } else {
            // special case k-1
            if(prokrustean_ext.prokrustean.get_right_cnt(i)>1){
                // divergence multiple -> convergence does not matter
                // maximal case 5
                if(!turn_on.has_value() || turn_on.value()==5){
                    cnt+=prokrustean_ext.prokrustean.get_right_cnt(i);
                }
                stats[4]+=prokrustean_ext.prokrustean.get_right_cnt(i);
            } else if(prokrustean_ext.prokrustean.get_right_cnt(i)==1){
                
                if(prokrustean_ext.prokrustean.get_left_cnt(i)==0){
                    // cout << "tip at stratum of k-1 " << i << endl; divergence single -> convergence can work
                    /*stratum tip at stratum of k-1 && right ext count 1 is impossible:
                    If right ext count is 1, that means at least one extension was terminal (to make stratum)
                    and being left side tip means all occurrences are terminal -> meaning the string is k-1
                    */
                    assert(false); // leave this as assert
                    // cnt++;
                    stats[5]++;
                } else if(prokrustean_ext.prokrustean.get_left_cnt(i)>1){
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
            cnt += ext.prokrustean.get_left_cnt(stratum_id)!=1? 1 : 0;
        }
        if(ext.stratum__refracted_at_back_of_cover(stratum_id, k-1)){
            // divergence
            cnt += ext.prokrustean.get_right_cnt(stratum_id)>1? ext.prokrustean.get_right_cnt(stratum_id) : 0;
        }
    } else { // k-1 case
        if(ext.prokrustean.get_right_cnt(stratum_id)>1){
            // divergence implication
            cnt += ext.prokrustean.get_right_cnt(stratum_id);    
        } else if(ext.prokrustean.get_right_cnt(stratum_id)==1 &&  ext.prokrustean.get_left_cnt(stratum_id)>1){
            // convergence implication
            cnt += 1;
        }
    }
    return cnt;
}

void _count_tip_events_in_range_at_seq(int k, SeqId seq_id, ProkrusteanExtension &ext, vector<uint64_t> &partial_C){
    auto seq_size=ext.prokrustean.sequences__size[seq_id];
    if(seq_size<k){
        return;
    }
    // tip -> available until reaching sequence length,
    // available from a valid refracted region starts - no stratified at pos 0 with L
    auto first_rgn=ext.seq__first_stratified(seq_id);
    if(first_rgn.has_value() && first_rgn.value().from==0 && first_rgn.value().size()>=k-1){
        bool one_stratified_covers_all_string=first_rgn.value().size()==seq_size; 
        if(one_stratified_covers_all_string){
            // there exists a sequence exactly matching with this.
            // then the common repeat(=stratum) governs the maximal unitig pattern
            return;
        }
        // from L refracted exists
        auto L=first_rgn.value().size()+2;
        if(L<partial_C.size()){
            partial_C[L]++;
        }
    }

    if(seq_size+2<partial_C.size()){
        partial_C[seq_size+2]--;
    }
}

void _count_divergence_events_in_range_at_stratum(int k, StratumId stratum_id, ProkrusteanExtension &ext, vector<uint64_t> &partial_C){
    auto stratum_size=ext.prokrustean.stratums__size[stratum_id];
    if(stratum_size<k-1){
        return;
    }
    if(ext.prokrustean.get_right_cnt(stratum_id)<=1){
        return;
    }
    // divergence
    // available from a valid refracted region ends - no stratified at last pos with L
    auto last_rgn=ext.stratum__last_stratified(stratum_id);
    if(last_rgn.has_value() && last_rgn.value().to==stratum_size && last_rgn.value().size()>=k-1){
        // from L refracted exists
        auto L=last_rgn.value().size()+2; // because last_region.value() is k-1 when L=last_region.value()+1
        if(L<partial_C.size()){
            partial_C[L]+=ext.prokrustean.get_right_cnt(stratum_id);
        }
    }
    // until it becomes less than k-2
    if(stratum_size+2<partial_C.size()){
        partial_C[stratum_size+2]-=ext.prokrustean.get_right_cnt(stratum_id);
    }
}

void _count_tip_events_in_range_at_stratum(int k, StratumId stratum_id, ProkrusteanExtension &ext, vector<uint64_t> &partial_C){
    auto stratum_size=ext.prokrustean.stratums__size[stratum_id];
    if(stratum_size<k-1){
        return;
    }
    if(ext.prokrustean.get_left_cnt(stratum_id)!=0){
        return;
    }
    // tip -> available until reaching sequence length,
    // available from a valid refracted region starts - no stratified at pos 0 with L
    auto first_rgn=ext.stratum__first_stratified(stratum_id);
    if(first_rgn.has_value() && first_rgn.value().from==0 && first_rgn.value().size()>=k-1){
        // from L refracted exists
        auto L=first_rgn.value().size()+2; // because last_region.value() is k-1 when L=last_region.value()+1
        if(L<partial_C.size()){
            partial_C[L]++;
        }
    } 
    // until it becomes less than k-1 (not k-2 unlike divergence)
    if(stratum_size+1<partial_C.size()){
        partial_C[stratum_size+1]--;
    }
}

void _count_convergence_events_in_range_at_stratum(int k, StratumId stratum_id, ProkrusteanExtension &ext, vector<uint64_t> &partial_C){
    auto stratum_size=ext.prokrustean.stratums__size[stratum_id];
    if(stratum_size<k-1){
        return;
    }
    if(ext.prokrustean.get_left_cnt(stratum_id)<=1){
        return;
    }
    // convergence
    // available from a valid refracted region starts - no stratified at first pos with L
    auto first_rgn=ext.stratum__first_stratified(stratum_id);
    if(first_rgn.has_value() && first_rgn.value().from==0 && first_rgn.value().size()>=k-1){
        // from L refracted exists
        auto L=first_rgn.value().size()+2; // because last_region.value() is k-1 when L=last_region.value()+1
        if(L<partial_C.size()){
            partial_C[L]++;
        }
    }
    
    // special case considering k-1
    if(ext.prokrustean.get_right_cnt(stratum_id)==1){
        // until it becomes less than k-1
        if(stratum_size+2<partial_C.size()){
            partial_C[stratum_size+2]--;
        }
    } else {
        // until it becomes less than k-2
        if(stratum_size+1<partial_C.size()){
            partial_C[stratum_size+1]--;
        }
    }
}


void count_maximal_unitigs_range_of_k(uint64_t from, uint64_t to, ProkrusteanExtension &ext, vector<uint64_t> &output){
    // Definitions:
    assert(from>1 && from<=to);
    output.clear();
    output.resize(to+1, 0);

    vector<uint64_t> partial_C(to+1, 0);
    
    for(int i=0; i<ext.prokrustean.sequence_count; i++){
        output[from]+=_count_max_unitig_start_at_seq(from, i, ext);
        _count_tip_events_in_range_at_seq(from, i, ext, partial_C);
    }

    for(int i=0; i<ext.prokrustean.stratum_count; i++){
        output[from]+=_count_max_unitig_start_at_stratum(from, i, ext);
        _count_tip_events_in_range_at_stratum(from, i, ext, partial_C);
        _count_divergence_events_in_range_at_stratum(from, i, ext, partial_C);
        _count_convergence_events_in_range_at_stratum(from, i, ext, partial_C);
    }
    
    for(int k=from+1; k<to+1; k++){
        output[k]=output[k-1] + partial_C[k];
    }
}

void count_maximal_unitigs_range_of_k_parallel(uint64_t from, uint64_t to, ProkrusteanExtension &ext, vector<uint64_t> &output, int thread_cnt){
    // Definitions:
    assert(from>1 && from<=to);
    output.clear();
    output.resize(to+1, 0);

    vector<vector<uint64_t>> outputs(thread_cnt);
    vector<vector<uint64_t>> partial_Cs(thread_cnt);
    for(int i=0; i<thread_cnt; i++){
        outputs[i].resize(to+1,0);
        partial_Cs[i].resize(to+1,0);
    }

    vector<future<void>> futures;
    auto func_ = [](ProkrusteanExtension &ext, uint64_t from, uint8_t thread_idx, uint8_t thread_cnt, vector<uint64_t> &output, vector<uint64_t> &partial_C) {
        Vertex vertex;
        for(int i=thread_idx; i<ext.prokrustean.sequence_count; i+=thread_cnt){
            output[from]+=_count_max_unitig_start_at_seq(from, i, ext);
            _count_tip_events_in_range_at_seq(from, i, ext, partial_C);
        }

        for(int i=thread_idx; i<ext.prokrustean.stratum_count; i+=thread_cnt){
            output[from]+=_count_max_unitig_start_at_stratum(from, i, ext);
            _count_tip_events_in_range_at_stratum(from, i, ext, partial_C);
            _count_divergence_events_in_range_at_stratum(from, i, ext, partial_C);
            _count_convergence_events_in_range_at_stratum(from, i, ext, partial_C);
        }
    };
    for(int i=0; i<thread_cnt; i++){futures.push_back(
        std::async(std::launch::async, func_, ref(ext), from, i, thread_cnt, ref(outputs[i]), ref(partial_Cs[i]))
    );}
    for (auto &f : futures) {f.wait();}
    
    vector<uint64_t> partial_C(to+1, 0);
    for(int i=0; i<thread_cnt; i++){
        for(int k=0; k<to+1; k++){
            output[k]+=outputs[i][k];
            partial_C[k]+=partial_Cs[i][k];
        }
    }
    for(int k=from+1; k<to+1; k++){
        output[k]=output[k-1] + partial_C[k];
    }

}
    
#endif