// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file construct.hpp
 * \brief construct.hpp contains methods to construct indexes (compressed suffix arrays and trees).
 * \author Simon Gog
 */

#ifndef INCLUDED_SDSL_CONSTRUCT
#define INCLUDED_SDSL_CONSTRUCT

#include <iosfwd>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <type_traits>

#include "config.hpp"
#include "int_vector.hpp"
#include "int_vector_buffer.hpp"
#include "int_vector_mapper.hpp"
#include "io.hpp"
#include "memory_tracking.hpp"
#include "ram_fs.hpp"
#include "sdsl_concepts.hpp"
#include "util.hpp"

namespace sdsl
{

template <class int_vector>
bool contains_no_zero_symbol(int_vector const & text, std::string const & file)
{
    for (int_vector_size_type i = 0; i < text.size(); ++i)
    {
        if ((uint64_t)0 == text[i])
        {
            throw std::logic_error(std::string("Error: File \"") + file + "\" contains zero symbol.");
            return false;
        }
    }
    return true;
}

template <class int_vector>
void append_zero_symbol(int_vector & text)
{
    text.resize(text.size() + 1);
    text[text.size() - 1] = 0;
}

template <class t_index>
void construct(t_index & idx, std::string file, uint8_t num_bytes = 0, bool move_input = false)
{
    tMSS file_map;
    cache_config config;
    if (is_ram_file(file))
    {
        config.dir = "@";
        config.delete_data = move_input;
    }
    construct(idx, file, config, num_bytes);
}

template <class t_index, class t_data>
void construct_im(t_index & idx, t_data && data, uint8_t num_bytes = 0)
{
    std::string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
    store_to_file(data, tmp_file);
    construct(idx, tmp_file, num_bytes, std::is_rvalue_reference<t_data &&>::value);
    ram_fs::remove(tmp_file);
}

//! Constructs an index object of type t_index for a text stored on disk.
/*!
 * \param idx       	t_index object.  Any sdsl suffix array of suffix tree.
 * \param file      	Name of the text file. The representation of the file
 *                  	is dependent on the next parameter.
 * \
 * \param num_bytes 	If `num_bytes` equals 0, the file format is a serialized
 *				    	int_vector<>. Otherwise the file is interpreted as sequence
 *                  	of `num_bytes`-byte integer stored in big endian order.
 */
template <class t_index>
void construct(t_index & idx, std::string const & file, cache_config & config, uint8_t num_bytes = 0)
{
    // delegate to CSA or CST construction
    typename t_index::index_category index_tag;
    construct(idx, file, config, num_bytes, index_tag);
}

// Specialization for WTs
template <class t_index>
void construct(t_index & idx, std::string const & file, cache_config & config, uint8_t num_bytes, wt_tag)
{
    auto event = memory_monitor::event("construct wavelet tree");
    if ((t_index::alphabet_category::WIDTH == 8 and num_bytes <= 1)
        or (t_index::alphabet_category::WIDTH == 0 and num_bytes != 'd'))
    {
        int_vector_buffer<t_index::alphabet_category::WIDTH> text_buf(file,
                                                                      std::ios::in,
                                                                      1024 * 1024,
                                                                      num_bytes * 8,
                                                                      (bool)num_bytes);
        idx = t_index(text_buf.begin(), text_buf.end(), config.dir);
    }
    else
    {
        int_vector<t_index::alphabet_category::WIDTH> text;
        load_vector_from_file(text, file, num_bytes);
        std::string tmp_key = util::to_string(util::pid()) + "_" + util::to_string(util::id());
        std::string tmp_file_name = cache_file_name(tmp_key, config);
        store_to_file(text, tmp_file_name);
        util::clear(text);
        {
            int_vector_buffer<t_index::alphabet_category::WIDTH> text_buf(tmp_file_name);
            idx = t_index(text_buf.begin(), text_buf.end(), config.dir);
        }
        sdsl::remove(tmp_file_name);
    }
}


} // end namespace sdsl
#endif
