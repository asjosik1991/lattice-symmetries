// Copyright (c) 2019-2020, Tom Westerhout
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include "lattice_symmetries/lattice_symmetries.h"
#include <immintrin.h>
#include <cstdint>

namespace lattice_symmetries {

inline constexpr int batch_size = 8;
static_assert(4 * sizeof(__m128i) == batch_size * sizeof(uint64_t));
static_assert(2 * sizeof(__m256i) == batch_size * sizeof(uint64_t));

auto benes_forward(uint64_t& x, uint64_t const masks[], unsigned size,
                   uint16_t const deltas[]) noexcept -> void;

auto benes_forward(uint64_t x[batch_size], uint64_t const (*masks)[batch_size], unsigned size,
                   uint16_t const deltas[]) noexcept -> void;

auto benes_forward(ls_bits512& x, ls_bits512 const masks[], unsigned size,
                   uint16_t const deltas[]) noexcept -> void;

auto search_sorted(uint64_t const* first, uint64_t const* last, uint64_t key) noexcept
    -> uint64_t const*;

namespace detail {
    auto benes_forward_sse2(uint64_t x[batch_size], uint64_t const (*masks)[batch_size],
                            unsigned size, uint16_t const deltas[]) noexcept -> void;
    auto benes_forward_avx(uint64_t x[batch_size], uint64_t const (*masks)[batch_size],
                           unsigned size, uint16_t const deltas[]) noexcept -> void;
    auto benes_forward_avx2(uint64_t x[batch_size], uint64_t const (*masks)[batch_size],
                            unsigned size, uint16_t const deltas[]) noexcept -> void;

    auto benes_forward_sse2(void* x, uint64_t const (*masks)[batch_size], unsigned size,
                            uint16_t const deltas[]) noexcept -> void;
    auto benes_forward_avx(void* x, uint64_t const (*masks)[batch_size], unsigned size,
                           uint16_t const deltas[]) noexcept -> void;
    auto benes_forward_avx2(void* x, uint64_t const (*masks)[batch_size], unsigned size,
                            uint16_t const deltas[]) noexcept -> void;

    auto benes_forward_512_sse2(ls_bits512& x, ls_bits512 const masks[], unsigned size,
                                uint16_t const deltas[]) noexcept -> void;
    auto benes_forward_512_avx(ls_bits512& x, ls_bits512 const masks[], unsigned size,
                               uint16_t const deltas[]) noexcept -> void;
    auto benes_forward_512_avx2(ls_bits512& x, ls_bits512 const masks[], unsigned size,
                                uint16_t const deltas[]) noexcept -> void;

    auto search_sorted_avx2(uint64_t const* first, uint64_t const* last, uint64_t key) noexcept
        -> uint64_t const*;
    auto search_sorted_avx(uint64_t const* first, uint64_t const* last, uint64_t key) noexcept
        -> uint64_t const*;
    auto search_sorted_sse2(uint64_t const* first, uint64_t const* last, uint64_t key) noexcept
        -> uint64_t const*;

} // namespace detail

} // namespace lattice_symmetries
