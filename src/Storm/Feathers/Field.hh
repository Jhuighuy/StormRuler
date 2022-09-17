/**
 *    ______             __     __  _____ _____
 *   / __/ /____ _____  / /__  /  |/  / // / _ \
 *  _\ \/  '_/ // / _ \/  '_/ / /|_/ / _  / // /
 * /___/_/\_\\_,_/_//_/_/\_\ /_/  /_/_//_/____/
 *
 * Copyright (c) 2019 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once
#ifndef FIELD_HH_
#define FIELD_HH_

#include "SkunkBase.hh"

namespace Storm {

#define FEATHERS_ALLOCA(T, size) \
    static_cast<T*>(alloca(sizeof(T)*(size)))

template<typename data_t>
class tGenericSubField {
private:
    friend class tGenericSubField<const data_t>;
    size_t m_num_vars;
    data_t* m_elements;

public:
    tGenericSubField(size_t num_vars, data_t* elements):
        m_num_vars(num_vars), m_elements(elements) {
    }

    tGenericSubField( // NOLINT(google-explicit-constructor)
           const tGenericSubField<std::remove_const_t<data_t>>& other):
        m_num_vars(other.m_num_vars), m_elements(other.m_elements) {
    }

    template<typename data_u = data_t,
        typename = std::enable_if_t<!std::is_const_v<data_u>>>
    tGenericSubField& operator=(const std::initializer_list<data_u>& other) {
        if (other.size() == 0) {
            std::fill_n(m_elements, m_num_vars, data_u{});
        } else {
            std::copy(other.begin(), other.end(), m_elements);
        }
        return *this;
    }

    template<typename data_u = data_t,
        typename = std::enable_if_t<!std::is_const_v<data_u>>>
    tGenericSubField& operator=(const tGenericSubField<data_u>& other) {
        std::copy_n(other.data(), m_num_vars, m_elements);
        return *this;
    }
    template<typename data_u = data_t,
        typename = std::enable_if_t<!std::is_const_v<data_u>>>
    tGenericSubField& operator=(const tGenericSubField<const data_u>& other) {
        std::copy_n(other.data(), m_num_vars, m_elements);
        return *this;
    }

    template<typename type_u = data_t>
    std::enable_if_t<!std::is_const_v<type_u>> fill(const type_u& value) {
        std::fill_n(m_elements, m_num_vars, value);
    }

    auto data() const {
        return m_elements;
    }

    auto& operator[](size_t i) const {
        return m_elements[i];
    }
};

using tScalarSubField = tGenericSubField<real_t>;
using tVectorSubField = tGenericSubField<vec3_t>;
using tMatrixSubField = tGenericSubField<mat3_t>;

using tScalarConstSubField = tGenericSubField<const real_t>;
using tVectorConstSubField = tGenericSubField<const vec3_t>;
using tMatrixConstSubField = tGenericSubField<const mat3_t>;

#define FEATHERS_TMP_SCALAR_FIELD(name, num_vars) \
    tScalarSubField name((num_vars), FEATHERS_ALLOCA(real_t, (num_vars)))

template<typename type_t, typename component_type_t = type_t,
    size_t num_components = sizeof(type_t) / sizeof(component_type_t)>
class tGenericField {
private:
    size_t m_num_vars;
    std::vector<component_type_t> m_elements;

public:
    tGenericField(size_t num_vars, size_t num_elements):
        m_num_vars(num_vars), m_elements(num_components*num_vars*num_elements) {
    }

    auto operator[](size_t element_index) {
        return tGenericSubField<type_t>(m_num_vars,
            reinterpret_cast<type_t*>(&m_elements[num_components*m_num_vars*element_index]));
    }
    auto operator[](size_t element_index) const {
        return tGenericSubField<const type_t>(m_num_vars,
            reinterpret_cast<const type_t*>(&m_elements[num_components*m_num_vars*element_index]));
    }

    void swap(tGenericField& other) {
        std::swap(m_num_vars, other.m_num_vars);
        std::swap(m_elements, other.m_elements);
    }
}; // class tGenericField

using tScalarField = tGenericField<real_t>;
using tVectorField = tGenericField<vec3_t, real_t>;
using tMatrixField = tGenericField<mat3_t, real_t>;

struct sFieldDesc {
    const char* name;
    size_t var_index;
    tScalarField* scalar;
}; // struct sFieldDesc

} // namespace feathers

#endif // FIELD_HH_
