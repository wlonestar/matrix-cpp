//
// Created by wjl on 2022/3/4.
//

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#pragma once

#include <cstdarg>
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

namespace matrix {

    using Type = enum { Identity, Zero, One, NegaOne };

    template<typename T>
    class matrix {
    private:
        int _row{};
        int _col{};
        vector<vector<T>> _data;

        void init(int row, int col);
        void copyfrom(const matrix<T> &m);
        bool elemeq(const vector<vector<T>> &data);
        bool rcvaild(int row, int col);
        matrix<T> &opposite();
    public:
        explicit matrix(int row = 1, int col = 1, Type type = Zero);
        matrix(int row, int col, T e...);
        matrix(const matrix &m);
        ~matrix();
        void resize(int row, int col, Type type = Zero);
        void set(int row, int col, const T &value);
        T get(int row, int col);
        [[nodiscard]] int rank() const;
        T det() const;
        vector<T> row(int row);
        vector<T> col(int col);
        bool operator==(const matrix<T> &m);
        matrix<T> &operator=(const matrix<T> &m);
        matrix<T> &operator+();
        matrix<T> &operator+=(const matrix<T> &m);
        matrix<T> &operator-();
        matrix<T> &operator-=(const matrix<T> &m);
        matrix<T> &operator*(T num);
        matrix<T> &operator*=(T num);
        matrix<T> &operator/(T num);
        matrix<T> &operator/=(T num);
        void mult(matrix<T> &m);
        matrix<T> transform();
        T cofactor(int row, int col) const;
        T algebraic_cofactor(int row, int col) const;
        matrix<T> cosubmatrix();
        matrix<T> adjoint();
        matrix<T> inverse();
        [[nodiscard]] bool is_invertible() const;
        [[nodiscard]] bool is_symmetric() const;
        void partitioned4(matrix<T> &m11, matrix<T> &m12, matrix<T> &m21, matrix<T> &m22);
        void print_row(int row);
        void print_col(int col);
        void print();
    };

    template<typename T>
    void matrix<T>::init(int row, int col) {
        _row = row;
        _col = col;
        _data.resize(_row);
        for (int i = 0; i < _data.size(); i++) {
            _data[i].resize(_col);
        }
    }

    template<typename T>
    void matrix<T>::copyfrom(const matrix<T> &m) {
        init(m._row, m._col);
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                _data[i][j] = m._data[i][j];
            }
        }
    }

    template<typename T>
    bool matrix<T>::elemeq(const vector<vector<T>> &data) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                if (_data[i][j] != data[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    template<typename T>
    bool matrix<T>::rcvaild(int row, int col) {
        if (row > _row || col > _col || row < 0 || col < 0) {
            return false;
        }
        return true;
    }

    template<typename T>
    matrix<T> &matrix<T>::opposite() {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                _data[i][j] = -_data[i][j];
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T>::matrix(int row, int col, Type type) {
        init(row, col);
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                if (type == Identity) {
                    i == j ? _data[i][j] = 1 : _data[i][j] = 0;
                } else if (type == Zero) {
                    _data[i][j] = 0;
                } else if (type == One){
                    _data[i][j] = 1;
                } else if (type == NegaOne) {
                    _data[i][j] = -1;
                }
            }
        }
    }

    template<typename T>
    matrix<T>::matrix(int row, int col, T e...) {
        init(row, col);
        va_list args;
        va_start(args, e);
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                _data[i][j] = va_arg(args, T);
            }
        }
        va_end(args);
    }

    template<typename T>
    matrix<T>::matrix(const matrix &m) {
        copyfrom(m);
    }

    template<typename T>
    matrix<T>::~matrix() {
        _data.resize(0);
    }

    template<typename T>
    void matrix<T>::resize(int row, int col, Type type) {
        _data.resize(0);
        init(row, col);
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                if (type == Identity) {
                    i == j ? _data[i][j] = 1 : _data[i][j] = 0;
                } else if (type == Zero) {
                    _data[i][j] = 0;
                } else if (type == One){
                    _data[i][j] = 1;
                } else if (type == NegaOne) {
                    _data[i][j] = -1;
                }
            }
        }
    }

    template<typename T>
    void matrix<T>::set(int row, int col, const T &value) {
        if (rcvaild(row, col)) {
            _data[row][col] = value;
        }
    }

    template<typename T>
    T matrix<T>::get(int row, int col) {
        if (rcvaild(row, col)) {
            return _data[row][col];
        }
        return -1;
    }

    template<typename T>
    int matrix<T>::rank() const {
        return 0;
    }

    template<typename T>
    T matrix<T>::det() const {
        if (_row == _col) {
            if (_row == 1) {
                return _data[0][0];
            } else if (_row == 2) {
                T a = _data[0][0]; T b = _data[0][1];
                T c = _data[1][0]; T d = _data[1][1];
                return (a * d - b * c);
            } else if (_row == 3) {
                T a = _data[0][0]; T b = _data[0][1]; T c = _data[0][2];
                T d = _data[1][0]; T e = _data[1][1]; T f = _data[1][2];
                T g = _data[2][0]; T h = _data[2][1]; T k = _data[2][2];
                return (a * e * k + c * d * h + b * f * g - c * e *g - b * d * k - a * f * h);
            } else {
                T sum = 0;
                for (int i = 0; i < _row; i++) {
                    sum += algebraic_cofactor(i, 0);
                }
                return sum;
            }
        }
        return -1;
    }

    template<typename T>
    vector<T> matrix<T>::row(int row) {
        vector<T> v;
        if (row > _row || row < 0) {
            v.push_back(-1);
            return v;
        }
        return _data[row];
    }

    template<typename T>
    vector<T> matrix<T>::col(int col) {
        if (col > _col || col < 0) {
            vector<T> v;
            v.push_back(-1);
            return v;
        }
        vector<T> v;
        for (int i = 0; i < _row; i++) {
            v.push_back(_data[i][col]);
        }
        return v;
    }

    template<typename T>
    bool matrix<T>::operator==(const matrix<T> &m) {
        return _row == m._row && _col == m._col && elemeq(m._data);
    }

    template<typename T>
    matrix<T> &matrix<T>::operator=(const matrix<T> &m) {
        copyfrom(m);
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator+() {
        return opposite();
    }

    template<typename T>
    matrix<T> &matrix<T>::operator+=(const matrix<T> &m) {
        if (_row == m._row && _col == m._col) {
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _col; j++) {
                    _data[i][j] += m._data[i][j];
                }
            }
            return *this;
        }
    }

    template<typename T>
    matrix<T> &matrix<T>::operator-() {
        return opposite();
    }

    template<typename T>
    matrix<T> &matrix<T>::operator-=(const matrix<T> &m) {
        if (_row == m._row && _col == m._col) {
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _col; j++) {
                    _data[i][j] -= m._data[i][j];
                }
            }
            return *this;
        }
    }

    template<typename T>
    matrix<T> &matrix<T>::operator*(T num) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                _data[i][j] *= num;
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator*=(T num) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                _data[i][j] *= num;
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator/(T num) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                _data[i][j] /= num;
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator/=(T num) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                _data[i][j] /= num;
            }
        }
        return *this;
    }

    template<typename T>
    void matrix<T>::mult(matrix<T> &m) {
        if (_col == m._row) {
            matrix<T> matrix(_row, m._col, Zero);
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < m._col; j++) {
                    for (int k = 0; k < _col; k++) {
                        matrix._data[i][j] += _data[i][k] * m._data[k][j];
                    }
                }
            }
            copyfrom(matrix);
        }
    }

    template<typename T>
    matrix<T> matrix<T>::transform() {
        matrix<T> matrix(_col, _row);
        matrix.init(_col, _row);
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                matrix._data[j][i] = _data[i][j];
            }
        }
        return matrix;
    }

    template<typename T>
    T matrix<T>::cofactor(int row, int col) const {
        matrix<T> matrix(_row - 1, _col - 1);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                matrix._data[i][j] = _data[i][j];
            }
        }
        for (int i = 0; i < row; i++) {
            for (int j = col + 1; j < _col; j++) {
                matrix._data[i][j - 1] = _data[i][j];
            }
        }
        for (int i = row + 1; i < _row; i++) {
            for (int j = 0; j < col; j++) {
                matrix._data[i - 1][j] = _data[i][j];
            }
        }
        for (int i = row + 1; i < _row; i++) {
            for (int j = col + 1; j < _col; j++) {
                matrix._data[i - 1][j - 1] = _data[i][j];
            }
        }
        return matrix.det();
    }

    template<typename T>
    T matrix<T>::algebraic_cofactor(int row, int col) const {
        return pow(-1, row + col + 2) * _data[row][col] * cofactor(row, col);
    }

    template<typename T>
    matrix<T> matrix<T>::cosubmatrix() {
        matrix<T> matrix(_row, _col);
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                matrix._data[i][j] = algebraic_cofactor(i, j);
            }
        }
        return matrix;
    }

    template<typename T>
    matrix<T> matrix<T>::adjoint() {
        return cosubmatrix().transform();
    }

    template<typename T>
    matrix<T> matrix<T>::inverse() {
        if (det() != 0) {
            matrix<T> matrix(_row, _col);
            ::matrix::matrix<T> adj = adjoint();
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _col; j++) {
                    matrix._data[i][j] = adj._data[i][j] / det();
                }
            }
            return matrix;
        }
    }

    template<typename T>
    bool matrix<T>::is_invertible() const {
        return det() != 0;
    }

    template<typename T>
    bool matrix<T>::is_symmetric() const {
        matrix<T> m(*this);
        matrix<T> n(*this);
        n = n.transform();
        return m == n;
    }

    template<typename T>
    void matrix<T>::partitioned4(matrix<T> &m11, matrix<T> &m12, matrix<T> &m21,
                                 matrix<T> &m22) {
        int r1 = _row / 2;
        int r2 = _row - _row / 2;
        int c1 = _col / 2;
        int c2 = _col - _col / 2;
        m11.init(r1, c1);
        m12.init(r1, c2);
        m21.init(r2, c1);
        m22.init(r2, c2);
        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c1; j++) {
                m11._data[i][j] = _data[i][j];
            }
        }
        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c2; j++) {
                m12._data[i][j] = _data[i][j + c1];
            }
        }
        for (int i = 0; i < r2; i++) {
            for (int j = 0; j < c1; j++) {
                m21._data[i][j] = _data[i + r1][j];
            }
        }
        for (int i = 0; i < r2; i++) {
            for (int j = 0; j < c2; j++) {
                m22._data[i][j] = _data[i + r1][j + c1];
            }
        }
    }

    template<typename T>
    void matrix<T>::print_row(int row) {
        vector<T> v = _data[row];
        cout << "----------\n";
        for (int i = 0; i < v.size(); i++) {
            cout << v[i] << " ";
        }
        cout << "\n----------\n";
    }

    template<typename T>
    void matrix<T>::print_col(int col) {
        vector<T> v(col);
        for (int i = 0; i < _row; i++) {
            cout << "|" << _data[i][col] << "|\n";
        }
    }

    template<typename T>
    void matrix<T>::print() {
        cout << "/----------\\\n";
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _col; j++) {
                cout << _data[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\\----------/\n";
    }

}

#endif //MATRIX_MATRIX_H
