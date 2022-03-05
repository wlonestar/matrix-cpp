//
// Created by wjl on 2022/3/4.
//

#include "matrix.h"

template<typename T>
void matrix_constructor_test() {
    cout << "\n==== test matrix constructor:\n";
    cout << "row: 1 & col: 2\n";
    matrix::matrix<T> matrix(1, 2, matrix::One);
    matrix.print();
    matrix.resize(3, 1, matrix::NegaOne);
    cout << "row: 3 & col: 1\n";
    matrix.print();
    matrix.resize(3, 3, matrix::Identity);
    cout << "row: 3 & col: 3\n";
    matrix.print();
    cout << "==== done\n";
}

template<typename T>
void matrix_getter_setter_test() {
    cout << "\n==== test matrix getter and setter methods: \n";
    cout << "matrix row:2 & col: 2\n";
    matrix::matrix<T> matrix(2, 2, matrix::One);
    matrix.print();
    cout << "get element at (row = 0, col = 1): " << matrix.get(0, 1) << "\n";
    cout << "set element at (row = 1, col = 1) = 23: \n";
    matrix.set(1, 1, 23);
    matrix.print();
    cout << "==== done\n";
}

template<typename T>
void matrix_get_row_col_test() {
    cout << "\n==== test matrix get row and col: \n";
    cout << "matrix row:3 & col: 4\n";
    matrix::matrix<T> matrix(3, 4, matrix::One);
    matrix.set(2, 1, 12);
    matrix.print();
    cout << "get matrix row(2):\n";
    vector<T> row = matrix.row(2);
    for (int i = 0; i < row.size(); i++) {
        cout << row[i] << " ";
    }
    cout << "\n";
    cout << "get matrix col(1):\n";
    vector<T> col = matrix.col(1);
    for (int i = 0; i < col.size(); i++) {
        cout << col[i] << " ";
    }
    cout << "\n";
    cout << "print matrix row(2):\n";
    matrix.print_row(2);
    cout << "print matrix col(2):\n";
    matrix.print_col(2);
    cout << "==== done\n";
}

template<typename T>
void matrix_calculation_test() {
    cout << "\n==== test matrix get row and col: \n";
    cout << "matrix addition:\n";
    matrix::matrix<T> matrix1(2, 2, matrix::Identity);
    matrix::matrix<T> matrix2(2, 2, matrix::NegaOne);
    matrix1.print();
    cout << "+\n";
    matrix2.print();
    cout << "=\n";
    matrix1 += matrix2;
    matrix1.print();

    cout << "matrix subtraction:\n";
    matrix1.resize(2, 2, matrix::Identity);
    matrix2.resize(2, 2, matrix::NegaOne);
    matrix1.print();
    cout << "-\n";
    matrix2.print();
    cout << "=\n";
    matrix1 -= matrix2;
    matrix1.print();

    cout << "matrix and number multiplication:\n";
    matrix1.resize(2, 2, matrix::Identity);
    matrix1.print();
    cout << "*\n";
    cout << "3\n";
    cout << "=\n";
    matrix1 *= 3;
    matrix1.print();

    cout << "matrix and matrix multiplication:\n";
    matrix1.resize(2, 3, matrix::One);
    matrix2.resize(3, 3, matrix::One);
    matrix1.print();
    cout << "*\n";
    matrix2.print();
    cout << "=\n";
    matrix1.mult(matrix2);
    matrix1.print();
    cout << "==== done\n";
}

template<typename T>
void matrix_transform_test() {
    matrix::matrix<T> matrix(3, 2, 0.0,
                             1, 9,
                             4, 5,
                             7, 3);
    matrix.print();
    matrix = matrix.transform();
    matrix.print();
}

template<typename T>
void matrix_adjoint_matrix_test() {
    matrix::matrix<T> matrix(3, 3, 0,
                             1, 2, 3,
                             4, 4, 6,
                             7, 8, 9);
    matrix.print();
    cout << matrix.det() << "\n";
    matrix::matrix<T> matrix1 = matrix.adjoint();
    matrix1.print();
    cout << matrix.is_invertible() << "\n";
    matrix::matrix<T> matrix2 = matrix.inverse();
    matrix2.print();
    cout << matrix.is_symmetric() << "\n";
}

template<typename T>
void matrix_block_test() {
    matrix::matrix<T> matrix(5, 5, 0,
                             1, 2, 3, 4, 5,
                             6, 7, 8, 9, 0,
                             0, 9, 8, 7, 6,
                             5, 4, 3, 2, 1,
                             1, 2, 3, 4, 5);
    matrix.print();
    matrix::matrix<T> m11;
    matrix::matrix<T> m12;
    matrix::matrix<T> m21;
    matrix::matrix<T> m22;
    matrix.partitioned4(m11, m12, m21, m22);
    m11.print();
    m12.print();
    m21.print();
    m22.print();
}

template<typename T>
void matrix_test() {
    matrix_constructor_test<T>();
    matrix_getter_setter_test<T>();
    matrix_get_row_col_test<T>();
    matrix_calculation_test<T>();
    matrix_transform_test<T>();
    matrix_adjoint_matrix_test<T>();
    matrix_block_test<T>();
}

int main() {
    matrix_test<int>();
    return 0;
}
