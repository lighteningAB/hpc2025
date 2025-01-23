#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

std::vector<double> vecadd(std::vector<double> a, std::vector<double> b)
{
    // allocating dynamic memory
    std::vector<double> returnable;
    for (int i = 0; i < a.size(); i++)
    {
        returnable.push_back(a.at(i) + b.at(i));
    }
    return returnable;
}

std::vector<double> dot(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> returnable;
    for (int i = 0; i < a.size(); i++)
    {
        returnable.push_back(a.at(i) * b.at(i));
    }
    return returnable;
}

std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b){
    if (a[0].size() != b.size()){
        std::cout << "Error: matrices are not compatible for multiplication" << std::endl;
        exit(1);
    }
    std::vector<std::vector<double>> returnable;
    // iterate over the rows of a
    for (int i = 0; i < a.size(); i++){
        std::vector<double> row = a[i];
        std::vector<double> addrow;
        //access each column of b
        for (int j =0; j<b[0].size(); j++){
            double sum = 0;
            //iterate down the column of b
            for (int k = 0; k < b.size(); k++){
                sum += row[k] * b[k][j];
            }
            addrow.push_back(sum);
        }
        returnable.push_back(addrow);
    }
    return returnable;
}
void print_vector(const std::vector<double>& vec) {
    for (double val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void print_matrix(const std::vector<std::vector<double>>& mat) {
    for (const auto& row : mat) {
        print_vector(row);
    }
}
int main()
{
    std::vector<double>* vec1 = new std::vector<double>{1.0, 2.0, 3.0};
    std::vector<double>* vec2 = new std::vector<double>{4.0, 5.0, 6.0};
    std::vector<double> result_vecadd = vecadd(*vec1, *vec2);
    std::cout << "Result of vecadd: ";
    print_vector(result_vecadd);

    // Test dot
    std::vector<double> result_dot = dot(*vec1, *vec2);
    std::cout << "Result of dot: ";
    print_vector(result_dot);

    // Test multiply
    std::vector<std::vector<double>>* mat1 = new std::vector<std::vector<double>>{
        {1.0, 2.0},
        {3.0, 4.0}
    };
    std::vector<std::vector<double>>* mat2 = new std::vector<std::vector<double>>{
        {5.0, 6.0},
        {7.0, 8.0}
    };
    std::vector<std::vector<double>> result_multiply = multiply(*mat1, *mat2);
    std::cout << "Result of multiply: " << std::endl;
    print_matrix(result_multiply);

    // Clean up dynamic memory
    delete vec1;
    delete vec2;
    delete mat1;
    delete mat2;
    return 0;
}