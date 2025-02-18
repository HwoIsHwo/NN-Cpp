#include <iostream>
#include <fstream>
#include "matrix.h"



matrix::matrix(unsigned int rows_in, unsigned int columns_in) //constructor with parameters
{
    rows = rows_in;
    columns = columns_in;
    data = new float[rows*columns];
}



matrix::matrix(unsigned int rows_in, unsigned int columns_in, float* data_in)//constructor with parameters
{
    rows = rows_in;
    columns = columns_in;
    data = new float[rows*columns];

    int i;
    for(i=0; i<rows*columns; i++)
    {
        data[i] = data_in[i];
    }
}



matrix::matrix(const matrix &m) //copy constructor
{
    rows = m.rows;
    columns = m.columns;
    data = new float[rows*columns];

    int i;
    for(i=0; i<rows*columns; i++)
    {
        data[i] = m.data[i];
    }
}



void matrix::print() //matrix output
{
    int i;
    int ii;
    for(i=0; i<rows; i++)
    {
        for(ii=0; ii<this->columns; ii++)
        {
            std::cout << data[i*columns+ii] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}



void matrix::file_read(char* filename) //read data (not raws and columns amont) from file
{
    std::ifstream file;
    file.open(filename);

    if(file.is_open()) //if file is opened
    {
        int i;
        for(i=0; i<rows*columns; i++) //trying to read each element
        {
            if(file.eof()) break; // go out if it's the end of file
            file >> data[i];
        }
    }else
    {
        std::cout << "File wasn't opened" << std::endl;
    }

    file.close();
}



void matrix::file_write(char* filename) //write data in file
{
    std::ofstream file;
    file.open(filename);

    if(file.is_open()) //if file is opened
    {
        int i;
        for(i=0; i<rows*columns; i++) //write each element
        {
            file << data[i] << std::endl;
        }
    }else
    {
        std::cout << "File wasn't opened" << std::endl;
    }

    file.close();
}



void matrix::write_data(int i, int j, float d)
{
    data[i*columns+j] = d;
}



float matrix::read_data(int i, int j)
{
    return data[i*columns+j];
}



matrix& matrix::operator=(const matrix& m) //operator overloadng '=' matrix = matrix
{
    if(&m != this) //check for writing in itself
    {
        delete [] data; //delete obj
        data = nullptr;

        rows = m.rows;
        columns = m.columns;
        data = new float[rows*columns];

        int i;
        for(i=0; i<rows*columns; i++) //copy data
        {
            data[i] = m.data[i];
        }
    }else
    {
        std::cout << "Rewriting matrix by itself" << std::endl;
    }

    return *this;
}



matrix matrix::operator+(const matrix& m) //operator overloadng '+' matrix + matrix
{
    matrix A{*this}; //output matrix

    if(rows==m.rows && columns==m.columns) //check for dimensional conformity
    {
        int i;
        for(i=0; i<rows*columns; i++)
        {
            A.data[i] += m.data[i];
        }
    }else
    {
        std::cout << "Dimensions aren't equal" << std::endl;
    }

    return A;
}



matrix matrix::operator+(float n) //operator overloadng '+' matrix + number
{
    matrix A{*this}; //output matrix

    int i;
    for(i=0; i<rows*columns; i++)
    {
        A.data[i] += n;
    }

    return A;
}



matrix matrix::operator*(const matrix &m) //operator overloadng '*' matrix * matrix
{
    matrix A(rows, m.columns); //output matrix

    if(columns == m.rows) //check for dimensional conformity
    {
        int i;
        int ii;
        int iii;
        float element;
        for(i=0; i<A.rows; i++) //for each row of A matrix
        {
            for(ii=0; ii<A.columns; ii++) //for each element in row
            {
                element = 0;
                for(iii=0; iii<columns; iii++) //calculation of product i-row of left matrix and ii-column of right matrix
                {
                    element = element + data[i*columns+iii] * m.data[iii*m.columns+ii];
                }
                A.data[i*A.columns+ii] = element;
            }
        }
    }else
    {
        std::cout << "Columns amount of left matrix != Rows amount of right matrix" << std::endl;
    }

    return A;
}



matrix matrix::operator*(float n) //operator overloadng '*' matrix * number
{
    matrix A{*this}; //output matrix

    int i;
    for(i=0; i<rows*columns; i++) //copy data
    {
        A.data[i] *= n;
    }

    return A;
}



matrix matrix::operator-(const matrix& m) //operator overloadng '-' matrix - matrix
{
    matrix A{*this}; //output matrix

    if(rows==m.rows && columns==m.columns) //check for dimensional conformity
    {
        int i;
        for(i=0; i<rows*columns; i++)
        {
            A.data[i] -= m.data[i];
        }
    }else
    {
        std::cout << "Dimensions aren't equal" << std::endl;
    }

    return A;
}



matrix matrix::operator-(float n) //operator overloadng '-' matrix - number
{
    matrix A{*this}; //output matrix

    int i;
    for(i=0; i<rows*columns; i++)
    {
        A.data[i] -= n;
    }

    return A;
}



matrix matrix::transpose() //matrix transpose
{
    matrix A(columns, rows); //output matrix

    int i;
    int ii;
    for(i=0; i<rows; i++)
    {
        for(ii=0; ii<columns; ii++)
        {
            A.data[ii*rows+i] = data[i*columns+ii];
        }
    }

    return A;
}



matrix matrix::minor(int row, int column) //minor matrix of element in 'row' and 'column'
{
    if(rows==columns) //matrix is square
    {
        matrix A(rows-1, columns-1);

        int i;
        int ii;
        int cnt = 0;
        for(i=0; i<rows; i++)
        {
            if(i==row) continue;

            for(ii=0; ii<columns; ii++)
            {
                if(ii==column) continue;

                A.data[cnt] = data[i*columns+ii];
                cnt++;
            }
        }

        return A;
    }else
    {
        std::cout << "Matrix isn't square" << std::endl;
        matrix A{*this};

        return A;
    }
}



float matrix::deter() //determinant of matrix
{
    float D = 0;

    if(rows==columns) //matrix is square
    {
        if(rows==2)//exit point if matrix is 2x2
        {
            D = data[0] * data[3] - data[2] * data[1];
            return D;
        }

        matrix S(1, columns); //matrix of sings [1 -1 1 -1 1 -1 ...]
        int a = -1;
        int i;
        for(i=0; i<columns; i++)
        {
            a = a * (-1);
            S.data[i] = a;
        }

        for(i=0; i<columns; i++) //determinant by the first row
        {
            D = D + S.data[i] * data[i] * minor(0, i).deter();
        }
    }else
    {
        std::cout << "Matrix isn't square" << std::endl;
    }

    return D;
}



matrix matrix::operator!() //inverse matrix
{
    matrix A{*this}; //output matrix

    if(rows==columns) //matrix is square
    {
        float d = deter();
        if(d!=0)
        {
            matrix M{*this}; //matrix of minors
            int i;
            int ii;
            for(i=0; i<rows; i++)
            {
                for(ii=0; ii<columns; ii++)
                {
                    M.data[i*columns+ii] = minor(i, ii).deter();
                }
            }

            matrix S{*this}; //matrix of sings
            int a;
            int b = 1;
            for(i=0; i<rows; i++)
            {
                a = 1;
                for(ii=0; ii<columns; ii++)
                {
                    S.data[i*columns+ii] = a * b;
                    a = a * (-1);
                }
                b = b * (-1);
            }

            for(i=0; i<rows; i++) //matrix of algebraic complements
            {
                for(ii=0; ii<columns; ii++)
                {
                    A.data[i*columns+ii] = M.data[i*columns+ii] * S.data[i*columns+ii];
                }
            }

            A = A.transpose();
            A = A * (1 / d);
        }else
        {
            std::cout << "Determinant = 0, it's impossible to get iverse matrix" << std::endl;
        }
    }else
    {
        std::cout << "Matrix isn't square" << std::endl;
    }

    return A;
}



void matrix::eye()
{
    int i;
    int ii;
    for(i=0; i<rows; i++)
    {
        for(ii=0; ii<columns; ii++)
        {
            if(i==ii)
            {
                data[i*columns+ii] = 1;
            }else
            {
                data[i*columns+ii] = 0;
            }
        }
    }
}



matrix::~matrix() //destructor
{
    delete [] data;
    data = nullptr;
}
