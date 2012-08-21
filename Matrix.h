#pragma once

template<class T>
class Matrix
{
  T* data;
  int rows, cols;

public:
  Matrix(int r, int c) : rows(r), cols(c) 
  {
    data = new T[r * c];
  }
  Matrix(int r) : rows(r), cols(r)
  {
    data = new T[r * r];
  }
  ~Matrix()
  {
    delete [] data;
  }
  T& operator()(int r, int c)
  {
    return data[c + r * cols];
  }
  const T& operator()(int r, int c) const
  {
    return data[c + r * cols];
  }
  int GetRows() const
  {
    return rows;
  }
  int GetCols() const
  {
    return cols;
  }
  void FillWithValue(const T& val)
  {
    fill(data, data + (rows * cols), val);
  }
  int GetMinRowElementColumn(int row, int start, int stop) const
  {
    T min = (*this)(row, start);
    int col = start;
    for (int i = start + 1; i < stop; ++i)
    {
      if ((*this)(row, i) < min)
      {
        min = (*this)(row, i);
        col = i;
      }
    }
    return col;
  }
};
