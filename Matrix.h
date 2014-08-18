#pragma once
#include "Vector.h"
#include <algorithm>

class InvalidRowException
{
	virtual const char* what() const throw()
	{
		return "Tried to access invalid row";
	}
};

class InvalidColumnException
{
	virtual const char* what() const throw()
	{
		return "Tried to access invalid column";
	}
};

template <typename CT, int R, int C>
class Matrix
{
	template <typename, int> friend class Vector;
	template <typename, int, int> friend class Matrix;

private:
	CT _Data[R][C];
public:
	Matrix(const Vector<CT, C>(&data)[R])
	{
		for (int i = 0; i < R; i++)
		{
			memcpy(_Data[i], data, sizeof(CT) * C);
		}
	}

	Matrix(const CT(&data)[R][C])
	{
		for (int i = 0; i < R; i++)
		{
			memcpy(_Data[i], data[i], sizeof(CT) * C);
		}
	}

	Matrix()
	{
		for (int i = 0; i < R; i++)
		{
			memset(_Data[i], 0, sizeof(CT) * C);
		}
	}

	static Matrix<CT, R, C> Identity()
	{
		Matrix<CT, R, C> result;
		static_assert(R == C, "Cannot create a non-square identity matrix");
		for (int n = 0; n < R; n++)
		{
			result._Data[n][n] = 1;
		}
		return result;
	}

	//Type safe accessors
	template<int RI, int CI>
	CT At()
	{
		static_assert(RI >= 0 && RI < R, "Specified row exceeds matrix bounds");
		static_assert(CI >= 0 && CI < C, "Specified column exceeds matrix bounds");
		return _Data[RI][CI];
	}

	//Type safe accessors
	template<int RI, int CI>
	Matrix<CT, R, C> At(CT n)
	{
		static_assert(RI >= 0 && RI < R, "Specified row exceeds matrix bounds");
		static_assert(CI >= 0 && CI < C, "Specified column exceeds matrix bounds");
		CT data[R][C];
		memcpy(data, _Data, sizeof(CT) * R * C);
		data[RI][CI] = n;
		return Matrix<CT, R, C>(data);
	}

	template<int RI>
	Vector<CT, C> Row()
	{
		static_assert(RI >= 0 && RI < R, "Specified row exceeds matrix bounds");
		return Vector<CT, C>(_Data[RI]);
	}

	template<int CI>
	Vector<CT, R> Column()
	{
		static_assert(CI >= 0 && CI < C, "Specified column exceeds matrix bounds");
		CT data[R];
		for (int i = 0; i < R; i++)
		{
			data[i] = _Data[i][CI];
		}
		return data;
	}

	Vector<CT, (R < C) ? R : C> Diagonal()
	{
		const int size = (R < C) ? R : C;
		CT data[size];
		for (int i = 0; i < size; i++)
		{
			data[i] = _Data[i][i];
		}
		return Vector<CT, size>(data);
	}

	Vector<CT, (R < C) ? R : C> AntiDiagonal()
	{
		const int size = (R < C) ? R : C;
		CT data[size];
		for (int i = 0; i < size; i++)
		{
			data[i] = _Data[size - i - 1][i];
		}
		return Vector<CT, size>(data);
	}

	template<int RI>
	Matrix<CT, R - 1, C> WithoutRow()
	{
		static_assert(RI >= 0 && RI < R, "Specified row exceeds matrix bounds");
		CT data[R - 1][C];
		int index = 0;
		for (int r = 0; r < R; r++)
		{
			if (r != RI) 
			{
				for (int c = 0; c < C; c++)
				{
					data[index][c] = _Data[r][c];
				}
				index++;
			}
		}
		return Matrix<CT, R - 1, C>(data);
	}

	template<int CI>
	Matrix<CT, R, C - 1> WithoutColumn()
	{
		static_assert(CI >= 0 && CI < C, "Specified column exceeds matrix bounds");
		CT data[R][C - 1];
		int index = 0;
		for (int c = 0; c < C; c++)
		{
			if (c != CI)
			{
				for (int r = 0; r < R; r++)
				{
					data[r][index] = _Data[r][c];
				}
				index++;
			}
		}
		return Matrix<CT, R, C - 1>(data);
	}

	template<int RI1, int RI2>
	Matrix<CT, R, C> SwapRows()
	{
		static_assert(RI1 >= 0 && RI2 < R, "First specified row exceeds matrix bounds");
		static_assert(RI2 >= 0 && RI2 < R, "Second specified row exceeds matrix bounds");
		CT data[R][C];
		for (int r = 0; r < R; r++)
		{
			int rIndex = (r == RI1) ? RI2 : (r == RI2) ? RI1 : r;
			for (int c = 0; c < C; c++)
			{
				data[r][c] = _Data[rIndex][c];
			}
		}
		return Matrix<CT, R, C>(data);
	}

	template<int CI1, int CI2>
	Matrix<CT, R, C> SwapColumns()
	{
		static_assert(CI1 >= 0 && CI2 < C, "First specified column exceeds matrix bounds");
		static_assert(CI2 >= 0 && CI2 < C, "Second specified column exceeds matrix bounds");
		CT data[R][C];
		for (int c = 0; c < C; c++)
		{
			int cIndex = (c == CI1) ? CI2 : (c == CI2) ? CI1 : c;
			for (int r = 0; r < R; r++)
			{
				data[r][c] = _Data[r][cIndex];
			}
		}
		return Matrix<CT, R, C>(data);
	}

	CT At(int ri, int ci)  
	{
		if (!(ri >= 0 && ri < R)) { throw InvalidRowException(); }
		if (!(ci >= 0 && ci < C)) { throw InvalidColumnException(); }
		return _Data[ri][ci];
	}

	Matrix<CT, R, C> At(int ri, int ci, CT value)
	{
		if (!(ri >= 0 && ri < R)) { throw InvalidRowException(); }
		if (!(ci >= 0 && ci < C)) { throw InvalidColumnException(); }
		CT data[R][C];
		memcpy(data, _Data, sizeof(CT) * R * C);
		data[ri][ci] = value;
		return Matrix<CT, R, C>(data);
	}

	Vector<CT, C> Row(int ri)
	{
		if (!(ri >= 0 && ri < R)) { throw InvalidRowException(); }
		return Vector<CT, C>(_Data[ri]);
	}

	Vector<CT, R> Column(int ci)
	{
		if (!(ci >= 0 && ci < C)) { throw InvalidColumnException(); }
		CT data[R];
		for (int i = 0; i < R; i++)
		{
			data[i] = _Data[i][ci];
		}
		return data;
	}

	Matrix<CT, R - 1, C> WithoutRow(int ri)
	{
		if (!(ri >= 0 && ri < R)) { throw InvalidRowException(); }
		CT data[R - 1][C];
		int index = 0;
		for (int r = 0; r < R; r++)
		{
			if (r != ri)
			{
				for (int c = 0; c < C; c++)
				{
					data[index][c] = _Data[r][c];
				}
				index++;
			}
		}
		return Matrix<CT, R - 1, C>(data);
	}

	Matrix<CT, R, C - 1> WithoutColumn(int ci)
	{
		if (!(ci >= 0 && ci < C)) { throw InvalidColumnException(); }
		CT data[R][C - 1];
		int index = 0;
		for (int c = 0; c < C; c++)
		{
			if (c != ci)
			{
				for (int r = 0; r < R; r++)
				{
					data[r][index] = _Data[r][c];
				}
				index++;
			}
		}
		return Matrix<CT, R, C - 1>(data);
	}

	Matrix<CT, R, C> SwapRows(int ri1, int ri2)
	{
		if (!(ri1 >= 0 && ri2 < R)) { throw InvalidRowException(); }
		if (!(ri2 >= 0 && ri2 < R)) { throw InvalidRowException(); }
		CT data[R][C];
		for (int r = 0; r < R; r++)
		{
			int rIndex = (r == ri1) ? ri2 : (r == ri2) ? ri1 : r;
			for (int c = 0; c < C; c++)
			{
				data[r][c] = _Data[rIndex][c];
			}
		}
		return Matrix<CT, R, C>(data);
	}

	Matrix<CT, R, C> SwapColumns(int ci1, int ci2)
	{
		if (!(ci1 >= 0 && ci2 < C)) { throw InvalidColumnException(); }
		if (!(ci2 >= 0 && ci2 < C)) { throw InvalidColumnException(); }
		CT data[R][C];
		for (int c = 0; c < C; c++)
		{
			int cIndex = (c == ci1) ? ci2 : (c == ci2) ? ci1 : c;
			for (int r = 0; r < R; r++)
			{
				data[r][c] = _Data[r][cIndex];
			}
		}
		return Matrix<CT, R, C>(data);
	}

	Matrix<CT, C, R> Transpose()
	{
		CT data[C][R];

		for (int r = 0; r < R; r++)
		{
			for (int c = 0; c < C; c++)
			{
				data[c][r] = _Data[r][c];
			}
		}

		return Matrix<CT, C, R>(data);
	}

	template<int C2>
	Matrix<CT, R, C2> operator*(Matrix<CT, C, C2> other)
	{
		CT data[R][C2];

		for (int r = 0; r < R; r++)
		{
			for (int c = 0; c < C2; c++)
			{
				data[r][c] = Row(r) & other.Column(c);
			}
		}

		return Matrix<CT, R, C2>(data);
	}

	Vector<CT, R> operator*(Vector<CT, R> other)
	{
		float data[R][1];
		for (int n = 0; n < R; n++)
		{
			data[n][0] = other.At(n);
		}

		return ((*this) * Matrix<CT, R, 1>(data)).Column<0>();
	}

	Matrix<CT, R, C> operator+(Matrix<CT, R, C> other)
	{
		CT data[R][C2];

		for (int r = 0; r < R; r++)
		{
			for (int c = 0; c < C2; c++)
			{
				data[r][c] = _Data[r][c] + other._Data[r][c];
			}
		}

		return Matrix<CT, R, C2>(data);
	}

	Matrix<CT, R, C> operator-(Matrix<CT, R, C> other)
	{
		CT data[R][C2];

		for (int r = 0; r < R; r++)
		{
			for (int c = 0; c < C2; c++)
			{
				data[r][c] = _Data[r][c] - other._Data[r][c];
			}
		}

		return Matrix<CT, R, C2>(data);
	}

	Vector<CT, R> GaussianElimination(Vector<CT, R> values)
	{
		Matrix<CT, R, R> mat(_Data);
		Matrix<CT, R, 1> solutions = Matrix<CT, 1, R>({ values }).Transpose();
		Matrix<CT, R, 1> results;

		for (int diagonal = 0; diagonal < R; diagonal++)
		{
			int pivot = diagonal;
			CT pivotValue = mat.At(diagonal, diagonal);
			for (int row = diagonal + 1; row < R; row++)
			{
				CT t = abs(mat.At(row, diagonal));
				if (t > pivotValue)
				{
					pivotValue = t;
					pivot = row;
				}
			}	

			mat = mat.SwapRows(diagonal, pivot);
			solutions = solutions.SwapRows(diagonal, pivot);

			for (int row = diagonal + 1; row < R; row++)
			{
				CT ratio = mat.At(row, diagonal) / mat.At(diagonal, diagonal); 
				for (int col = diagonal + 1; col < R; col++)
				{
					mat = mat.At(row, col, mat.At(row, col) - ratio * mat.At(diagonal, col));
				}
				mat = mat.At(row, diagonal, CT());
				solutions = solutions.At(row, 0, solutions.At(row, 0) - ratio * solutions.At(diagonal, 0));
			}
		}

		for (int row = R - 1; row >= 0; row--)
		{
			CT t = solutions.At(row, 0);
			for (int v = R - 1; v > row; v--)
			{
				t -= results.At(v, 0) * mat.At(row, v);
			}
			results = results.At(row, 0, t / mat.At(row, row));
		}

		return Vector<CT, R>(results.Transpose()._Data[0]);
	}

	int Rows()
	{
		return R;
	}

	int Columns()
	{
		return C;
	}

	bool Square()
	{
		return R == C;
	}
};

template<typename CT, int R, int C>
std::ostream& operator<<(std::ostream &strm, Matrix<CT, R, C> &m)
{
	for (int i = 0; i < R; i++)
	{
		strm << "[" << m.Row(i) << "]" << std::endl;
	}
	return strm;
}