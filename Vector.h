#pragma once
#include <memory>
#include <iostream>

class InvalidVectorElementException
{
	virtual const char* what() const throw()
	{
		return "Invalid vector element";
	}
};

template <typename CT, int S>
class Vector
{
	template <typename, int> friend class Vector;

private:
	CT _Data[S];

public:
	//Named element getters
	CT X() { static_assert(S >= 1, "Vector must have 1 or more elements to access X"); return _Data[0]; }
	CT Y() { static_assert(S >= 2, "Vector must have 2 or more elements to access Y"); return _Data[1]; }
	CT Z() { static_assert(S >= 3, "Vector must have 3 or more elements to access Z"); return _Data[2]; }
	CT W() { static_assert(S >= 4, "Vector must have 4 or more elements to access W"); return _Data[3]; }

	Vector<CT, S> X(CT n) { static_assert(S >= 1, "Vector must have 1 or more elements to access X"); CT data[S]; memcpy(data, _Data, sizeof(CT)* S); data[0] = n; return Vector<CT, S>(data); }
	Vector<CT, S>  Y(CT n) { static_assert(S >= 2, "Vector must have 2 or more elements to access Y"); CT data[S]; memcpy(data, _Data, sizeof(CT)* S); data[1] = n; return Vector<CT, S>(data); }
	Vector<CT, S>  Z(CT n) { static_assert(S >= 3, "Vector must have 3 or more elements to access Z"); CT data[S]; memcpy(data, _Data, sizeof(CT)* S); data[2] = n; return Vector<CT, S>(data); }
	Vector<CT, S>  W(CT n) { static_assert(S >= 4, "Vector must have 4 or more elements to access W"); CT data[S]; memcpy(data, _Data, sizeof(CT)* S); data[3] = n; return Vector<CT, S>(data); }

	//Static element getter (type-safe)
	template <int N> CT At() { static_assert(N < S, "Vector must have N or more elements to access N"); return _Data[N]; }
	template <int N> CT At(CT n) { static_assert(N < S, "Vector must have N or more elements to access N"); CT data[S]; memcpy(data, _Data, sizeof(CT)* S); data[N] = n; return Vector<CT, S>(data); }

	//Dynamic element getter (bounds-checked at runtime)
	CT At(int n) { if (n >= 0 && n < S) { return _Data[n]; } throw InvalidVectorElementException(); }

	Vector() { ZeroData(); }
	Vector(Vector<CT, S> &other) { memcpy(_Data, other._Data, sizeof(CT) * S); } //Type safe constructor from an identically statically-sized vector
	Vector(const CT(&data)[S]) { ZeroData();  memcpy(_Data, data, sizeof(CT)* S); } //Type safe constructor from a fixed const-sized array

	static Vector<CT, S> From(CT* data) { Vector<CT, S> newVec; memcpy(newVec._Data, data, sizeof(CT)* S); return newVec;  } //Unsafe dynamic constructor. Cannot be bounds-checked, so try to avoid using it!

	//Vector-Vector arithmetic
#define VECTOR_OP(FNAME, OP) Vector<CT, S> FNAME (Vector<CT, S> other) { CT __result[S]; for(int i=0;i<S;i++) { CT x = _Data[i]; CT y = other._Data[i]; __result[i] = OP; } return Vector<CT, S>(__result); }
	VECTOR_OP(operator+, x + y)
	VECTOR_OP(operator-, x - y)
	VECTOR_OP(operator*, x * y)
	VECTOR_OP(operator/, x / y)
#undef VECTOR_OP

	//Vector-Scalar arithmetic
#define SCALAR_OP(FNAME, OP) Vector<CT, S> FNAME (CT other) { CT __result[S]; for(int i=0;i<S;i++) { CT x = _Data[i]; CT y = other; __result[i] = OP; } return Vector<CT, S>(__result); }
	SCALAR_OP(operator*, x * y)
	SCALAR_OP(operator/, x / y)
#undef SCALAR_OP

	CT Sum()
	{
		CT r = 0;
		for (int i = 0; i < S; i++)
		{
			r += _Data[i];
		}
		return r;
	}

	CT Product()
	{
		CT r = S == 0 ? 0 : 1;
		for (int i = 0; i < S; i++)
		{
			r *= _Data[i];
		}
		return r;
	}

	CT MagnitudeSquared()
	{
		return (*this * *this).Sum();
	}

	CT Magnitude()
	{
		return sqrt(MagnitudeSquared());
	}

	Vector<CT, S> Normalized()
	{
		return (*this) / Magnitude();
	}

	static CT Dot(Vector<CT, S> first, Vector<CT, S> second)
	{
		return (first * second).Sum();
	}

	static Vector<CT, S> Cross(Vector<CT, S> first, Vector<CT, S> second)
	{
		static_assert(S == 3, "Cannot take cross-product of a non 3 dimensional vector");
		CT data[3];
		data[0] = first.At<1>() * second.At<2>() - first.At<2>() * second.At<1>();
		data[1] = first.At<2>() * second.At<0>() - first.At<0>() * second.At<2>();
		data[2] = first.At<0>() * second.At<1>() - first.At<1>() * second.At<0>();
		return Vector<CT, S>(data);
	}

	//No good dot product operators :-(
	CT operator&(Vector<CT, S> other)
	{
		return Dot(*this, other);
	}

	//No good cross product operators either :-(
	Vector<CT, S> operator|(Vector<CT, S> other)
	{
		return Cross(*this, other);
	}

private:
	void ZeroData()
	{
		memset(_Data, 0, sizeof(CT) * S);
	}
};

template<int S> using FVector = Vector<float, S>;
template<int S> using DVector = Vector<double, S>;

template<typename CT, int S>
std::ostream& operator<<(std::ostream &strm, Vector<CT, S> &v)
{
	strm << "(";
	for (int i = 0; i < S; i++)
	{
		strm << v.At(i) << ((i == (S - 1)) ? "" : ", ");
	}
	strm << ")";
	return strm;
}