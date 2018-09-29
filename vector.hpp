#ifndef VECTOR
#define VECTOR

#include<iostream>
#include<initializer_list>
#include<cmath>
#include<stdexcept>
#include<string>
#include<typeinfo>

template<class T, int size>
class Vector
{
private:
	T v[size] = { };
public:
// Constructors ----------------------------------------------------//
	Vector() { };

	Vector(const Vector &vector);

	Vector(const T *list);

	Vector(const std::initializer_list<T> &list);


// Accessing Indices -----------------------------------------------//
	const T& operator[](int index) const { return v[index]; }


// Explicit Type Conversions ---------------------------------------//
	template<class U>
	explicit operator Vector<U, size>() const
	{
		U *t0 = new U[size];
		for(int i = 0; i < size; i++)
			t0[i] = (U)v[i];
		Vector<U, size> t1(t0);
		delete[] t0;
		return t1;
	}
	
	template<int sizeU>
	explicit operator Vector<T, sizeU>() const
	{
		T *t0 = new T[sizeU];
		for(int i = 0; i < sizeU; i++)
			t0[i] = i < size ? v[i] : 0;
		Vector<T, sizeU> t1(t0);
		delete[] t0;
		return t1;
	}

	template<class U, int sizeU>
	explicit operator Vector<U, sizeU>() const
	{
		U *t0 = new U[sizeU];
		for(int i = 0; i < sizeU; i++)
			t0[i] = (U)(i < size ? v[i] : 0);
		Vector<U, sizeU> t1(t0);
		delete[] t0;
		return t1;
	}


// Addition --------------------------------------------------------//
	Vector operator+(const Vector &other) const;
	Vector operator+() const { return Vector(*this); }
	Vector& operator+=(const Vector &other);


// Subtraction -----------------------------------------------------//
	Vector operator-(const Vector &other) const;
	Vector operator-() const;
	Vector& operator-=(const Vector &other);


// Scalar Multiplication -------------------------------------------//	
	Vector operator*(const T &m) const;
	Vector& operator*=(const T &m);


// Scalar Division -------------------------------------------------//
	Vector operator/(const T &d) const;
	Vector& operator/=(const T &d);


// Scalar Product --------------------------------------------------//
	T dot(const Vector &other) const;


// Cross Product ---------------------------------------------------//
	T cross(const Vector<T, 2> &other) const;

	Vector<T, 3> cross(const Vector<T, 3> &other) const;


// Vector Length ---------------------------------------------------//
	float norm() const;

	float norm2() const;


// Unit Vector -----------------------------------------------------//
	Vector unit() const;
};


// Constructors ----------------------------------------------------//
template<class T, int size>
Vector<T, size>::Vector(const Vector &vector)
{
	for (int i = 0; i < size; i++)
		v[i] = vector[i];
}

template<class T, int size>
Vector<T, size>::Vector(const T *list)
{
	for(int i = 0; i < size; i++)
		v[i] = list[i];
}

template<class T, int size>
Vector<T, size>::Vector(const std::initializer_list<T> &list)
{
	int i = 0;
	for (const T &e : list)
		v[i++] = e;
}
//------------------------------------------------------------------//


// Addition --------------------------------------------------------//
template<class T, int size>
Vector<T, size> Vector<T, size>::operator+(
		const Vector &other) const
{
	Vector t;
	for (int i = 0; i < size; i++)
		t.v[i] = v[i] + other[i];
	return t;
}

template<class T, int size>
Vector<T, size>& Vector<T, size>::operator+=(
		const Vector &other)
{
	for (int i = 0; i < size; i++)
		v[i] += other[i];
	return *this;
}
//------------------------------------------------------------------//


// Subtraction -----------------------------------------------------//
template<class T, int size>
Vector<T, size> Vector<T, size>::operator-(
		const Vector &other) const
{
	Vector t;
	for (int i = 0; i < size; i++)
		t.v[i] = v[i] - other[i];
	return t;
}

template<class T, int size>
Vector<T, size> Vector<T, size>::operator-() const
{
	Vector t;
	for (int i = 0; i < size; i++)
		t.v[i] = -v[i];
	return t;
}

template<class T, int size>
Vector<T, size>& Vector<T, size>::operator-=(
		const Vector &other)
{
	for (int i = 0; i < size; i++)
		v[i] -= other[i];
	return *this;
}
//------------------------------------------------------------------//


// Scalar Multiplication -------------------------------------------//
template<class U, class T, int size>
Vector<T, size> operator*(const U &m0,
		const Vector<T, size> &vector)
{
	T *t0 = new T[size];
	T m1 = (T)m0;
	for(int i = 0; i < size; i++)
		t0[i] = m1 * vector[i];
	Vector<T, size> t1(t0);
	delete[] t0;
	return t1;
}

template<class T, int size>
Vector<T, size> Vector<T, size>::operator*(const T &m) const
{
	Vector t;
	for (int i = 0; i < size; i++)
		t.v[i] = v[i] * m;
	return t;
}

template<class T, int size>
Vector<T, size>& Vector<T, size>::operator*=(const T &m)
{
	for (int i = 0; i < size; i++)
		v[i] *= m;
	return *this;
}
//------------------------------------------------------------------//


// Scalar Division -------------------------------------------------//
template<class U, class T, int size>
Vector<T, size> operator/(const U &d0,
		const Vector<T, size> &vector)
{
	T *t0 = new T[size];
	T d1 = (T)d0;
	for(int i = 0; i < size; i++)
		t0[i] = d1 / vector[i];
	Vector<T, size> t1(t0);
	delete[] t0;
	return t1;
}

template<class T, int size>
Vector<T, size> Vector<T, size>::operator/(const T &d) const
{
	Vector t;
	for (int i = 0; i < size; i++)
		t.v[i] = v[i] / d;
	return t;
}

template<class T, int size>
Vector<T, size>& Vector<T, size>::operator/=(const T &d)
{
	for (int i = 0; i < size; i++)
		v[i] /= d;
	return *this;
}
//------------------------------------------------------------------//


// Scalar Product --------------------------------------------------//
template<class T, int size>
T Vector<T, size>::dot(const Vector &other) const
{
	T t = 0;
	for (int i = 0; i < size; i++)
		t += v[i] * other[i];
	return t;
}
//------------------------------------------------------------------//


// Cross Product ---------------------------------------------------//
template<class T, int size>
T Vector<T, size>::cross(const Vector<T, 2> &other) const
{
	if (size != 2)
		throw std::domain_error(std::string("Vector cross product
				not defined for Vector<")
				+ typeid(T).name() + ", " + std::to_string(size)
				+ ">, Vector<" + typeid(T).name() + ", 2>");
	return v[0] * other[1] - v[1] * other[0];
}

template<class T, int size>
Vector<T, 3> Vector<T, size>::cross(const Vector<T, 3> &other) const
{
	if (size != 3)
		throw std::domain_error(std::string("Vector cross product
				not defined for Vector<")
				+ typeid(T).name() + ", " + std::to_string(size)
				+ ">, Vector<" + typeid(T).name() + ", 3>");
	return Vector<T, 3> { v[1] * other[2] - v[2] * other[1],
						  v[2] * other[0] - v[0] * other[2],
						  v[0] * other[1] - v[1] * other[0] };
}
//------------------------------------------------------------------//


// Vector Length ---------------------------------------------------//
template<class T, int size>
float Vector<T, size>::norm() const
{
	float t = 0;
	for(int i = 0; i < size; i++)
		t += v[i] * v[i];
	return sqrt(t);
}

template<class T, int size>
float Vector<T, size>::norm2() const
{
	float t = 0;
	for(int i = 0; i < size; i++)
		t += v[i] * v[i];
	return t;
}
//------------------------------------------------------------------//


// Unit Vector -----------------------------------------------------//
template<class T, int size>
Vector<T, size> Vector<T, size>::unit() const
{
	return (*this) / norm();
}
//------------------------------------------------------------------//


// Printing Vector -------------------------------------------------//
template<class T, int size>
std::ostream& operator<<(std::ostream &out,
		const Vector<T, size> &vector)
{
	out << "(";
	for(int i = 0; i < size; i++)
		out << vector[i] << (i == size - 1 ? ")" : ", ");
	return out;
}
//------------------------------------------------------------------//


#endif //VECTOR
