/*=====================================================================================*/
/*! \file		utilities.h
	\author		peterkau
	\brief		Declaration of class Vector2T
 */
/*=====================================================================================*/

#ifndef UTILITIES_T_H
#define UTILITIES_T_H

#include <iostream>
#include <string>
#include <map>
using std::map;

//! Simple 2D vector class
template<typename T>
class Vector2T
{
public:
	//----------------------------------------------------------- constructors
	Vector2T() { m_xy[0] = (T)0; m_xy[1] = (T)0;}
	Vector2T(T v) { m_xy[0] = v; m_xy[1] = v; }
	Vector2T(const T x, const T y)  { m_xy[0] = x; m_xy[1] = y; }
	Vector2T(const Vector2T<T> &other) { *this = other; }
	
	//----------------------------------------------------------- assignment
	Vector2T<T> & operator=( const Vector2T<T> &other ) { 
		m_xy[0]= other.m_xy[0]; m_xy[1] = other.m_xy[1]; return *this; }

	//----------------------------------------------------------- element access
	//!  read-write
	T& operator[](int i) { assert(i>=0 && i<dim()); return m_xy[i]; }
	// NT
	const T& x() const { return m_xy[0]; }
	const T& y() const { return m_xy[1]; }

	//!  read-only
	const T& operator[](int i) const { assert(i>=0 && i<dim()); return m_xy[i]; }

	bool operator==(const Vector2T<T> &other) const
	{
		if(x()!=other.x()) return false;
		if(y()!=other.y()) return false;
		return true;
	}

	bool operator!=(const Vector2T<T> &other) const
	{
		return !(*this == other);
	}

	//----------------------------------------------------------- basic arithmetic operations
	//----------------------------------------------------------- scalar arithmetic operations
	Vector2T<T> operator*( T s ) const
	{
		return Vector2T<T>(m_xy[0]*s, m_xy[1]*s);
	}

	Vector2T<T>& operator*=( T s )
	{
		m_xy[0] = m_xy[0] * s;
		m_xy[1] = m_xy[1] * s;
		return *this;
	}

	Vector2T<T> operator/( T s ) const
	{
		return Vector2T<T>(m_xy[0]/s, m_xy[1]/s);
	}

	Vector2T<T>& operator/=( T s )
	{
		m_xy[0] = m_xy[0] / s;
		m_xy[1] = m_xy[1] / s;
		return *this;
	}

	//----------------------------------------------------------- vector arithmetic operations
	Vector2T<T> operator+( const Vector2T<T> &v2 ) const
	{
		return Vector2T<T>(x()+v2.x(),y()+v2.y());
	}

	Vector2T<T>& operator+=( const Vector2T<T> &v )
	{
		m_xy[0] += v.x();
		m_xy[1] += v.y();
		return *this;
	}

	Vector2T<T> operator-( const Vector2T<T> &v2 ) const
	{
		return Vector2T<T>(x()-v2.x(),y()-v2.y());
	}

	Vector2T<T>& operator-=( const Vector2T<T> &v )
	{
		m_xy[0] -= v.x();
		m_xy[1] -= v.y();
		return *this;
	}

	Vector2T<T> operator-() const
	{
		return Vector2T<T>(-m_xy[0], -m_xy[1]);
	}

	// inner product (dot product)
	T operator|(const Vector2T<T> &v2) const
	{
		return (*this)[0] * v2[0] + (*this)[1] * v2[1];
	}

	//----------------------------------------------------------- other
	Vector2T<T> normalized() const
	{
		return *this / length();
	}

	T norm() const
	{
		return length();
	}

	T sqrnorm() const
	{
		return squaredLength();
	}

	// NT
	T length() const
	{
		return sqrt(x()*x() + y()*y());
	}

	T squaredLength() const
	{
		return (x()*x() + y()*y());
	}

	static int dim() { return 2; }

private:
	T m_xy[2];
};

// scalar * vector
template<class T>
Vector2T<T> operator*(T const &s, Vector2T<T> const &v)
{
	return v*s;
}

// helper for printing a vector
template<class T>
std::ostream& operator<<( std::ostream& os, const  Vector2T<T> &v2 )
{
	char buf[256];
	//snprintf(buf,256,"[%f,%f]", (double)v2[0],(double)v2[1]);
	sprintf(buf,"[%f,%f]", (double)v2[0],(double)v2[1]);
	os << std::string(buf); 
	//os << '[' << i[0] << ", " << i[1] << ", " << i[2] << ']';
	return os;
}

//=====================================================================================*/

//! Simple 2x2 matrix class
template <typename Scalar>
class Matrix2x2T
{
public:
	//! Empty constructor - object remains uninitialized.
	Matrix2x2T() {}

	//! Constructor using matrix elements a00 ... a11
	Matrix2x2T(Scalar a00, Scalar a01, Scalar a10, Scalar a11) {
		m_elem[0][0] = a00; m_elem[0][1] = a01;
		m_elem[1][0] = a10; m_elem[1][1] = a11;
	}

	//! Copy-constructor
	Matrix2x2T(const Matrix2x2T &p) {
		m_elem[0][0] = p.m_elem[0][0];
		m_elem[0][1] = p.m_elem[0][1];
		m_elem[1][0] = p.m_elem[1][0];
		m_elem[1][1] = p.m_elem[1][1];
	}

	//! Returns 2
	static int rows() { return 2; }

	//! Returns 2
	static int cols() { return 2; }

	Scalar &operator() (int i, int j) {
		return m_elem[i][j];
	}

	const Scalar &operator() (int i, int j) const {
		return m_elem[i][j];
	}

	//! +=-Operator
	Matrix2x2T &operator+=(const Matrix2x2T &p) {
		m_elem[0][0] += p.m_elem[0][0];
		m_elem[0][1] += p.m_elem[0][1];
		m_elem[1][0] += p.m_elem[1][0];
		m_elem[1][1] += p.m_elem[1][1];
		return *this;
	}

	//! -=-Operator
	Matrix2x2T &operator-=(const Matrix2x2T &p) {
		m_elem[0][0] -= p.m_elem[0][0];
		m_elem[0][1] -= p.m_elem[0][1];
		m_elem[1][0] -= p.m_elem[1][0];
		m_elem[1][1] -= p.m_elem[1][1];
		return *this;
	}

	//! /=-Operator
	Matrix2x2T &operator/=(Scalar s) {
		m_elem[0][0] /= s;
		m_elem[0][1] /= s;
		m_elem[1][0] /= s;
		m_elem[1][1] /= s;
		return *this;
	}

	//! *=-Operator : Matrix * Scalar
	Matrix2x2T &operator*=(Scalar s) {
		m_elem[0][0] *= s;
		m_elem[0][1] *= s;
		m_elem[1][0] *= s;
		m_elem[1][1] *= s;
		return *this;
	}

	//! *=-Operator : Matrix * Matrix
	Matrix2x2T &operator*=(const Matrix2x2T &p) {
		*this = *this * p;
		return *this;
	}

	//! +-Operator
	Matrix2x2T operator+(const Matrix2x2T &p) const {
		return Matrix2x2T(
			m_elem[0][0] + p.m_elem[0][0],
			m_elem[0][1] + p.m_elem[0][1],
			m_elem[1][0] + p.m_elem[1][0],
			m_elem[1][1] + p.m_elem[1][1]);
	}

	//! --Operator
	Matrix2x2T operator-(const Matrix2x2T &p) const {
		return Matrix2x2T(
			m_elem[0][0] - p.m_elem[0][0],
			m_elem[0][1] - p.m_elem[0][1],
			m_elem[1][0] - p.m_elem[1][0],
			m_elem[1][1] - p.m_elem[1][1]);
	}

	//!-Operator
	Matrix2x2T operator/(Scalar s) const {
		return Matrix2x2T(
			m_elem[0][0] / s, m_elem[0][1] / s,
			m_elem[1][0] / s, m_elem[1][1] / s);
	}

	//! *-Operator : Matrix * Scalar
	Matrix2x2T operator*(Scalar s) const {
		return Matrix2x2T(
			m_elem[0][0] * s, m_elem[0][1] * s,
			m_elem[1][0] * s, m_elem[1][1] * s);
	}

	//! friend *-Operator : Scalar * Matrix
	friend Matrix2x2T operator*(Scalar s, const Matrix2x2T &p) {
		return Matrix2x2T(
			p.m_elem[0][0] * s, p.m_elem[0][1] * s,
			p.m_elem[1][0] * s, p.m_elem[1][1] * s);
	}

	//! *-Operator : Matrix * Vector
	Vector2T<Scalar> operator*(const Vector2T<Scalar> &vec) const {
		return Vector2T<Scalar>(
			m_elem[0][0]*vec[0] + m_elem[0][1]*vec[1],
			m_elem[1][0]*vec[0] + m_elem[1][1]*vec[1]);
	}

	//! *-Operator : Matrix * Matrix
	Matrix2x2T operator*(const Matrix2x2T &p) const {
		Matrix2x2T result;
		result.m_elem[0][0] = m_elem[0][0]*p.m_elem[0][0] + m_elem[0][1]*p.m_elem[1][0];
		result.m_elem[0][1] = m_elem[0][0]*p.m_elem[0][1] + m_elem[0][1]*p.m_elem[1][1];
		result.m_elem[1][0] = m_elem[1][0]*p.m_elem[0][0] + m_elem[1][1]*p.m_elem[1][0];
		result.m_elem[1][1] = m_elem[1][0]*p.m_elem[0][1] + m_elem[1][1]*p.m_elem[1][1];

		return result;
	}

	//! Sets all elements of the matrix to zero
	void setZero() {
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 2; j++)
				m_elem[i][j] = 0;
	}

	//! Sets the matrix to the 2x2 identity matrix
	void setIdentity() {
		for(int i = 0; i < 2; i++)
		{
			for(int j = 0; j < 2; j++)
			{
				m_elem[i][j] = 0;
			}

			m_elem[i][i] = 1;
		}
	}

	//! Returns the determinant of the 2x2 matrix
	Scalar det() const {
		return m_elem[0][0]*m_elem[1][1] - m_elem[0][1]*m_elem[1][0];
	}

	//! Returns the trace of the 2x2 matrix
	Scalar trace() const {
		return m_elem[0][0] + m_elem[1][1];
	}

	//! Returns the transposed matrix
	Matrix2x2T transposed() const {
		return Matrix2x2T(
			m_elem[0][0], m_elem[1][0],
			m_elem[0][1], m_elem[1][1]);
	}

	//! Returns the inverse matrix
	Matrix2x2T inverse() const {
		Scalar a = (*this)(0, 0);
		Scalar b = (*this)(0, 1);
		Scalar c = (*this)(1, 0);
		Scalar d = (*this)(1, 1);

		Scalar det = a*d - b*c;

		Matrix2x2T M;
		M(0, 0) = d;
		M(0, 1) = -b;
		M(1, 0) = -c;
		M(1, 1) = a;

		return M * (1.0f / det);
	}

private:
	//! Matrix elements
	Scalar m_elem[2][2];
};

// helper for printing a vector
template<class T>
std::ostream& operator<<( std::ostream& os, const  Matrix2x2T<T> &mat )
{
	char buf[256];
	//snprintf(buf,256,"{%f,%f,%f,%f}", (double)mat(0,0), (double)mat(1,0), (double)mat(0,1), (double)mat(1,1)   );
	sprintf(buf,"{%f,%f,%f,%f}", (double)mat(0,0), (double)mat(1,0), (double)mat(0,1), (double)mat(1,1)   );
	os << std::string(buf); 
	return os;
}

//=====================================================================================*/

/*! A symmetric sparse matrix, using dynamic data structures to allow fill-ins. Slow but flexible.
	NOTE: Only the lower-triagonal elements of the matrix are stored and only
	those entries can be accessed! */
template<class T>
class SparseSymmetricDynamicRowMatrixT
{
public:
	/*! Constructor. The matrix will have \c numRowsCols rows and \c numRowsCols columns. 
		Only the lower-triagonal elements of the symmetric matrix will be stored. */
	SparseSymmetricDynamicRowMatrixT(size_t numRowsCols) {
		m_numCols = numRowsCols;
		m_rowData.resize(numRowsCols);
	}

	/*! Constructor. The matrix will have 0 rows and 0 columns.
		Only the lower-triagonal elements of the matrix will be stored. */
	SparseSymmetricDynamicRowMatrixT() {
		m_numCols = 0;
	}

	//! Removes all elements and resizes the matrix to size (0, 0)
	void Clear() {
		m_numCols = 0;
		m_rowData.clear();
	}

	virtual ~SparseSymmetricDynamicRowMatrixT() {}

	//! Resizes the matrix to \c numRowsCols rows and \c numRowsCols columns and sets all elements to zero
	void ClearResize(size_t numRowsCols) {
		Clear();
		m_numCols = numRowsCols;
		m_rowData.resize(numRowsCols);
	}

	//! Computes: b = this * x
	void MultVector(const vector<T> &x, vector<T> &b) const {
		// Set b to zero:
		for(size_t i=0; i<b.size(); i++)
			b[i] = 0;

		size_t nrows = GetNumRows();

		for(size_t row=0; row<nrows; row++)
		{
			const map<size_t, T> &rowData = m_rowData[row];

			T rowSum = 0;

			for(typename map<size_t, T>::const_iterator iter = rowData.begin(); iter != rowData.end(); iter++)
			{
				size_t col = iter->first;
				assert(col <= row);

				T val = iter->second;

				rowSum += val * x[col];

				if(col < row)
					b[col] += val * x[row];
			}

			b[row] += rowSum;
		}
	}

	/*! Modifies the matrix and the vector \c b so the linear system 'this * x = b' will have the solution \c value at index \c idx. The row
		\c idx of this matrix must be occupied. */
	void FixSolution(std::vector<T> &b, size_t idx, T value) {
		size_t n = (size_t)b.size();
		assert(GetNumCols() == n);
		assert(GetNumRows() == n);
		assert(idx >= 0 && idx < n);

		map<size_t, T> &rowData = m_rowData[idx];

		// Need at least one entry in row 'idx':
		assert(rowData.size() > 0);

		for(typename map<size_t, T>::iterator iter = rowData.begin(); iter != rowData.end(); iter++)
		{
			size_t col = iter->first;
			assert(col <= idx);

			b[col] -= iter->second * value;

			if(col == idx)
				iter->second = 1;
			else
				iter->second = 0;
		}

		b[idx] = value;

		for(size_t i=idx+1; i<n; i++)
		{
			T oldValue = GetAt(i, idx);
			if(oldValue != 0)
			{
				b[i] -= oldValue * value;
				GetAt(i, idx) = 0;
			}
		}
	}

	/*! Returns the element at (\c row, \c col). Only elements in the lower-triagonal part
		of the matrix can be accessed! */
	const T &operator()(size_t row, size_t col) const {
		return GetAt(row, col);
	}

	/*! Returns the element at (\c row, \c col). Only elements in the lower-triagonal part
		of the matrix can be accessed! */
	T &operator()(size_t row, size_t col) {
		return GetAt(row, col);
	}

	/*! Returns the element at (\c row, \c col). Only elements in the lower-triagonal part
		of the matrix can be accessed! */
	const T &GetAt(size_t row, size_t col) const {
		assert(row >= 0 && row < GetNumRows());
		assert(col >= 0 && col < GetNumCols());

		assert(row >= col);

		const map<size_t, T> &rowData = m_rowData[row];

		typename map<size_t, T>::const_iterator iter = rowData.find(col);
		if(iter == rowData.end())
			return m_zero;

		return iter->second;
	}

	/*! Returns the element at (\c row, \c col). Only elements in the lower-triagonal part
		of the matrix can be accessed! */
	T &GetAt(size_t row, size_t col) {
		assert(row >= 0 && row < GetNumRows());
		assert(col >= 0 && col < GetNumCols());

		assert(row >= col);

		map<size_t, T> &rowData = m_rowData[row];

		typename map<size_t, T>::iterator iter = rowData.find(col);
		if(iter == rowData.end())
		{
			// Q_ASSERT(FALSE);	// Element does not exist!
			// return *(T*)NULL;

			// fill-in:
			rowData[col] = 0;
			iter = rowData.find(col);
			assert(iter != rowData.end());
		}

		return iter->second;
	}

	//! Returns the element at (\c row, \c col). Elements from the upper-triagonal part are returned as well.
	const T &GetAtFull(size_t row, size_t col) const {
		if(row < col)
		{
			// Swap row and col:
			size_t tmp = row;
			row = col;
			col = tmp;
		}

		return GetAt(row, col);
	}

	//! Returns the number of rows of the matrix
	size_t GetNumRows() const { return m_rowData.size(); }

	//! Returns the number of columns of the matrix
	size_t GetNumCols() const { return m_numCols; }

	//! Returns the number of non-zero entries. The values above the diagonal are not counted.
	size_t GetNumNonZero() const {
		size_t nnz = 0;
		for(size_t i=0; i<m_rowData.size(); i++)
		{
			nnz += m_rowData[i].size();
		}

		return nnz;
	}

private:
	size_t m_numCols;
	vector<map<size_t, T> > m_rowData;

	static T m_zero;
};

template<class T>
T SparseSymmetricDynamicRowMatrixT<T>::m_zero = 0;

//=====================================================================================*/

//! Diagonally preconditioned conjugate gradient solver
template<class T>
class SparseLinSolverPCGT
{
public:
	/*! Solves the linear system 'matA * x = b' for the unknown vector \c x. \c matA must be symmetric, positive-definite.
		\c residual defines the desired accuracy of the solution and \c maxIterations is the maximum number of
		iterations to perform or -1 for an infinite amount of iterations. */
	bool SolveLinearSystem(
			SparseSymmetricDynamicRowMatrixT<T> &matA, vector<T> &x, const vector<T> &b, 
			T residual, int maxIterations) {

		size_t n = matA.GetNumRows();

		if(x.size() != n || b.size() != n)
		{
			std::cerr<<"SolveLinearSystem x and b dont match! "<<x.size() <<" "<<b.size()<<" n="<<n<<" \n";
			assert(false);	// Dimension mismatch!
			return false;
		}

		vector<T> precond(n);
		vector<T> r(n);
		vector<T> d(n);
		vector<T> q(n);
		vector<T> s(n);

		// Prepare preconditioner:
		for(size_t i=0; i<n; i++)
			precond[i] = 1 / matA(i, i);

		// From: Jonathan Richard Shewchuk, "An Introduction to the Conjugate Gradient Method
		// Without the Agonizing Pain", Edition 1 1/4, August 4, 1994
		// http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf:

		// r = b - A*x
		matA.MultVector(x, r);
		for(size_t i=0; i<n; i++)
			r[i] = b[i] - r[i];

		// d = M^-1*r
		for(size_t i=0; i<n; i++)
			d[i] = precond[i] * r[i];

		// deltaNew = r^T*d
		T deltaNew = dotProd(r, d);

		// delta0 = deltaNew
		// Note: Using 1.0 here
		// to interpret residual as absolute error
		T delta0 = 1.0; //deltaNew;

		int iter = 0;
		while(maxIterations == -1 || iter < maxIterations)
		{
			if(deltaNew <= residual*residual*delta0)
				break;

			// q = A*d
			matA.MultVector(d, q);

			// alpha = deltaNew / (d^T*q)
			T alpha = deltaNew / dotProd(d, q);

			// x += alpha*d
			for(size_t i=0; i<n; i++)
				x[i] += alpha*d[i];

			// NOTE: if i%50 == 0: r = b - A*x

			// r -= alpha*q;
			for(size_t i=0; i<n; i++)
				r[i] -= alpha*q[i];

			// s = M^-1*r
			for(size_t i=0; i<n; i++)
				s[i] = precond[i] * r[i];

			// deltaOld = deltaNew
			T deltaOld = deltaNew;

			// deltaNew = r^T*s
			deltaNew = dotProd(r, s);

			// beta = deltaNew / deltaOld
			T beta = deltaNew / deltaOld;

			// d = s + beta*d:
			for(size_t i=0; i<n; i++)
				d[i] = s[i] + beta*d[i];

			iter++;
			std::cout<<"PCG, iter="<<iter<<", deltaNew="<<sqrt(deltaNew)<< " vs "<< (residual) <<"\n";
		}

		return true;
	}

private:
	//! Computes the dot product of two vectors
	static T dotProd(const vector<T> &a, const vector<T> &b) {
		T v = 0;

		assert(a.size() == b.size());

		for(size_t i=0; i<a.size(); i++)
			v += a[i] * b[i];

		return v;
	}
};

//=================================================================================================================
//===================== Vector3 =================================================================================
//=================================================================================================================

template<typename T>
class Vector3T
{
public:
	//----------------------------------------------------------- constructors
	Vector3T() { m_xyz[0] = (T)0; m_xyz[1] = (T)0; m_xyz[2] = (T)0;}
	explicit Vector3T(T v) { m_xyz[0] = v; m_xyz[1] = v; m_xyz[2] = v;}
	Vector3T(T x, T y, T z)  { m_xyz[0] = x; m_xyz[1] = y; m_xyz[2] = z; }
	explicit Vector3T(const T* xyz) { m_xyz[0] = xyz[0]; m_xyz[1] = xyz[1]; m_xyz[2] = xyz[2]; }
	Vector3T(const Vector3T<T> &other) { *this = other; }

	//----------------------------------------------------------- cast internal data type
	template <typename Scalar2>
	explicit Vector3T(const Vector3T<Scalar2> &m) {
		for(int r=0; r<m.dim(); r++)
		{
			(*this)[r] = static_cast<T> (m[r]);
		}
	}

	//----------------------------------------------------------- assignment
	Vector3T<T> & operator=( const Vector3T<T> &other ) { m_xyz[0]= other.m_xyz[0]; m_xyz[1] = other.m_xyz[1]; m_xyz[2] = other.m_xyz[2]; return *this; }

	//----------------------------------------------------------- element access
	// read-write
	T& operator[](int i) { assert(i>=0 && i<dim()); return m_xyz[i]; }
	T& x() { return m_xyz[0]; }
	T& y() { return m_xyz[1]; }
	T& z() { return m_xyz[2]; }

	// read-only
	const T& operator[](int i) const { assert(i>=0 && i<dim()); return m_xyz[i]; }
	const T& x() const { return m_xyz[0]; }
	const T& y() const { return m_xyz[1]; }
	const T& z() const { return m_xyz[2]; }

	// write only
	void setX(T x) { m_xyz[0] = x; }
	void setY(T y) { m_xyz[1] = y; }
	void setZ(T z) { m_xyz[2] = z; }

	void xyz(T* xyz) const { xyz[0] = m_xyz[0]; xyz[1] = m_xyz[1]; xyz[2] = m_xyz[2]; }
	void setXYZ(const T* xyz) { m_xyz[0] = xyz[0]; m_xyz[1] = xyz[1]; m_xyz[2] = xyz[2]; }

	//----------------------------------------------------------- logic operations
	bool isNullVector() const { return (m_xyz[0] == (T)0) && (m_xyz[1] == (T)0) && (m_xyz[2] == (T)0); }

	bool operator==(const Vector3T<T> &other) const {
		if(x()!=other.x()) return false;
		if(y()!=other.y()) return false;
		if(z()!=other.z()) return false;
		return true;
	}

	bool operator!=(const Vector3T<T> &other) const {
		return !(*this == other);
	}

	//! Lexicographical comparison
	bool operator<(const Vector3T<T> &other) const {
		for(int i=0; i<3; i++)
		{
			if((*this)[i] != other[i])
				return (*this)[i] < other[i];
		}

		return false;
	}

	//! Lexicographical comparison
	bool operator>(const Vector3T<T> &other) const {
		for(int i=0; i<3; i++)
		{
			if((*this)[i] != other[i])
				return (*this)[i] > other[i];
		}

		return false;
	}

	//----------------------------------------------------------- basic arithmetic operations
	//----------------------------------------------------------- scalar arithmetic operations
	Vector3T<T> operator*( T s ) const {
		return Vector3T<T>(m_xyz[0]*s, m_xyz[1]*s, m_xyz[2]*s);
	}

	Vector3T<T>& operator*=( T s ) {
		m_xyz[0] *= s;
		m_xyz[1] *= s;
		m_xyz[2] *= s;
		return *this;
	}

	Vector3T<T> operator/( T s ) const {
		return Vector3T<T>(m_xyz[0]/s, m_xyz[1]/s, m_xyz[2]/s);
	}

	Vector3T<T>& operator/=( T s ) {
		m_xyz[0] /= s;
		m_xyz[1] /= s;
		m_xyz[2] /= s;
		return *this;
	}

	//----------------------------------------------------------- vector arithmetic operations
	Vector3T<T> operator+( const Vector3T<T> &v2 ) const {
		return Vector3T<T>(x()+v2.x(),y()+v2.y(),z()+v2.z());
	}

	Vector3T<T>& operator+=( const Vector3T<T> &v ) {
		m_xyz[0] += v.x();
		m_xyz[1] += v.y();
		m_xyz[2] += v.z();
		return *this;
	}

	Vector3T<T> operator-( const Vector3T<T> &v2 ) const {
		return Vector3T<T>(x()-v2.x(),y()-v2.y(),z()-v2.z());
	}

	Vector3T<T>& operator-=( const Vector3T<T> &v ) {
		m_xyz[0] -= v.x();
		m_xyz[1] -= v.y();
		m_xyz[2] -= v.z();
		return *this;
	}

	Vector3T<T> operator-() const {
		return Vector3T<T>(-m_xyz[0], -m_xyz[1], -m_xyz[2]);
	}

	Vector3T<T> componentWiseProduct(const Vector3T<T> &v2) const {
		return Vector3T<T>(x()*v2.x(), y()*v2.y(), z()*v2.z());
	}

	Vector3T<T> componentWiseDivision(const Vector3T<T> &v2) const {
		return Vector3T<T>(x()/v2.x(), y()/v2.y(), z()/v2.z());
	}

	//! Inner product (dot product)
	T operator*(  const Vector3T<T> &v2 ) const {
		return x()*v2.x()+y()*v2.y()+z()*v2.z();
	}

	//! Inner product (dot product)
	T operator|(const Vector3T<T> &v2) const {
		return (*this) * v2;
	}

	//! Cross product
	Vector3T<T> cross( const Vector3T<T> &v2 ) const {
		return Vector3T<T>(y()*v2.z()-v2.y()*z(), z()*v2.x()-v2.z()*x(), x()*v2.y()-v2.x()*y() );
	}

	//! Cross product
	Vector3T<T> operator%(const Vector3T<T> &v2) const {
		return cross(v2);
	}

	//----------------------------------------------------------- other
	void normalize() {
		*this /= length();
	}

	Vector3T<T> normalized() const {
		return *this / length();
	}

	T length() const {
		return sqrt(squaredLength());
	}

	T norm() const {
		return length();
	}

	T squaredLength() const
	{
		return m_xyz[0]*m_xyz[0]+m_xyz[1]*m_xyz[1]+m_xyz[2]*m_xyz[2];
	}

	T sqrnorm() const {
		return squaredLength();
	}

	//! Returns an arbitrary vector which is orthogonal to this vector
	Vector3T<T> orthogonalVector() const {
		Vector3T<T> vecNorm = normalized();

		const int dim = 3;

		int minIdx = -1;
		T minDot = 0;

		for(int i=0; i<dim; i++)
		{
			// Create i-th basis vector:
			Vector3T<T> basisVec;
			for(int j=0; j<dim; j++)
				basisVec[j] = (j==i);

			T dotProd = vecNorm | basisVec;

			if(minIdx == -1 || fabs(dotProd) < fabs(minDot))
			{
				minIdx = i;
				minDot = dotProd;
			}
		}

		// Create minIdx-th basis vector:
		Vector3T<T> basisVec;
		for(int j=0; j<dim; j++)
			basisVec[j] = (j==minIdx);

		return basisVec - vecNorm * minDot;
	}

	//! Returns the minimal component
	T minComponent() const {
		T m = (*this)[0];

		for(int i=1; i<3; i++)
		{
			if((*this)[i] < m)
				m = (*this)[i];
		}

		return m;
	}

	//! Returns the maximal component
	T maxComponent() const {
		T m = (*this)[0];

		for(int i=1; i<3; i++)
		{
			if((*this)[i] > m)
				m = (*this)[i];
		}

		return m;
	}

	//! Component-wise minimum
	Vector3T<T> componentWiseMin(const Vector3T<T> &v2) const {
		Vector3T<T> v(*this);
		v.minimize(v2);
		return v;
	}

	//! Component-wise maximum
	Vector3T<T> componentWiseMax(const Vector3T<T> &v2) const {
		Vector3T<T> v(*this);
		v.maximize(v2);
		return v;
	}

	//! Minimizes values: same as *this = min(*this, v2)
	void minimize(const Vector3T<T> &v2) {
		for(int i=0; i<3; i++)
		{
			if(v2[i] < (*this)[i])
				(*this)[i] = v2[i];
		}
	}

	//! Maximizes values: same as *this = max(*this, v2)
	void maximize(const Vector3T<T> &v2) {
		for(int i=0; i<3; i++)
		{
			if(v2[i] > (*this)[i])
				(*this)[i] = v2[i];
		}
	}

	//! Compute the reflection of this vector with \c n
	Vector3T<T> reflectionAt( const Vector3T<T> &n ) const {
		return *this-2*(n*(*this))*n;
	}
	
    //! Returns a pointer to the raw vector data of length 3
	const T *data() const { return m_xyz; }

	//! Cast to (const T *). See \c data().
	operator const T *() const { return data(); }

//== STATIC MEMBERS =========================================================
public:
	static int dim() { return 3; }

	static const Vector3T<T> &zero() {
		static const Vector3T<T> zeroVec(0,0,0);
		return zeroVec;
	}

private:
	T m_xyz[3];
};

//== NON MEMBER DEFINITION =========================================================
// scalar * vector
template<class T>
Vector3T<T> operator*(T const &s, Vector3T<T> const &v)
{
	return v*s;
}

//== TYPEDEFS =========================================================
typedef Vector3T<float> Vector3f;
typedef Vector3T<double> Vector3d;




//=================================================================================================================
//===================== Matrix3x3 =================================================================================
//=================================================================================================================

template<typename Scalar>
class Matrix3x3T
{
public:
	//! Empty constructor - object remains uninitialized.
	Matrix3x3T() {}

	//! Copy-Constructor
	Matrix3x3T(const Matrix3x3T &p) {
		m_elem[0][0] = p.m_elem[0][0];
		m_elem[0][1] = p.m_elem[0][1];
		m_elem[0][2] = p.m_elem[0][2];
		m_elem[1][0] = p.m_elem[1][0];
		m_elem[1][1] = p.m_elem[1][1];
		m_elem[1][2] = p.m_elem[1][2];
		m_elem[2][0] = p.m_elem[2][0];
		m_elem[2][1] = p.m_elem[2][1];
		m_elem[2][2] = p.m_elem[2][2];
	}

	//! constructor using individual matrix element
	Matrix3x3T(Scalar x00, Scalar x01, Scalar x02, Scalar x10, Scalar x11, Scalar x12, Scalar x20, Scalar x21, Scalar x22){
		m_elem[0][0] = x00;
		m_elem[0][1] = x01;
		m_elem[0][2] = x02;
		m_elem[1][0] = x10;
		m_elem[1][1] = x11;
		m_elem[1][2] = x12;
		m_elem[2][0] = x20;
		m_elem[2][1] = x21;
		m_elem[2][2] = x22;
	}

	Scalar& operator() (int _i, int _j) {
		return m_elem[_i][_j];
	}

	const Scalar& operator() (int _i, int _j) const {
		return m_elem[_i][_j];
	}

	//! +=-Operator
	Matrix3x3T& operator+=( const Matrix3x3T& p ) {
		m_elem[0][0] += p.m_elem[0][0];
		m_elem[0][1] += p.m_elem[0][1];
		m_elem[0][2] += p.m_elem[0][2];
		m_elem[1][0] += p.m_elem[1][0];
		m_elem[1][1] += p.m_elem[1][1];
		m_elem[1][2] += p.m_elem[1][2];
		m_elem[2][0] += p.m_elem[2][0];
		m_elem[2][1] += p.m_elem[2][1];
		m_elem[2][2] += p.m_elem[2][2];
		return *this;
	}

	//! -=-Operator
	Matrix3x3T& operator-=( const Matrix3x3T& p ) {
		m_elem[0][0] -= p.m_elem[0][0];
		m_elem[0][1] -= p.m_elem[0][1];
		m_elem[0][2] -= p.m_elem[0][2];
		m_elem[1][0] -= p.m_elem[1][0];
		m_elem[1][1] -= p.m_elem[1][1];
		m_elem[1][2] -= p.m_elem[1][2];
		m_elem[2][0] -= p.m_elem[2][0];
		m_elem[2][1] -= p.m_elem[2][1];
		m_elem[2][2] -= p.m_elem[2][2];
		return *this;
	}

	//! /=-Operator
	Matrix3x3T& operator*=( Scalar s ) {
		m_elem[0][0] *= s;
		m_elem[0][1] *= s;
		m_elem[0][2] *= s;
		m_elem[1][0] *= s;
		m_elem[1][1] *= s;
		m_elem[1][2] *= s;
		m_elem[2][0] *= s;
		m_elem[2][1] *= s;
		m_elem[2][2] *= s;
		return *this;
	}

	//! *=-Operator : Matrix * Matrix
	Matrix3x3T& operator*=( const Matrix3x3T& p ) {
		return ( *this = *this * p );
	}

	//! +-Operator
	Matrix3x3T  operator+( const Matrix3x3T& p ) const {
		return Matrix3x3T( m_elem[0][0] + p.m_elem[0][0],
			m_elem[0][1] + p.m_elem[0][1],
			m_elem[0][2] + p.m_elem[0][2],
			m_elem[1][0] + p.m_elem[1][0],
			m_elem[1][1] + p.m_elem[1][1],
			m_elem[1][2] + p.m_elem[1][2],
			m_elem[2][0] + p.m_elem[2][0],
			m_elem[2][1] + p.m_elem[2][1],
			m_elem[2][2] + p.m_elem[2][2] );
	}

	//! --Operator
	Matrix3x3T  operator-( const Matrix3x3T& p ) const {
		return Matrix3x3T( m_elem[0][0] - p.m_elem[0][0],
			m_elem[0][1] - p.m_elem[0][1],
			m_elem[0][2] - p.m_elem[0][2],
			m_elem[1][0] - p.m_elem[1][0],
			m_elem[1][1] - p.m_elem[1][1],
			m_elem[1][2] - p.m_elem[1][2],
			m_elem[2][0] - p.m_elem[2][0],
			m_elem[2][1] - p.m_elem[2][1],
			m_elem[2][2] - p.m_elem[2][2] );
	}

	//! -Operator
	Matrix3x3T  operator/( Scalar s ) const {
		return Matrix3x3T( m_elem[0][0] / s, m_elem[0][1] / s, m_elem[0][2] / s,
			m_elem[1][0] / s, m_elem[1][1] / s, m_elem[1][2] / s,
			m_elem[2][0] / s, m_elem[2][1] / s, m_elem[2][2] / s );
	}

	//! *-Operator : Matrix * Scalar
	Matrix3x3T  operator*( Scalar s ) const {
		return Matrix3x3T( m_elem[0][0] * s, m_elem[0][1] * s, m_elem[0][2] * s,
			m_elem[1][0] * s, m_elem[1][1] * s, m_elem[1][2] * s,
			m_elem[2][0] * s, m_elem[2][1] * s, m_elem[2][2] * s );
	}

	// friend *-Operator : Scalar * Matrix
	friend Matrix3x3T operator*( Scalar s, const Matrix3x3T& p) {
		return Matrix3x3T( p.m_elem[0][0] * s, p.m_elem[0][1] * s, p.m_elem[0][2] * s,
			p.m_elem[1][0] * s, p.m_elem[1][1] * s, p.m_elem[1][2] * s,
			p.m_elem[2][0] * s, p.m_elem[2][1] * s, p.m_elem[2][2] * s );
	}

	//! *-Operator : Matrix * Vector
	Vector3T<Scalar> operator*( const Vector3T<Scalar>& vec ) const {
		return Vector3T<Scalar>( m_elem[0][0]*vec[0] +
			m_elem[0][1]*vec[1] +
			m_elem[0][2]*vec[2],
			m_elem[1][0]*vec[0] +
			m_elem[1][1]*vec[1] +
			m_elem[1][2]*vec[2],
			m_elem[2][0]*vec[0] +
			m_elem[2][1]*vec[1] +
			m_elem[2][2]*vec[2] );
	}

	// Computes ((*this)^T * vec)
	Matrix3x3T  operator*( const Matrix3x3T& p ) const {
		Matrix3x3T result;
		result.m_elem[0][0] = (m_elem[0][0]*p.m_elem[0][0] +
			m_elem[0][1]*p.m_elem[1][0] +
			m_elem[0][2]*p.m_elem[2][0] );
		result.m_elem[0][1] = ( m_elem[0][0]*p.m_elem[0][1] +
			m_elem[0][1]*p.m_elem[1][1] +
			m_elem[0][2]*p.m_elem[2][1] );
		result.m_elem[0][2] = ( m_elem[0][0]*p.m_elem[0][2] +
			m_elem[0][1]*p.m_elem[1][2] +
			m_elem[0][2]*p.m_elem[2][2] );
		result.m_elem[1][0] = ( m_elem[1][0]*p.m_elem[0][0] +
			m_elem[1][1]*p.m_elem[1][0] +
			m_elem[1][2]*p.m_elem[2][0] );
		result.m_elem[1][1] = ( m_elem[1][0]*p.m_elem[0][1] +
			m_elem[1][1]*p.m_elem[1][1] +
			m_elem[1][2]*p.m_elem[2][1] );
		result.m_elem[1][2] = ( m_elem[1][0]*p.m_elem[0][2] +
			m_elem[1][1]*p.m_elem[1][2] +
			m_elem[1][2]*p.m_elem[2][2] );
		result.m_elem[2][0] = ( m_elem[2][0]*p.m_elem[0][0] +
			m_elem[2][1]*p.m_elem[1][0] +
			m_elem[2][2]*p.m_elem[2][0] );
		result.m_elem[2][1] = ( m_elem[2][0]*p.m_elem[0][1] +
			m_elem[2][1]*p.m_elem[1][1] +
			m_elem[2][2]*p.m_elem[2][1] );
		result.m_elem[2][2] = ( m_elem[2][0]*p.m_elem[0][2] +
			m_elem[2][1]*p.m_elem[1][2] +
			m_elem[2][2]*p.m_elem[2][2] );
		return result;
	}

	// Computes ((*this)^T * p)
	//! transform vector (x',y',0) = A * (x,y,0)
	void zero() {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				m_elem[i][j] = (Scalar)0;
	}


	//! build identity matrix
	void identity() {
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
				m_elem[i][j] = (Scalar)0;
			m_elem[i][i] = (Scalar)1;
		}
	}

	//! returns true if all matrix m_elements are exactly equal to zero
	//! determinant of 3x3 Matrix
	Scalar det() const
	{
		return ( (m_elem[0][1]*m_elem[1][2] - m_elem[0][2]*m_elem[1][1]) * m_elem[2][0]
		+ (m_elem[0][2]*m_elem[1][0] - m_elem[0][0]*m_elem[1][2]) * m_elem[2][1]
		+ (m_elem[0][0]*m_elem[1][1] - m_elem[0][1]*m_elem[1][0]) * m_elem[2][2] );
	}

	//! trace of a 3x3 Matrix
	Scalar trace() const {
		return m_elem[0][0] + m_elem[1][1] + m_elem[2][2];
	}

	//! Transposed Matrix
	Matrix3x3T transposed() const {
		return( Matrix3x3T( m_elem[0][0], m_elem[1][0], m_elem[2][0],
			m_elem[0][1], m_elem[1][1], m_elem[2][1],
			m_elem[0][2], m_elem[1][2], m_elem[2][2] ) );
	}

	void transpose() {
		Scalar a = m_elem[1][0];
		m_elem[1][0] = m_elem[0][1];
		m_elem[0][1] = a;

		a = m_elem[2][0];
		m_elem[2][0] = m_elem[0][2];
		m_elem[0][2] = a;

		a = m_elem[2][1];
		m_elem[2][1] = m_elem[1][2];
		m_elem[1][2] = a;
	}

	Matrix3x3T<Scalar> inverse() const {
		Matrix3x3T<Scalar> M;

		computeInverseUnsave(M);

		return M;
	}

	//! Returns the inverse of the matrix, or sets \c ok to \c false if the matrix is (nearly) singular
	Matrix3x3T<Scalar> inverse(bool &ok) const {
		Matrix3x3T<Scalar> M;

		ok = computeInverse(M, std::numeric_limits<Scalar>::epsilon());

		return M;
	}

	//! Inverts the matrix or returns \c false if the matrix is (nearly) singular
	bool invert() {
		bool ok = false;
		*this = inverse(ok);
		return ok;
	}

	const Scalar *data() const { return &m_elem[0][0]; }

	//! Cast to (const Scalar *). See \c data().
	operator const Scalar *() const { return data(); }

private:
	bool computeInverse(Matrix3x3T<Scalar> &M, Scalar eps) const {
		assert(&M != this);

		Scalar d =  - m_elem[0][0]*m_elem[1][1]*m_elem[2][2]
		+ m_elem[0][0]*m_elem[1][2]*m_elem[2][1]
		+ m_elem[1][0]*m_elem[0][1]*m_elem[2][2]
		- m_elem[1][0]*m_elem[0][2]*m_elem[2][1]
		- m_elem[2][0]*m_elem[0][1]*m_elem[1][2]
		+ m_elem[2][0]*m_elem[0][2]*m_elem[1][1];

		if( fabs(d) <= eps )
		{
			// Need to make sure something gets written to \c M:
			M.setValues(std::numeric_limits<Scalar>::quiet_NaN());

			return false;
		}

		// d = 1/d;
		M(0,0) = (m_elem[1][2]*m_elem[2][1] - m_elem[1][1]*m_elem[2][2]) / d;
		M(0,1) = (m_elem[0][1]*m_elem[2][2] - m_elem[0][2]*m_elem[2][1]) / d;
		M(0,2) = (m_elem[0][2]*m_elem[1][1] - m_elem[0][1]*m_elem[1][2]) / d;
		M(1,0) = (m_elem[1][0]*m_elem[2][2] - m_elem[1][2]*m_elem[2][0]) / d;
		M(1,1) = (m_elem[0][2]*m_elem[2][0] - m_elem[0][0]*m_elem[2][2]) / d;
		M(1,2) = (m_elem[0][0]*m_elem[1][2] - m_elem[0][2]*m_elem[1][0]) / d;
		M(2,0) = (m_elem[1][1]*m_elem[2][0] - m_elem[1][0]*m_elem[2][1]) / d;
		M(2,1) = (m_elem[0][0]*m_elem[2][1] - m_elem[0][1]*m_elem[2][0]) / d;
		M(2,2) = (m_elem[0][1]*m_elem[1][0] - m_elem[0][0]*m_elem[1][1]) / d;

		return true;
	}

	void computeInverseUnsave(Matrix3x3T<Scalar> &M) const {
		assert(&M != this);

		Scalar d =  - m_elem[0][0]*m_elem[1][1]*m_elem[2][2]
		+ m_elem[0][0]*m_elem[1][2]*m_elem[2][1]
		+ m_elem[1][0]*m_elem[0][1]*m_elem[2][2]
		- m_elem[1][0]*m_elem[0][2]*m_elem[2][1]
		- m_elem[2][0]*m_elem[0][1]*m_elem[1][2]
		+ m_elem[2][0]*m_elem[0][2]*m_elem[1][1];

		// d = 1/d;
		M(0,0) = (m_elem[1][2]*m_elem[2][1] - m_elem[1][1]*m_elem[2][2]) / d;
		M(0,1) = (m_elem[0][1]*m_elem[2][2] - m_elem[0][2]*m_elem[2][1]) / d;
		M(0,2) = (m_elem[0][2]*m_elem[1][1] - m_elem[0][1]*m_elem[1][2]) / d;
		M(1,0) = (m_elem[1][0]*m_elem[2][2] - m_elem[1][2]*m_elem[2][0]) / d;
		M(1,1) = (m_elem[0][2]*m_elem[2][0] - m_elem[0][0]*m_elem[2][2]) / d;
		M(1,2) = (m_elem[0][0]*m_elem[1][2] - m_elem[0][2]*m_elem[1][0]) / d;
		M(2,0) = (m_elem[1][1]*m_elem[2][0] - m_elem[1][0]*m_elem[2][1]) / d;
		M(2,1) = (m_elem[0][0]*m_elem[2][1] - m_elem[0][1]*m_elem[2][0]) / d;
		M(2,2) = (m_elem[0][1]*m_elem[1][0] - m_elem[0][0]*m_elem[1][1]) / d;
	}

private:
	// matrix elements
	Scalar m_elem[3][3];
};

typedef Matrix3x3T<float>  Matrix3x3f;
typedef Matrix3x3T<double> Matrix3x3d;


#endif
