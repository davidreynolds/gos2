package r3

/*
type Matrix struct {
	m [3][3]float64
}
*/

type Matrix [3][3]float64

func MatrixFromCols(v0, v1, v2 Vector) Matrix {
	m := Matrix{}
	m.Set(
		v0.X, v1.X, v2.X,
		v0.Y, v1.Y, v2.Y,
		v0.Z, v1.Z, v2.Z,
	)
	return m
}

func Matrix9(
	m00, m01, m02,
	m10, m11, m12,
	m20, m21, m22 float64) Matrix {
	m := Matrix{}
	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02

	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12

	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	return m
}

func (m Matrix) Transpose() Matrix {
	return Matrix9(
		m[0][0], m[1][0], m[2][0],
		m[0][1], m[1][1], m[2][1],
		m[0][2], m[1][2], m[2][2],
	)
}

func (m *Matrix) Set(
	m00, m01, m02,
	m10, m11, m12,
	m20, m21, m22 float64) Matrix {
	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02

	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12

	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	return *m
}

func (m Matrix) Col(i int) Vector {
	return Vector{m[0][i], m[1][i], m[2][i]}
}

func (m *Matrix) SetCol(i int, v Vector) {
	m[0][i] = v.X
	m[1][i] = v.Y
	m[2][i] = v.Z
}

func (m Matrix) MulVector(v Vector) Vector {
	return Vector{
		m[0][0]*v.X + m[0][1]*v.Y + m[0][2]*v.Z,
		m[1][0]*v.X + m[1][1]*v.Y + m[1][2]*v.Z,
		m[2][0]*v.X + m[2][1]*v.Y + m[2][2]*v.Z,
	}
}
