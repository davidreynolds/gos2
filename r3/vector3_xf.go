package r3

import (
	"github.com/davidreynolds/gos2/exactfloat"
)

type Vector3_xf struct {
	X, Y, Z exactfloat.ExactFloat
}

func Vector3_xf_FromVector(v Vector) Vector3_xf {
	return Vector3_xf{
		X: exactfloat.NewExactFloat(v.X),
		Y: exactfloat.NewExactFloat(v.Y),
		Z: exactfloat.NewExactFloat(v.Z),
	}
}

func (a Vector3_xf) CrossProd(b Vector3_xf) Vector3_xf {
	return Vector3_xf{
		a.Y.Mul(b.Z).Sub(a.Z.Mul(b.Y)),
		a.Z.Mul(b.X).Sub(a.X.Mul(b.Z)),
		a.X.Mul(b.Y).Sub(a.Y.Mul(b.X)),
	}
}

func (a Vector3_xf) DotProd(b Vector3_xf) exactfloat.ExactFloat {
	x := a.X.Mul(b.X)
	y := a.Y.Mul(b.Y)
	z := a.Z.Mul(b.Z)
	return x.Add(y).Add(z)
}

func (a Vector3_xf) Mul(m exactfloat.ExactFloat) Vector3_xf {
	return Vector3_xf{a.X.Mul(m), a.Y.Mul(m), a.Z.Mul(m)}
}

func (a Vector3_xf) ApproxEqual(b Vector3_xf) bool {
	epsilon := exactfloat.NewExactFloat(1e-14)
	return exactfloat.Abs(a.X.Sub(b.X)).LessThan(epsilon) &&
		exactfloat.Abs(a.Y.Sub(b.Y)).LessThan(epsilon) &&
		exactfloat.Abs(a.Z.Sub(b.Z)).LessThan(epsilon)
}
