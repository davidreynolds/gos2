package s2

import (
	"math"

	"github.com/davidreynolds/gos2/r3"
	"github.com/davidreynolds/gos2/s1"
)

// Point represents a point on the unit sphere as a normalized 3D vector.
//
// Points are guaranteed to be close to normal in the sense that the norm of any points will be very close to 1.
//
// Fields should be treated as read-only. Use one of the factory methods for creation.
type Point struct {
	r3.Vector
}

const (
	maxDetError = 0.8e-15 // 14 * (2**-54)
)

// PointFromCoords creates a new normalized point from coordinates.
//
// This always returns a valid point. If the given coordinates can not be normalized the origin point will be returned.
func PointFromCoords(x, y, z float64) Point {
	if x == 0 && y == 0 && z == 0 {
		return OriginPoint()
	}
	return Point{r3.Vector{x, y, z}.Normalize()}
}

// OriginPoint returns a unique "origin" on the sphere for operations that need a fixed
// reference point. In particular, this is the "point at infinity" used for
// point-in-polygon testing (by counting the number of edge crossings).
//
// It should *not* be a point that is commonly used in edge tests in order
// to avoid triggering code to handle degenerate cases (this rules out the
// north and south poles). It should also not be on the boundary of any
// low-level S2Cell for the same reason.
func OriginPoint() Point {
	return Point{r3.Vector{0.00456762077230, 0.99947476613078, 0.03208315302933}}
}

// PointCross returns a Point that is orthogonal to both p and op. This is similar to
// p.Cross(op) (the true cross product) except that it does a better job of
// ensuring orthogonality when the Point is nearly parallel to op, it returns
// a non-zero result even when p == op or p == -op and the result is a Point,
// so it will have norm 1.
//
// It satisfies the following properties (f == PointCross):
//
//   (1) f(p, op) != 0 for all p, op
//   (2) f(op,p) == -f(p,op) unless p == op or p == -op
//   (3) f(-p,op) == -f(p,op) unless p == op or p == -op
//   (4) f(p,-op) == -f(p,op) unless p == op or p == -op
func (p Point) PointCross(op Point) Point {
	// NOTE(dnadasi): In the C++ API the equivalent method here was known as "RobustCrossProd",
	// but PointCross more accurately describes how this method is used.
	x := p.Add(op.Vector).Cross(op.Sub(p.Vector))

	if x.ApproxEqual(r3.Vector{0, 0, 0}) {
		// The only result that makes sense mathematically is to return zero, but
		// we find it more convenient to return an arbitrary orthogonal vector.
		return Point{p.Ortho()}
	}

	return Point{x.Normalize()}
}

// CCW returns true if the points A, B, C are strictly counterclockwise,
// and returns false if the points are clockwise or collinear (i.e. if they are all
// contained on some great circle).
//
// Due to numerical errors, situations may arise that are mathematically
// impossible, e.g. ABC may be considered strictly CCW while BCA is not.
// However, the implementation guarantees the following:
//
//   If CCW(a,b,c), then !CCW(c,b,a) for all a,b,c.
func CCW(a, b, c Point) bool {
	// NOTE(dnadasi): In the C++ API the equivalent method here was known as "SimpleCCW",
	// but CCW seems like a fine name at least until the need for a RobustCCW is demonstrated.

	// We compute the signed volume of the parallelepiped ABC. The usual
	// formula for this is (A ⨯ B) · C, but we compute it here using (C ⨯ A) · B
	// in order to ensure that ABC and CBA are not both CCW. This follows
	// from the following identities (which are true numerically, not just
	// mathematically):
	//
	//     (1) x ⨯ y == -(y ⨯ x)
	//     (2) -x · y == -(x · y)
	return c.Cross(a.Vector).Dot(b.Vector) > 0
}

// The following function returns the sign of the determinant of three points
// A, B, C under a model where every possible Point is slightly perturbed by
// a unique infinitesmal amount such that no three perturbed points are
// colinear and no four points are coplanar. The perturbations are so small
// that they do not change the sign of any determinant that was non-zero
// before the perturbations, and therefore can be safely ignored unless the
// determinant of three points is exactly zero (using multiple-precision
// arithmetic).
//
// Since the symbolic perturbation of a given point is fixed (i.e., the
// perturbation is the same for all calls to this method and does not depend
// on the other two arguments), the results of this method are always
// self-consistent. It will never return results that would correspond to an
// "impossible" configuration of non-degenerate points.
//
// Requirements:
//   The 3x3 determinant of A, B, C must be exactly zero.
//   The points must be distinct, with A < B < C in lexicographic order.
//
// Returns:
//   +1 or -1 according to the sign of the determinant after the symbolic
// perturbations are taken into account.
//
// Reference:
//   "Simulation of Simplicity" (Edelsbrunner and Muecke, ACM Transactions on
//   Graphics, 1990).
//
func SymbolicallyPerturbedCCW(a, b, c, b_cross_c r3.Vector3_xf) int {
	// This method requires that the points are sorted in lexicographically
	// increasing order. This is because every possible Point has its own
	// symbolic perturbation such that if A < B then the symbolic
	// perturbation for A is much larger than the perturbation for B.
	//
	// Alternatively, we could sort the points in this method and keep track
	// of the sign of the permutation, but it is more efficient to do this
	// before converting the inputs to the multi-precision representation,
	// and this also lets us re-use the result of the cross product B x C.

	// Every input coordinate x[i] is assigned a symbolic perturbation
	// dx[i]. We then compute the sign of the determinant of the perturbed
	// points, i.e.
	//               | a[0]+da[0]  a[1]+da[1]  a[2]+da[2] |
	//               | b[0]+db[0]  b[1]+db[1]  b[2]+db[2] |
	//               | c[0]+dc[0]  c[1]+dc[1]  c[2]+dc[2] |
	//
	// The perturbations are chosen such that
	//
	//   da[2] > da[1] > da[0] > db[2] > db[1] > db[0] > dc[2] > dc[1] > dc[0]
	//
	// where each perturbation is so much smaller than the previous one that
	// we don't even need to consider it unless the coefficients of all
	// previous perturbations are zero. In fact, it is so small that we
	// don't need to consider it unless the coefficients of all products of
	// the previous perturbations are zero. For example, we don't need to
	// consider the coefficient of db[1] unless the coefficient of
	// db[2]*da[0] is zero.
	//
	// The following code simply enumerates the coefficients of the
	// perturbations (and products of perturbations) that appear in the
	// determinant above, in order of decreasing perturbation magnitude.
	// The first non-zero coefficient determines the sign of the result.
	// The easiest way to enumerate the coefficients in the correct order
	// is to pretend that each perturbation is some tiny value "eps" raised
	// to a power of two:
	//
	// eps**    1      2      4      8     16     32     64     128    256
	//        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]  dc[2]  dc[1]  dc[0]
	//
	// Essentially we can then just count in binary and test the
	// corresponding subset of perturbations at each step. So for example,
	// we must test the coefficient of db[2]*da[0] before db[1] because
	// eps**12 > eps16.
	//
	// Of course, not all products of these perturbations appear in the
	// determinant above, since the determinant only contains the products
	// of elements in distinct rows and columns. Thus we don't need to
	// consider da[2]*da[1], db[1]*da[1], etc. Furthermore, sometimes
	// different pairs of perturbations have the same coefficient in the
	// determinant; for example, da[1]*db[0] and db[1]*da[0] have the
	// same coefficient (c[2]). Therefore we only need to test this
	// coefficient the first time we encounter it in the binary order above
	// (which will be db[1]*da[0]).
	//
	// The sequence of tests below also appears in Table 4-ii of the paper
	// referenced above, if you just want to look it up, with the following
	// translations: [a,b,c] -> [i,j,k] and [0,1,2] -> [1,2,3]. Also note
	// that some of the signs are different because the opposite cross
	// product is used (e.g., B x C rather than C x B).

	if det_sign := b_cross_c.Z.Sgn(); det_sign != 0 { // da[2]
		return det_sign
	}
	if det_sign := b_cross_c.Y.Sgn(); det_sign != 0 { // da[1]
		return det_sign
	}
	if det_sign := b_cross_c.X.Sgn(); det_sign != 0 { // da[0]
		return det_sign
	}

	first := c.X.Mul(a.Y)
	second := c.Y.Mul(a.X)

	if det_sign := first.Sub(second).Sgn(); det_sign != 0 { // db[2]
		return det_sign
	}
	if det_sign := c.X.Sgn(); det_sign != 0 { // db[2] * da[1]
		return det_sign
	}
	if det_sign := -(c.Y.Sgn()); det_sign != 0 { // db[2] * da[0]
		return det_sign
	}

	first = c.Z.Mul(a.X)
	second = c.X.Mul(a.Z)

	if det_sign := first.Sub(second).Sgn(); det_sign != 0 { // db[1]
		return det_sign
	}
	if det_sign := c.Z.Sgn(); det_sign != 0 { // db[1] * da[0]
		return det_sign
	}

	// The previous tests guarantee that C == (0, 0, 0).

	first = a.X.Mul(b.Y)
	second = a.Y.Mul(b.X)

	if det_sign := first.Sub(second).Sgn(); det_sign != 0 { // dc[2]
		return det_sign
	}
	if det_sign := -(b.X.Sgn()); det_sign != 0 { // dc[2] * da[1]
		return det_sign
	}
	if det_sign := b.Y.Sgn(); det_sign != 0 { // dc[2] * da[0]
		return det_sign
	}
	if det_sign := a.X.Sgn(); det_sign != 0 { // dc[2] * db[1]
		return det_sign
	}
	return 1 // dc[2] * db[1] * da[0]
}

func ExpensiveCCW(a, b, c Point) int {
	// Return zero if and only if two points are the same. This ensures (1).
	if a == b || b == c || c == a {
		return 0
	}

	// Sort the three points in lexicographic order, keeping track of the
	// sign of the permutation. (Each exchange inverts the sign of the
	// determinant).
	perm_sign := 1
	pa := a
	pb := b
	pc := c
	if pa.GreaterThan(pb.Vector) {
		pa, pb = pb, pa
		perm_sign = -perm_sign
	}
	if pb.GreaterThan(pc.Vector) {
		pb, pc = pc, pb
		perm_sign = -perm_sign
	}
	if pa.GreaterThan(pb.Vector) {
		pa, pb = pb, pa
		perm_sign = -perm_sign
	}

	// Construct multiple-precision versions of the sorted points and
	// compute their exact 3x3 determinant.
	xa := r3.Vector3_xf_FromVector(pa.Vector)
	xb := r3.Vector3_xf_FromVector(pb.Vector)
	xc := r3.Vector3_xf_FromVector(pc.Vector)
	xb_cross_xc := xb.CrossProd(xc)
	det := xa.DotProd(xb_cross_xc)

	// The precision of ExactFloat is high enough that the result should
	// always be exact (no rounding was performed).

	// If the exact determinant is non-zero, we're done.
	det_sign := det.Sgn()
	if det_sign == 0 {
		// Otherwise, we need to resort to symbolic perturbations to
		// resolve the sign of the determinant.
		det_sign = SymbolicallyPerturbedCCW(xa, xb, xc, xb_cross_xc)
	}
	return perm_sign * det_sign
}

func RobustCCW(a, b, c Point) int {
	// We don't need PointCross (RobustCrossProd) here because RobustCCW
	// does its own error estimation and calls ExpensiveCCW if there is
	// any uncertainty about the result.
	return RobustCCW2(a, b, c, Point{a.Cross(b.Vector)})
}

func RobustCCW2(a, b, c, a_cross_b Point) int {
	ccw := TriageCCW(a, b, c, a_cross_b)
	if ccw == 0 {
		ccw = ExpensiveCCW(a, b, c)
	}
	return ccw
}

func TriageCCW(a, b, c, a_cross_b Point) int {
	det := a_cross_b.Dot(c.Vector)
	if det > maxDetError {
		return 1
	}
	if det < -maxDetError {
		return -1
	}
	return 0
}

func OrderedCCW(a, b, c, o Point) bool {
	// The last inequality below is ">" rather than ">=" so that we return
	// true if A == B or B == C, and otherwise false if A == C. Recall that
	// RobustCCW(x,y,z) == -RobustCCW(z,y,x) for all x,y,z.
	sum := 0
	if RobustCCW(b, o, a) >= 0 {
		sum++
	}
	if RobustCCW(c, o, b) >= 0 {
		sum++
	}
	if RobustCCW(a, o, c) > 0 {
		sum++
	}
	return sum >= 2
}

// Distance returns the angle between two points.
func (p Point) Distance(b Point) s1.Angle {
	return p.Vector.Angle(b.Vector)
}

// ApproxEqual reports if the two points are similar enough to be equal.
func (p Point) ApproxEqual(other Point) bool {
	const epsilon = 1e-14
	return p.Vector.Angle(other.Vector) <= epsilon
}

// ApproxEqualWithin reports if the two points are similar enough to be equal.
func (p Point) ApproxEqualWithin(other Point, maxError float64) bool {
	return float64(p.Vector.Angle(other.Vector)) <= maxError
}

func TurnAngle(a, b, c Point) float64 {
	angle := b.PointCross(a).Vector.Angle(c.PointCross(b).Vector)
	if RobustCCW(a, b, c) > 0 {
		return float64(angle)
	}
	return float64(-angle)
}

// PointArea returns the area on the unit sphere for the triangle defined by the
// given points.
//
// This method is based on l'Huilier's theorem,
//
//   tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
//
// where E is the spherical excess of the triangle (i.e. its area),
//       a, b, c are the side lengths, and
//       s is the semiperimeter (a + b + c) / 2.
//
// The only significant source of error using l'Huilier's method is the
// cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
// *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
// to a relative error of about 1e-15 / E using Girard's formula, where E is
// the true area of the triangle. Girard's formula can be even worse than
// this for very small triangles, e.g. a triangle with a true area of 1e-30
// might evaluate to 1e-5.
//
// So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
// dmin = min(s-a, s-b, s-c). This basically includes all triangles
// except for extremely long and skinny ones.
//
// Since we don't know E, we would like a conservative upper bound on
// the triangle area in terms of s and dmin. It's possible to show that
// E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
// Using this, it's easy to show that we should always use l'Huilier's
// method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
// if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
// k3 is about 0.1. Since the best case error using Girard's formula
// is about 1e-15, this means that we shouldn't even consider it unless
// s >= 3e-4 or so.
func PointArea(a, b, c Point) float64 {
	sa := float64(b.Angle(c.Vector))
	sb := float64(c.Angle(a.Vector))
	sc := float64(a.Angle(b.Vector))
	s := 0.5 * (sa + sb + sc)
	if s >= 3e-4 {
		// Consider whether Girard's formula might be more accurate.
		dmin := s - math.Max(sa, math.Max(sb, sc))
		if dmin < 1e-2*s*s*s*s*s {
			// This triangle is skinny enough to use Girard's formula.
			ab := a.PointCross(b)
			bc := b.PointCross(c)
			ac := a.PointCross(c)
			area := math.Max(0.0, float64(ab.Angle(ac.Vector)-ab.Angle(bc.Vector)+bc.Angle(ac.Vector)))

			if dmin < s*0.1*area {
				return area
			}
		}
	}

	// Use l'Huilier's formula.
	return 4 * math.Atan(math.Sqrt(math.Max(0.0, math.Tan(0.5*s)*math.Tan(0.5*(s-sa))*
		math.Tan(0.5*(s-sb))*math.Tan(0.5*(s-sc)))))
}

func SignedArea(a, b, c Point) interface{} {
	return PointArea(a, b, c) * float64(RobustCCW(a, b, c))
}

func TrueCentroid(a, b, c Point) interface{} {
	angleA := b.Angle(c.Vector).Radians()
	angleB := c.Angle(a.Vector).Radians()
	angleC := a.Angle(b.Vector).Radians()
	ra, rb, rc := 1., 1., 1.
	if angleA != 0 {
		ra = angleA / math.Sin(angleA)
	}
	if angleB != 0 {
		rb = angleB / math.Sin(angleB)
	}
	if angleC != 0 {
		rc = angleC / math.Sin(angleC)
	}

	x := r3.Vector{a.X, b.X - a.X, c.X - a.X}
	y := r3.Vector{a.Y, b.Y - a.Y, c.Y - a.Y}
	z := r3.Vector{a.Z, b.Z - a.Z, c.Z - a.Z}
	r := r3.Vector{ra, rb - ra, rc - ra}
	v := r3.Vector{
		y.Cross(z).Dot(r),
		z.Cross(x).Dot(r),
		x.Cross(y).Dot(r),
	}
	return Point{v.Mul(0.5)}
}

func FrameFromPoint(z Point) (m r3.Matrix) {
	m.SetCol(2, z.Vector)
	m.SetCol(1, z.Ortho())
	m.SetCol(0, m.Col(1).Cross(z.Vector)) // Already unit-length.
	return
}

func PointFromFrame(m r3.Matrix, q Point) Point {
	return Point{m.MulVector(q.Vector)}
}

// TODO(dnadasi):
//   - Other CCW methods
//   - Area methods
//   - Centroid methods
