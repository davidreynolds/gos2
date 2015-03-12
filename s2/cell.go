package s2

import (
	"math"

	"github.com/davidreynolds/gos2/r1"
	"github.com/davidreynolds/gos2/r2"
	"github.com/davidreynolds/gos2/s1"
)

const (
	maxError = 1.0 / (1 << 51)
)

var (
	poleMinLat = math.Asin(math.Sqrt(1. / 3))
)

// Cell is an S2 region object that represents a cell. Unlike CellIDs,
// it supports efficient containment and intersection tests. However, it is
// also a more expensive representation.
type Cell struct {
	face        int8
	level       int8
	orientation int8
	id          CellID
	uv          r2.Rect
}

// CellFromCellID constructs a Cell corresponding to the given CellID.
func CellFromCellID(id CellID) Cell {
	c := Cell{}
	c.id = id
	f, i, j, o := c.id.faceIJOrientation()
	c.face = int8(f)
	c.level = int8(c.id.Level())
	c.orientation = int8(o)
	c.uv = ijLevelToBoundUV(i, j, int(c.level))
	return c
}

// CellFromPoint constructs a cell for the given Point.
func CellFromPoint(p Point) Cell {
	return CellFromCellID(cellIDFromPoint(p))
}

// CellFromLatLng constructs a cell for the given LatLng.
func CellFromLatLng(ll LatLng) Cell {
	return CellFromCellID(CellIDFromLatLng(ll))
}

// IsLeaf returns whether this Cell is a leaf or not.
func (c Cell) IsLeaf() bool {
	return c.level == maxLevel
}

// SizeIJ returns the CellID value for the cells level.
func (c Cell) SizeIJ() int {
	return sizeIJ(int(c.level))
}

func (c Cell) Id() CellID           { return c.id }
func (c Cell) AverageArea() float64 { return AverageArea(int(c.level)) }

func (c Cell) ApproxArea() float64 {
	// All cells at the first two levels have the same area.
	if c.level < 2 {
		return c.AverageArea()
	}
	v0 := c.Vertex(0)
	v1 := c.Vertex(1)
	v2 := c.Vertex(2)
	v3 := c.Vertex(3)

	// First, compute the approximate area of the cell when projected
	// perpendicular to its normal. The cross product of its diagonals
	// gives the normal, and the length of the normal is twice the
	// projected area.
	flatArea := 0.5 * v2.Sub(v0.Vector).Cross(v3.Sub(v1.Vector)).Norm()

	// Now, compensate for the curvature of the cell surface by pretending
	// that the cell is shaped like a spherical cap. The ratio of the area
	// of a spherical cap to the area of its projected disc turns out to be
	// 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc. For
	// example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
	// Here we set Pi*r*r == flatArea to find the equivalent disc.
	return flatArea * 2 / (1 + math.Sqrt(1-math.Min(M_1_PI*flatArea, 1.0)))
}

func (c Cell) ExactArea() float64 {
	v0 := c.Vertex(0)
	v1 := c.Vertex(1)
	v2 := c.Vertex(2)
	v3 := c.Vertex(3)
	return PointArea(v0, v1, v2) + PointArea(v0, v2, v3)
}

func (c Cell) CapBound() Cap {
	// Use th cell center in (u,v)-space as the cap axis. This vector
	// is very close to Center() and faster to compute. Neither one of
	// these vectors yields the bounding cap with minimal surface area,
	// but they are both pretty close.
	//
	// It's possible to show that the two vertices that are furthest from
	// the (u,v)-origin never determine the maximum cap size (this is a
	// possible future optimization).
	u := 0.5 * (c.uv.X.Lo + c.uv.X.Hi)
	v := 0.5 * (c.uv.Y.Lo + c.uv.Y.Hi)
	s2cap := CapFromCenterHeight(Point{faceUVToXYZ(int(c.face), u, v).Normalize()}, 0)
	for k := 0; k < 4; k++ {
		s2cap.AddPoint(c.Vertex(k))
	}
	return s2cap
}

func (c Cell) Subdivide(children *[]Cell) bool {
	if c.IsLeaf() {
		return false
	}
	ci := c.id.ChildBegin()
	for i := 0; i < 4; i++ {
		*children = append(*children, CellFromCellID(ci))
		ci = ci.Next()
	}
	return true
}

func AverageArea(level int) float64 {
	return AvgArea.Value(level)
}

func (c Cell) ijToUV(i, j int) (u, v float64) {
	switch i {
	case 0:
		u = c.uv.X.Lo
	default:
		u = c.uv.X.Hi
	}
	switch j {
	case 0:
		v = c.uv.Y.Lo
	default:
		v = c.uv.Y.Hi
	}
	return
}

func (c Cell) Latitude(i, j int) float64 {
	u, v := c.ijToUV(i, j)
	p := Point{faceUVToXYZ(int(c.face), u, v)}
	return latitude(p).Radians()
}

func (c Cell) Longitude(i, j int) float64 {
	u, v := c.ijToUV(i, j)
	p := Point{faceUVToXYZ(int(c.face), u, v)}
	return longitude(p).Radians()
}

// Vertex returns the k-th vertex of the cell (k = [0,3]) in CCW order
// (lower left, lower right, upper right, upper left in the UV plane).
func (c Cell) Vertex(k int) Point {
	return Point{c.VertexRaw(k).Normalize()}
}

func (c Cell) VertexRaw(k int) Point {
	return Point{faceUVToXYZ(int(c.face), c.uv.Vertices()[k].X, c.uv.Vertices()[k].Y)}
}

// Edge returns the inward-facing normal of the great circle passing through
// the CCW ordered edge from vertex k to vertex k+1 (mod 4).
func (c Cell) Edge(k int) Point {
	return Point{c.EdgeRaw(k).Normalize()}
}

func (c Cell) EdgeRaw(k int) Point {
	switch k {
	case 0:
		return Point{vNorm(int(c.face), c.uv.Y.Lo)} // Bottom
	case 1:
		return Point{uNorm(int(c.face), c.uv.X.Hi)} // Right
	case 2:
		return Point{vNorm(int(c.face), c.uv.Y.Hi).Mul(-1.0)} // Top
	default:
		return Point{uNorm(int(c.face), c.uv.X.Lo).Mul(-1.0)} // Left
	}
}

func (c Cell) CenterRaw() Point {
	return Point{c.id.rawPoint()}
}

func (c Cell) Center() Point {
	return Point{c.CenterRaw().Normalize()}
}

func (c Cell) MayIntersect(cell Cell) bool {
	return c.id.Intersects(cell.id)
}

func (c Cell) ContainsCell(cell Cell) bool {
	return c.id.Contains(cell.id)
}

func (c Cell) ContainsPoint(p Point) bool {
	// We can't just call xyzToFaceUV, because for points that lie on the
	// boundary between two faces (i.e. u or v is +1/-1) we need to return
	// true for both adjacent cells.
	var u, v float64
	if !faceXYZToUV(int(c.face), p, &u, &v) {
		return false
	}
	return u >= c.uv.X.Lo && u <= c.uv.X.Hi &&
		v >= c.uv.Y.Lo && v <= c.uv.Y.Hi
}

func (c Cell) RectBound() Rect {
	if c.level > 0 {
		// Except for cells at level 0, the latitude and longitude
		// extremes are attained at the vertices. Furthermore, the
		// latitude range is determined by one pair of diagonally
		// opposite vertices and the longitude range is determined by
		// the other pair.
		//
		// We first determine which corner (i,j) of the cell has the
		// largest absolute latitude. To maximize latitude, we want to
		// find the point in the cell that has the largest absolute
		// z-coordinate and the smallest absolute x- and y-coordinates.
		// To do this we look at each coordinate (u and v), and
		// determine whether we want to minimize or maximize that
		// coordinate based on the axis direction and the cell's (u,v)
		// quadrant.
		u := c.uv.X.Lo + c.uv.X.Hi
		v := c.uv.Y.Lo + c.uv.Y.Hi
		i, j := ijFromFaceZ(c.face, u, v)

		// We grow the bounds slightly to make sure that the bounding
		// rectangle also contains the normalized versions of the
		// vertices. Note that the maximum result magnitude is Pi, with
		// a floating-point exponent of 1. Therefore adding or
		// subtracting 2**-51 will always change the result.
		lat := r1.IntervalFromPointPair(c.Latitude(i, j), c.Latitude(1-i, 1-j))
		lat = lat.Expanded(maxError).Intersection(validRectLatRange)
		if lat.Lo == -M_PI_2 || lat.Hi == M_PI_2 {
			return Rect{lat, s1.FullInterval()}
		}

		lng := s1.IntervalFromPointPair(c.Longitude(i, 1-j), c.Longitude(1-i, j))
		return Rect{lat, lng.Expanded(maxError)}
	}

	// The 4 cells around the equator extend to +/-45 degrees latitude at
	// the midpoints of their top and bottom edges. The two cells covering
	// the poles extend down to +/-35.26 degrees at their vertices.
	//
	// The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
	switch c.face {
	case 0:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{-M_PI_4, M_PI_4},
		}
	case 1:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{M_PI_4, 3 * M_PI_4},
		}
	case 2:
		return Rect{
			r1.Interval{poleMinLat, M_PI_2},
			s1.Interval{-math.Pi, math.Pi},
		}
	case 3:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{3 * M_PI_4, -3 * M_PI_4},
		}
	case 4:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{-3 * M_PI_4, -M_PI_4},
		}
	default:
		return Rect{
			r1.Interval{-M_PI_2, -poleMinLat},
			s1.Interval{-math.Pi, math.Pi},
		}
	}
}

func ijFromFaceZ(face int8, u, v float64) (i, j int) {
	switch uAxis(int(face)).Z {
	case 0:
		i = bool2int(u < 0)
	default:
		i = bool2int(u > 0)
	}
	switch vAxis(int(face)).Z {
	case 0:
		j = bool2int(v < 0)
	default:
		j = bool2int(v > 0)
	}
	return
}

func bool2int(b bool) int {
	if b {
		return 1
	}
	return 0
}

// TODO(roberts, or $SOMEONE): Differences from C++, almost everything else still.
// Implement the accessor methods on the internal fields.
