package s2

import (
	"math"

	"github.com/davidreynolds/gos2/r1"
	"github.com/davidreynolds/gos2/s1"
)

// Indexing structure to efficiently compute intersections.
type LoopIndex struct {
	*EdgeIndex
	loop *Loop
}

func NewLoopIndex(loop *Loop) *LoopIndex {
	return &LoopIndex{
		EdgeIndex: NewEdgeIndex(),
		loop:      loop,
	}
}

func (idx *LoopIndex) EdgeFromTo(i int) (*Point, *Point) {
	return idx.loop.vertex(i), idx.loop.vertex(i + 1)
}

func (idx *LoopIndex) NumEdges() int {
	return len(idx.loop.vertices)
}

// A Loop represents a simple spherical polygon. It consists of a single
// chain of vertices where the first vertex is implicitly connected to the
// last. All loops are defined to have a CCW orientation, i.e. the interior
// of the polygon is on the left side of the edges. This implies that a
// clockwise loop enclosing a small area is interpreted to be a CCW loop
// enclosing a very large area.
//
// Loops are not allowed to have any duplicate vertices (whether adjacent or
// not), and non-adjacent edges are not allowed to intersect. Loops must have
// at least 3 vertices. Although these restrictions are not enforced in
// optimized code, you may get unexpected results if they are violated.
//
// Point containment is defined such that if the sphere is subdivided into
// faces (loops), every point is contained by exactly one face. This implies
// that loops do not necessarily contain all (or any) of their vertices.
type Loop struct {
	vertices      []Point
	bound         Rect
	origin_inside bool
	depth         int
	index         *LoopIndex

	// Map for speeding up FindVertex: We will compute a map from vertex to
	// index in the vertex array as soon as there have been enough calls.
	num_find_vertex_calls int
	vertex_to_index       map[Point]int
}

func NewLoopFromPath(vertices []Point) *Loop {
	loop := &Loop{
		vertices: make([]Point, len(vertices)),
		bound:    FullRect(),
		depth:    0,
		num_find_vertex_calls: 0,
	}
	loop.ResetMutableFields()
	copy(loop.vertices, vertices)
	loop.index = NewLoopIndex(loop)
	loop.InitOrigin()
	loop.InitBound()
	return loop
}

func NewLoopFromCell(cell Cell) *Loop {
	loop := &Loop{
		vertices: make([]Point, 4),
		bound:    cell.RectBound(),
		depth:    0,
		num_find_vertex_calls: 0,
	}
	for i := 0; i < 4; i++ {
		loop.vertices[i] = cell.Vertex(i)
	}
	loop.index = NewLoopIndex(loop)
	loop.InitOrigin()
	loop.InitBound()
	return loop
}

func (l *Loop) IsValid() bool {
	// Loops must have at least 3 vertices.
	if len(l.vertices) < 3 {
		return false
	}
	// All vertices must be unit length.
	for _, v := range l.vertices {
		if !v.IsUnit() {
			return false
		}
	}
	// Loops are not allowed to have any duplicate vertices.
	vmap := map[Point]int{}
	for i, v := range l.vertices {
		if _, ok := vmap[v]; ok {
			return false
		}
		vmap[v] = i
	}
	// Non-adjacent edges are not allowed to intersect.
	crosses := false
	PredictAdditionalCalls(l.index, len(l.vertices))
	it := NewEdgeIndexIterator(l.index)
	for i := 0; i < len(l.vertices); i++ {
		crosser := NewEdgeCrosser(l.vertex(i), l.vertex(i+1), l.vertex(0))
		prevIndex := -2
		for it.GetCandidates(*l.vertex(i), *l.vertex(i + 1)); !it.Done(); it.Next() {
			ai := it.Index()
			if ai > i+1 {
				if prevIndex != ai {
					crosser.RestartAt(l.vertex(ai))
				}
				crosses = crosser.RobustCrossing(l.vertex(ai+1)) > 0
				prevIndex = ai + 1
				if crosses {
					break
				}
			}
		}
		if crosses {
			break
		}
	}
	return !crosses
}

func (l *Loop) Depth() int { return l.depth }

func (l *Loop) IsHole() bool { return (l.depth & 1) != 0 }
func (l *Loop) Sign() int {
	if l.IsHole() {
		return -1
	}
	return 1
}

func (l *Loop) Clone() *Loop {
	loop := &Loop{
		vertices:      make([]Point, len(l.vertices)),
		bound:         l.bound,
		origin_inside: l.origin_inside,
		depth:         l.depth,
		num_find_vertex_calls: 0,
	}
	copy(loop.vertices, l.vertices)
	loop.index = NewLoopIndex(loop)
	return loop
}

func (l *Loop) Bound() Rect {
	return l.bound
}

func (l *Loop) CapBound() Cap {
	return l.bound.CapBound()
}

func (l *Loop) FindVertex(p Point) int {
	l.num_find_vertex_calls++
	if len(l.vertices) < 10 || l.num_find_vertex_calls < 20 {
		// Exhaustive search
		for i := 1; i <= len(l.vertices); i++ {
			if *l.vertex(i) == p {
				return i
			}
		}
		return -1
	}
	if len(l.vertex_to_index) == 0 { // We haven't computed it yet.
		for i := len(l.vertices); i > 0; i-- {
			l.vertex_to_index[*l.vertex(i)] = i
		}
	}

	if idx, ok := l.vertex_to_index[p]; ok {
		return idx
	}
	return -1
}

func (l *Loop) Invert() {
	l.ResetMutableFields()
	// Reverse vertices.
	for i, j := 0, len(l.vertices)-1; i < j; i, j = i+1, j-1 {
		l.vertices[i], l.vertices[j] = l.vertices[j], l.vertices[i]
	}
	l.origin_inside = !l.origin_inside
	if l.bound.Lo().Lat > -(math.Pi/2) && l.bound.Hi().Lat < math.Pi/2 {
		// The complement of this loop contains both poles.
		l.bound = FullRect()
	} else {
		l.InitBound()
	}
}

type f_tri func(a, b, c Point) interface{}

// GetSurfaceIntegral sums "f_tri" over a collection of oriented triangles,
// possibly overlapping.  Let the sign of a triangle be +1 if it is CCW and -1
// otherwise, and let the sign of a point "x" be the sum of the signs of the
// triangles containing "x".  Then the collection of triangles is chosen
// such that either:
//
//  (1) Each point in the loop interior has sign +1, and sign 0 otherwise; or
//  (2) Each point in the loop exterior has sign -1, and sign 0 otherwise.
//
// The triangles basically consist of a "fan" from vertex 0 to every loop
// edge that does not include vertex 0.  These triangles will always satisfy
// either (1) or (2).  However, what makes this a bit tricky is that
// spherical edges become numerically unstable as their length approaches
// 180 degrees.  Of course there is not much we can do if the loop itself
// contains such edges, but we would like to make sure that all the triangle
// edges under our control (i.e., the non-loop edges) are stable.  For
// example, consider a loop around the equator consisting of four equally
// spaced points.  This is a well-defined loop, but we cannot just split it
// into two triangles by connecting vertex 0 to vertex 2.
//
// We handle this type of situation by moving the origin of the triangle fan
// whenever we are about to create an unstable edge.  We choose a new
// location for the origin such that all relevant edges are stable.  We also
// create extra triangles with the appropriate orientation so that the sum
// of the triangle signs is still correct at every point.

// The maximum length of an edge for it to be considered numerically stable.
// The exact value is fairly arbitrary since it depends on the stability of
// the "f_tri" function.  The value below is quite conservative but could be
// reduced further if desired.
type Summer interface {
	Add(value interface{})
}

type FloatSummer float64
type PointSummer Point

func (s *FloatSummer) Add(value interface{}) {
	*s += FloatSummer(value.(float64))
}

func (s *PointSummer) Add(value interface{}) {
	s.X += value.(Point).X
	s.Y += value.(Point).Y
	s.Z += value.(Point).Z
}

func (l *Loop) GetSurfaceIntegral(fn f_tri, sum Summer) {
	const kMaxLength = math.Pi - 1e-5
	origin := *l.vertex(0)
	for i := 1; i+1 < len(l.vertices); i++ {
		// Let V_i be vertex(i), let O be the current origin, and let length(A,B)
		// be the length of edge (A,B).  At the start of each loop iteration, the
		// "leading edge" of the triangle fan is (O,V_i), and we want to extend
		// the triangle fan so that the leading edge is (O,V_i+1).
		//
		// Invariants:
		//  1. length(O,V_i) < kMaxLength for all (i > 1).
		//  2. Either O == V_0, or O is approximately perpendicular to V_0.
		//  3. "sum" is the oriented integral of f over the area defined by
		//     (O, V_0, V_1, ..., V_i).
		if (*l.vertex(i + 1)).Angle(origin.Vector).Radians() > kMaxLength {
			// We are about to create an unstable edge, so choose a new origin O'
			// for the triangle fan.
			oldOrigin := origin
			switch {
			case origin == *l.vertex(0):
				// The following point is well-separated from V_i and V_0
				// (and therefore V_i+1 as well).
				origin = Point{(*l.vertex(0)).PointCross(*l.vertex(i)).Normalize()}
			case (*l.vertex(i)).Angle((*l.vertex(0)).Vector).Radians() < kMaxLength:
				// All edges of the triangle (O, V_0, V_i) are stable,
				// so we can revert to using V_0 as the origin.
				origin = *l.vertex(0)
			default:
				// (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0
				// are perpendicular. Therefore V_0.Cross(0) is approximately
				// perpendicular to all of {O, V_0, V_i, V_i+1}, and we can
				// choose this point O' as the new origin.
				origin = Point{(*l.vertex(0)).Cross(oldOrigin.Vector)}
				// Advance the edge (O, V_i) to (O', V_i).
				sum.Add(fn(*l.vertex(0), oldOrigin, origin))
			}
			// Advance the edge (O, V_i) to (O', V_i).
			sum.Add(fn(oldOrigin, *l.vertex(i), origin))
		}
		// Advance the edge (O, V_i) to (O, V_i+1).
		sum.Add(fn(origin, *l.vertex(i), *l.vertex(i + 1)))
	}
	// If the origin is not V_0, we need to sum one more triangle.
	if origin != *l.vertex(0) {
		// Advance the edge (O, V_n-1) to (O, V_0).
		sum.Add(fn(origin, *l.vertex(len(l.vertices) - 1), *l.vertex(0)))
	}
}

func (l *Loop) Area() float64 {
	var sum FloatSummer
	l.GetSurfaceIntegral(SignedArea, &sum)
	area := float64(sum)
	if area < 0 {
		area += 4 * math.Pi
	}
	return math.Max(0.0, math.Min(4*math.Pi, area))
}

func (l *Loop) Centroid() Point {
	var point PointSummer
	l.GetSurfaceIntegral(TrueCentroid, &point)
	return Point(point)
}

func (l *Loop) ResetMutableFields() {
	if l.index != nil {
		l.index.Reset()
	}
	l.num_find_vertex_calls = 0
	l.vertex_to_index = map[Point]int{}
}

func (l *Loop) InitOrigin() {
	// The bounding box does not need to be correct before calling this
	// function, but it must at least contain vertex(1) since we need to
	// do a Contains() test on this point below.
	//l.bound.Contains(LatLngFromPoint(*l.vertex(1)))

	// To ensure that every point is contained in exactly on face of a
	// subdivision of the sphere, all containment tests are done by counting
	// the edge crossings starting at a fixed point on there sphere
	// (Origin()). We need to know whether this point is inside or outside
	// of the loop. We do this by first guessing that it is outside, and
	// then seeing whether we get the correct containment result for
	// vertex 1. If the result is incorrect, the origin must be inside the
	// loop.
	//
	// A loop with consecutive vertices A,B,C contains vertex B iff the
	// fixed vector R = Ortho(B) is on the left side of the wedge ABC.
	// The test below is written so that B is inside if C=R but not if A=R.
	l.origin_inside = false // Initialize before calling Contains()
	v1_inside := OrderedCCW(Point{l.vertex(1).Ortho()}, *l.vertex(0), *l.vertex(2), *l.vertex(1))
	if v1_inside != l.Contains(*l.vertex(1)) {
		l.origin_inside = true
	}
}

func (l *Loop) InitBound() {
	// The bounding rectangle of a loop is not necessarily the same as the
	// bounding rectangle of its vertices. First, the loop may wrap entirely
	// around the sphere (e.g. a loop that defines two revolutions of a
	// candy-cane stripe). Second, the loop may include one or both poles.
	// Note that a small clockwise loop near the equator contains both
	// poles.

	bounder := NewRectBounder()
	for i := 0; i <= len(l.vertices); i++ {
		bounder.AddPoint(l.vertex(i))
	}

	b := bounder.Bound()
	// Note that we need to initialize l.bound with a temporary value since
	// Contains() does a bounding rectangle check before doing anything
	// else.
	l.bound = FullRect()
	if l.Contains(PointFromCoords(0, 0, 1)) {
		b = Rect{
			Lat: r1.Interval{b.Lat.Lo, math.Pi / 2},
			Lng: s1.FullInterval(),
		}
	}
	// If a loop contains the south pole, then either it wraps entirely
	// around the sphere (full longitude range), or it also contains the
	// north pole in which case b.Lng.IsFull() due to the test above.
	// Either way, we only need to do the south pole containment test if
	// b.Lng.IsFull()
	if b.Lng.IsFull() && l.Contains(PointFromCoords(0, 0, -1)) {
		b.Lat.Lo = -math.Pi / 2
	}
	l.bound = b
}

func (l *Loop) NumVertices() int    { return len(l.vertices) }
func (l *Loop) Vertex(i int) *Point { return l.vertex(i) }

func (l *Loop) vertex(i int) *Point {
	j := i - l.NumVertices()
	if j >= 0 {
		return &l.vertices[j]
	}
	return &l.vertices[i]
}

func (l *Loop) IsNormalized() bool {
	// Optimization: if the longitude span is less than 180 degrees, then
	// the loop covers less than half the sphere and is therefore
	// normalized.
	if l.bound.Lng.Length() < math.Pi {
		return true
	}
	// We allow some error so that hemispheres are always considered
	// normalized.
	return l.TurningAngle() >= -1e-14
}

func (l *Loop) Normalize() {
	if !l.IsNormalized() {
		l.Invert()
	}
}

// Return (first, dir) such that first..first+n*dir are valid indices.
func (l *Loop) CanonicalFirstVertex() (first, dir int) {
	first = 0
	n := l.NumVertices()
	for i := 1; i < n; i++ {
		if l.vertex(i).LessThan(l.vertex(first).Vector) {
			first = i
		}
	}
	if l.vertex(first + 1).LessThan(l.vertex(first + n - 1).Vector) {
		dir = 1
		// 0 <= first <= n-1, so (first+n*dir) <= 2*n-1.
	} else {
		dir = -1
		first += n
		// n <= first <= 2*n-1, so (first+n*dir) >= 0.
	}
	return
}

func (l *Loop) TurningAngle() float64 {
	// Don't crash even if the loop is not well-defined.
	if l.NumVertices() < 3 {
		return 0
	}
	// To ensure that we get the same result when the loop vertex order is
	// rotated, and that we get the same result with the opposite sign when
	// the vertices are reversed, we need to be careful to add up the
	// individual turn angles in a consistent order. In general, adding up
	// a set of numbers in a different order can change the sum due to
	// rounding errors.
	n := l.NumVertices()
	i, dir := l.CanonicalFirstVertex()
	angle := TurnAngle(*l.vertex((i + n - dir) % n), *l.vertex(i), *l.vertex((i + dir) % n))
	for n = n - 1; n > 0; n-- {
		i += dir
		angle += TurnAngle(*l.vertex(i - dir), *l.vertex(i), *l.vertex(i + dir))
	}
	return float64(dir) * angle
}

func (l *Loop) ContainsCell(cell Cell) bool {
	if !l.bound.ContainsPoint(cell.Center()) {
		return false
	}
	cell_loop := NewLoopFromCell(cell)
	return l.ContainsLoop(cell_loop)
}

func (a *Loop) ContainsLoop(b *Loop) bool {
	// For this loop A to contain the given loop B, all of the following
	// must be true:
	//
	//  (1) There are no edge crossings between A and B except at vertices.
	//
	//  (2) At every vertex that is shared between A and B, the local edge
	//      ordering implies that A contains B.
	//
	//  (3) If there are no shared vertices, then A must contain a vertex
	//      of B and B must not contain a vertex of A. (An arbitrary vertex
	//      may be chosen in each case.)
	//
	// The second part of (3) is necessary to detect the case of two loops
	// whose union is the entire sphere, i.e. two loops that contain each
	// other's boundaries but not each other's interiors.
	if !a.bound.ContainsRect(b.bound) {
		return false
	}

	// Unless there are shared vertices, we need to check whether A
	// contains a vertex of B. Since shared vertices are rare, it is more
	// efficient to do this test up front as a quick rejection test.
	if !a.Contains(*b.vertex(0)) && a.FindVertex(*b.vertex(0)) < 0 {
		return false
	}

	// Now check whether there are any edge crossings, and also check the
	// loop relationship at any shared vertices.
	var wedge ContainsWedgeProcessor
	if a.AreBoundariesCrossing(b, &wedge) || wedge.doesnt_contain {
		return false
	}

	// At this point we know that the boundaries of A and B do not
	// intersect, and that A contains a vertex of B. However, we still
	// need to check for the case mentioned above, where (A union B) is
	// the entire sphere. Normally this check is very cheap due to the
	// bounding box precondition.
	if a.bound.Union(b.bound).IsFull() {
		if b.Contains(*a.vertex(0)) && b.FindVertex(*a.vertex(0)) < 0 {
			return false
		}
	}
	return true
}

func (a *Loop) ContainsNested(b *Loop) bool {
	if !a.bound.ContainsRect(b.bound) {
		return false
	}
	// We are given that A and B do not share any edges, and that either
	// one loop contains the other or they do not intersect.
	m := a.FindVertex(*b.vertex(1))
	if m < 0 {
		// Since b.vertex(1) is not shared, we can check whether A
		// contains it.
		return a.Contains(*b.vertex(1))
	}
	// Check whether the edge order around b.vertex(1) is compatible with
	// A containing B.
	return WedgeContains(*a.vertex(m - 1), *a.vertex(m), *a.vertex(m + 1),
		*b.vertex(0), *b.vertex(2))
}

func (l *Loop) ContainsPoint(p Point) bool { return l.Contains(p) }

func (l *Loop) Contains(p Point) bool {
	if !l.bound.Contains(LatLngFromPoint(p)) {
		return false
	}
	inside := l.origin_inside
	origin := OriginPoint()
	crosser := NewEdgeCrosser(&origin, &p, l.vertex(0))
	vlen := l.NumVertices()
	if vlen < 2000 {
		for i := 1; i <= vlen; i++ {
			inside = inside != crosser.EdgeOrVertexCrossing(l.vertex(i))
		}
		return inside
	}

	it := NewEdgeIndexIterator(l.index)
	prev_index := -2
	for it.GetCandidates(origin, p); !it.Done(); it.Next() {
		ai := it.Index()
		if prev_index != ai-1 {
			crosser.RestartAt(l.vertex(ai))
		}
		prev_index = ai
		inside = inside != crosser.EdgeOrVertexCrossing(l.vertex(ai+1))
	}
	return inside
}

func (l *Loop) MayIntersect(cell Cell) bool {
	if !l.bound.Intersects(cell.RectBound()) {
		return false
	}
	return NewLoopFromCell(cell).Intersects(l)
}

type WedgeProcessor interface {
	ProcessWedge(a0, ab1, a2, b0, b2 Point) bool
}

// WedgeProcessor to be used to check if loop A intersects loop B.
// Intersects() then returns true when A and B have at least one pair
// of associated wedges that intersect.
type IntersectsWedgeProcessor struct {
	intersects bool
}

func (p *IntersectsWedgeProcessor) ProcessWedge(a0, ab1, a2, b0, b2 Point) bool {
	p.intersects = WedgeIntersects(a0, ab1, a2, b0, b2)
	return p.intersects
}

// WedgeProcessor to be used to check if loop A contains loop B.
// DoesntContain() then returns true if there is a wedge of B not contained
// in the associated wedge of A (and hence loop B is not contained in loop A).
type ContainsWedgeProcessor struct {
	doesnt_contain bool
}

func (p *ContainsWedgeProcessor) ProcessWedge(a0, ab1, a2, b0, b2 Point) bool {
	p.doesnt_contain = !WedgeContains(a0, ab1, a2, b0, b2)
	return p.doesnt_contain
}

// WedgeProcessor to be used to check if the interior of loop A contains the
// interior of loop B, or their boundaries cross each other (therefore they
// have a proper intersection). CrossesOrMayContain() then returns -1 if A
// crossed B, 0 if it is not possible for A to contain B, and 1 otherwise.
type ContainsOrCrossesProcessor struct {
	// True if any crossing on the boundary is discovered.
	has_boundary_crossing bool
	// True if A (B) has a strictly superwedge on a pair of wedges that
	// share a common center point.
	a_has_strictly_super_wedge bool
	b_has_strictly_super_wedge bool
	// True if there is a pair of disjoint wedges with a common center
	// point.
	has_disjoint_wedge bool
}

func (p *ContainsOrCrossesProcessor) CrossesOrMayContain() int {
	if p.has_boundary_crossing {
		return -1
	}
	if p.has_disjoint_wedge || p.b_has_strictly_super_wedge {
		return 0
	}
	return 1
}

func (p *ContainsOrCrossesProcessor) ProcessWedge(a0, ab1, a2, b0, b2 Point) bool {
	wedgeRelation := GetWedgeRelation(a0, ab1, a2, b0, b2)
	if wedgeRelation == WEDGE_PROPERLY_OVERLAPS {
		p.has_boundary_crossing = true
		return true
	}
	p.a_has_strictly_super_wedge = p.a_has_strictly_super_wedge || (wedgeRelation == WEDGE_PROPERLY_CONTAINS)
	p.b_has_strictly_super_wedge = p.b_has_strictly_super_wedge || (wedgeRelation == WEDGE_IS_PROPERLY_CONTAINED)
	if p.a_has_strictly_super_wedge && p.b_has_strictly_super_wedge {
		p.has_boundary_crossing = true
		return true
	}
	p.has_disjoint_wedge = p.has_disjoint_wedge || (wedgeRelation == WEDGE_IS_DISJOINT)
	return false
}

// This method checks all edges of this loop (A) for intersection against
// all edges of B. If there is any shared vertex, the wedges centered at this
// vertex are sent to wedge_processor.
//
// Returns true only when the edges intersect in the sense of RobustCrossing,
// returns false immediately when the wedge_processor returns true: this means
// the wedge processor knows the value of the property that the caller wants
// to compute, and no further inspection is needed. For instance, if the
// property is "loops intersect", then a wedge intersection is all it takes
// to return true.
//
// See Intersects().
func (a *Loop) AreBoundariesCrossing(b *Loop, wedge_processor WedgeProcessor) bool {
	PredictAdditionalCalls(a.index, len(b.vertices))
	it := NewEdgeIndexIterator(a.index)
	for j := 0; j < len(b.vertices); j++ {
		crosser := NewEdgeCrosser(b.vertex(j), b.vertex(j+1), b.vertex(0))
		prev_index := -2
		for it.GetCandidates(*b.vertex(j), *b.vertex(j + 1)); !it.Done(); it.Next() {
			ai := it.Index()
			if prev_index != ai-1 {
				crosser.RestartAt(a.vertex(ai))
			}
			prev_index = ai
			crossing := crosser.RobustCrossing(a.vertex(ai + 1))
			if crossing < 0 {
				continue
			}
			if crossing > 0 {
				return true
			}
			// We only need to check each shared vertex once, so we
			// only consider the case where
			// a.vertex(ai+1) == b.vertex(j+1).
			if *a.vertex(ai + 1) == *b.vertex(j + 1) &&
				wedge_processor.ProcessWedge(*a.vertex(ai), *a.vertex(ai + 1), *a.vertex(ai + 2),
					*b.vertex(j), *b.vertex(j + 2)) {
				return false
			}
		}
	}
	return false
}

func (a *Loop) Intersects(b *Loop) bool {
	// a.Intersects(b) if and only if !a.Complement().Contains(b).
	// This code is similar to Contains(), but is optimized for the case
	// where both loops enclose less than half of the sphere.

	// The largest of the two loops should be edgeindex'd.
	if len(b.vertices) > len(a.vertices) {
		return b.Intersects(a)
	}

	if !a.bound.Intersects(b.bound) {
		return false
	}

	// Unless there are shared vertices, we need to check whether A
	// contains a vertex of B. Since shared vertices are rare, it is more
	// efficient to do this test up front as a quick acceptance test.
	if a.Contains(*b.vertex(0)) && a.FindVertex(*b.vertex(0)) < 0 {
		return true
	}

	// Now check whether there are any edge crossings, and also check the
	// loop relationship at any shared vertices.
	var wedge_processor IntersectsWedgeProcessor
	if a.AreBoundariesCrossing(b, &wedge_processor) || wedge_processor.intersects {
		return true
	}

	// We know that A does not contain a vertex of B, and that there are
	// no edge crossings. Therefore the only way that A can intersect B is
	// if B entirely contains A. We can check this by testing whether B
	// contains an arbitrary non-shared vertex of A. Note that this check
	// is usually cheap because of the bounding box precondition.
	if b.bound.ContainsRect(a.bound) {
		if b.Contains(*a.vertex(0)) && b.FindVertex(*a.vertex(0)) < 0 {
			return true
		}
	}
	return false
}

func (a *Loop) ContainsOrCrosses(b *Loop) int {
	// There can be containment or crossing only if the bounds intersect.
	if !a.bound.Intersects(b.bound) {
		return 0
	}
	// Now check whether there are any edge crossings, and also check
	// the loop relationship at any shared vertices. Note that unlike
	// Contains() or Intersects(), we can't do a point containment test
	// as a shortcut because we need to detect whether there are any
	// edge crossings.
	var wedgeProcessor ContainsOrCrossesProcessor
	if a.AreBoundariesCrossing(b, &wedgeProcessor) {
		return -1
	}
	res := wedgeProcessor.CrossesOrMayContain()
	if res <= 0 {
		return res
	}

	// At this point we know that the boundaries do not intersect, and we
	// are given that (A union B) is a proper subset of the sphere.
	// Furthermore, either A contains B or there are no shared vertices
	// (due to the check above). So now we just need to distinguish the
	// case where A contains B from the case where B contains A or the
	// two loops are disjoint.
	if !a.bound.ContainsRect(b.bound) {
		return 0
	}
	if !a.Contains(*b.vertex(0)) && a.FindVertex(*b.vertex(0)) < 0 {
		return 0
	}
	return 1
}

func (a *Loop) BoundaryApproxEquals(b *Loop, maxError float64) bool {
	numVerts := len(a.vertices)
	if numVerts != len(b.vertices) {
		return false
	}
	for offset := 0; offset < numVerts; offset++ {
		if a.vertex(offset).ApproxEqualWithin(*b.vertex(0), maxError) {
			success := true
			for i := 0; i < numVerts; i++ {
				if !a.vertex(i+offset).ApproxEqualWithin(*b.vertex(i), maxError) {
					success = false
					break
				}
			}
			if success {
				return true
			}
		}
	}
	return false
}

func (a *Loop) BoundaryNear(b *Loop, maxError float64) bool {
	for offset := 0; offset < len(a.vertices); offset++ {
		if MatchBoundaries(a, b, offset, maxError) {
			return true
		}
	}
	return false
}

func MatchBoundaries(a, b *Loop, offset int, maxError float64) bool {
	pending := []IntPair{}
	done := map[IntPair]bool{}
	pending = append(pending, IntPair{0, 0})
	alen := len(a.vertices)
	blen := len(b.vertices)
	for len(pending) > 0 {
		back := len(pending) - 1
		i := pending[back].first
		j := pending[back].second
		pending = pending[:back]
		if i == alen && j == blen {
			return true
		}
		done[IntPair{i, j}] = true

		io := i + offset
		if io >= alen {
			io -= alen
		}

		if i < alen {
			if _, ok := done[IntPair{i + 1, j}]; !ok {
				if a.vertex(io+1).DistanceToEdge(
					*b.vertex(j),
					*b.vertex(j + 1)).Radians() <= maxError {
					pending = append(pending, IntPair{i + 1, j})
				}
			}
		}
		if j < blen {
			if _, ok := done[IntPair{i, j + 1}]; !ok {
				if b.vertex(j+1).DistanceToEdge(
					*a.vertex(io),
					*a.vertex(io + 1)).Radians() <= maxError {
					pending = append(pending, IntPair{i, j + 1})
				}
			}
		}
	}
	return false
}
