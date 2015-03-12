package s2

import (
	"container/heap"
	"sort"

	"github.com/davidreynolds/gos2/s1"
)

type LoopMap map[*Loop][]*Loop

type Polygon struct {
	loops       []*Loop
	bound       Rect
	hasHoles    bool
	numVertices int
}

func NewPolygonFromLoops(loops []*Loop) *Polygon {
	p := &Polygon{
		bound:    EmptyRect(),
		hasHoles: false,
	}
	p.Init(loops)
	return p
}

func NewPolygonFromLoop(loop *Loop) *Polygon {
	p := &Polygon{
		bound:       loop.Bound(),
		hasHoles:    false,
		numVertices: len(loop.vertices),
	}
	p.loops = append(p.loops, loop)
	return p
}

func NewPolygonFromCell(cell Cell) *Polygon {
	p := &Polygon{
		bound:       EmptyRect(),
		numVertices: 4,
	}
	loop := NewLoopFromCell(cell)
	p.bound = loop.Bound()
	p.loops = append(p.loops, loop)
	return p
}

func ContainsChild(a, b *Loop, loopMap LoopMap) bool {
	if a == b {
		return true
	}
	children := loopMap[a]
	for _, child := range children {
		if ContainsChild(child, b, loopMap) {
			return true
		}
	}
	return false
}

func (p *Polygon) NumLoops() int    { return len(p.loops) }
func (p *Polygon) Loop(k int) *Loop { return p.loops[k] }

func (p *Polygon) Init(loops []*Loop) {
	p.loops = make([]*Loop, len(loops))
	copy(p.loops, loops)
	for _, loop := range p.loops {
		p.numVertices += len(loop.vertices)
	}

	loopMap := LoopMap{}
	for _, loop := range p.loops {
		p.InsertLoop(loop, nil, loopMap)
	}

	// Reorder the loops in depth-first order.
	p.loops = p.loops[:0]
	p.InitLoop(nil, -1, loopMap)

	p.hasHoles = false
	p.bound = EmptyRect()

	for _, loop := range p.loops {
		if loop.Sign() < 0 {
			p.hasHoles = true
		} else {
			p.bound = p.bound.Union(loop.bound)
		}
	}
}

func (p *Polygon) CapBound() Cap {
	return p.bound.CapBound()
}

func (p *Polygon) RectBound() Rect {
	return p.bound
}

func (a *Polygon) ContainsOrCrosses(b *Loop) int {
	inside := false
	for _, loop := range a.loops {
		result := loop.ContainsOrCrosses(b)
		if result < 0 {
			return -1
		}
		if result > 0 {
			inside = inside != true
		}
	}
	if inside {
		return 1
	}
	return 0
}

func (a *Polygon) AnyLoopContains(b *Loop) bool {
	for _, loop := range a.loops {
		if loop.ContainsLoop(b) {
			return true
		}
	}
	return false
}

func (a *Polygon) ContainsAllShells(b *Polygon) bool {
	for _, loop := range b.loops {
		if loop.Sign() < 0 {
			continue
		}
		if a.ContainsOrCrosses(loop) <= 0 {
			return false
		}
	}
	return true
}

func (a *Polygon) ExcludesAllHoles(b *Polygon) bool {
	for _, loop := range b.loops {
		if loop.Sign() > 0 {
			continue
		}
		if a.ContainsOrCrosses(loop) != 0 {
			return false
		}
	}
	return true
}

func (a *Polygon) IntersectsAnyShell(b *Polygon) bool {
	// Return true if this polygon (A) intersects any shell of B.
	for _, loop := range b.loops {
		if loop.Sign() < 0 {
			continue
		}
		if a.IntersectsShell(loop) {
			return true
		}
	}
	return false
}

func (a *Polygon) IntersectsShell(b *Loop) bool {
	inside := false
	for _, loop := range a.loops {
		if loop.ContainsLoop(b) {
			inside = inside != true
		} else if !b.ContainsLoop(loop) && loop.Intersects(b) {
			// We definitely have an intersection if the loops
			// intersect AND one is not properly contained in the
			// other. If A is properly contained in a loop of B,
			// we don't know yet if it may be actually inside a
			// hole within B.
			return true
		}
	}
	return inside
}

func (a *Polygon) ContainsPolygon(b *Polygon) bool {
	if len(a.loops) == 1 && len(b.loops) == 1 {
		return a.loops[0].ContainsLoop(b.loops[0])
	}
	if !a.bound.ContainsRect(b.bound) {
		if !a.bound.Lng.Union(b.bound.Lng).IsFull() {
			return false
		}
	}
	if !a.hasHoles && !b.hasHoles {
		for _, loop := range b.loops {
			if !a.AnyLoopContains(loop) {
				return false
			}
		}
		return true
	}
	return a.ContainsAllShells(b) && b.ExcludesAllHoles(a)
}

func (a *Polygon) ContainsPoint(p Point) bool {
	if len(a.loops) == 1 {
		return a.loops[0].Contains(p)
	}
	if !a.bound.ContainsPoint(p) {
		return false
	}
	inside := false
	for _, loop := range a.loops {
		inside = inside != loop.Contains(p)
		if inside && !a.hasHoles { // Shells are disjoint.
			break
		}
	}
	return inside
}

func (p *Polygon) ContainsCell(cell Cell) bool {
	if len(p.loops) == 1 {
		return p.loops[0].ContainsCell(cell)
	}
	if !p.bound.ContainsPoint(cell.Center()) {
		return false
	}
	loop := NewLoopFromCell(cell)
	poly := NewPolygonFromLoop(loop)
	return p.ContainsPolygon(poly)
}

func (p *Polygon) MayIntersect(cell Cell) bool {
	if len(p.loops) == 1 {
		return p.loops[0].MayIntersect(cell)
	}
	if !p.bound.Intersects(cell.RectBound()) {
		return false
	}
	loop := NewLoopFromCell(cell)
	poly := NewPolygonFromLoop(loop)
	return p.Intersects(poly)
}

func (a *Polygon) Intersects(b *Polygon) bool {
	if len(a.loops) == 1 && len(b.loops) == 1 {
		return a.loops[0].Intersects(b.loops[0])
	}
	if !a.bound.Intersects(b.bound) {
		return false
	}
	if !a.hasHoles && !b.hasHoles {
		for _, l1 := range a.loops {
			for _, l2 := range b.loops {
				if l1.Intersects(l2) {
					return true
				}
			}
		}
		return false
	}

	// Otherwise if any shell of B is contained by an odd number of
	// loops of A, or any shell of A is contained by an odd number of
	// loops of B, or there is an intersection without containment,
	// then there is an intersection.
	return a.IntersectsAnyShell(b) || b.IntersectsAnyShell(a)
}

func (p *Polygon) Release(loops *[]*Loop) {
	if loops != nil {
		*loops = append(*loops, p.loops...)
	}
	p.loops = []*Loop{}
	p.bound = EmptyRect()
	p.hasHoles = false
}

func (p *Polygon) InitLoop(loop *Loop, depth int, loopMap LoopMap) {
	if loop != nil {
		loop.depth = depth
		p.loops = append(p.loops, loop)
	}
	for _, child := range loopMap[loop] {
		p.InitLoop(child, depth+1, loopMap)
	}
}

func (p *Polygon) InsertLoop(newLoop, parent *Loop, loopMap LoopMap) {
	for _, child := range loopMap[parent] {
		if child.ContainsNested(newLoop) {
			p.InsertLoop(newLoop, child, loopMap)
			return
		}
	}
	// Some of the children of the parent loop may now be children of the
	// new loop.
	for i := 0; i < len(loopMap[parent]); {
		child := loopMap[parent][i]
		if newLoop.ContainsNested(child) {
			loopMap[newLoop] = append(loopMap[newLoop], child)
			loopMap[parent] = append(loopMap[parent][:i], loopMap[parent][i+1:]...)
		} else {
			i++
		}
	}
	loopMap[parent] = append(loopMap[parent], newLoop)
}

type PointPair struct {
	first, second Point
}

type IntPair struct {
	first, second int
}

func AreLoopsValid(loops []*Loop) bool {
	// If a loop contains an edge AB, then no other loop may contain
	// AB or BA.
	if len(loops) > 1 {
		edges := make(map[PointPair]IntPair)
		for i, loop := range loops {
			for j := 0; j < len(loop.vertices); j++ {
				key := PointPair{*loop.vertex(j), *loop.vertex(j + 1)}
				if _, ok := edges[key]; !ok {
					edges[key] = IntPair{i, j}
					continue
				}
				return false
			}
		}
	}

	// Verify that no loop covers more than half of the sphere, and that
	// no two loops cross.
	for i, loop := range loops {
		if !loop.IsNormalized() {
			return false
		}
		for j := i + 1; j < len(loops); j++ {
			// This test not only checks for edge crossings, it
			// also detects cases where the two boundaries cross
			// at a shared vertex.
			if loop.ContainsOrCrosses(loops[j]) < 0 {
				return false
			}
		}
	}
	return true
}

func (p *Polygon) Parent(k int) int {
	depth := p.Loop(k).depth
	if depth == 0 {
		return -1
	}
	k -= 1
	for k >= 0 && p.Loop(k).depth >= depth {
		k--
	}
	return k
}

func (p *Polygon) IsNormalized() bool {
	var lastParent *Loop
	vertices := make(map[Point]struct{})
	for i := 0; i < p.NumLoops(); i++ {
		child := p.Loop(i)
		if child.depth == 0 {
			continue
		}
		parent := p.Loop(p.Parent(i))
		if parent != lastParent {
			vertices = make(map[Point]struct{})
			for j := 0; j < parent.NumVertices(); j++ {
				vertices[*parent.vertex(j)] = struct{}{}
			}
			lastParent = parent
		}
		count := 0
		for j := 0; j < child.NumVertices(); j++ {
			if _, ok := vertices[*child.vertex(j)]; ok {
				count++
			}
		}
		if count > 1 {
			return false
		}
	}
	return true
}

var intersectionTolerance = s1.Angle(1.5e-15)

func (p *Polygon) InitToIntersection(a, b *Polygon) {
	p.InitToIntersectionSloppy(a, b, intersectionTolerance)
}

func (p *Polygon) InitToIntersectionSloppy(a, b *Polygon, vertexMergeRadius s1.Angle) {
	if !a.bound.Intersects(b.bound) {
		return
	}

	// We want the boundary of A clipped to the interior of B,
	// plus the boundary of B clipped to the interior of A,
	// plus one copy of any directed edges that are in both boundaries.
	options := DIRECTED_XOR()
	options.vertex_merge_radius = vertexMergeRadius
	builder := NewPolygonBuilder(options)
	ClipBoundary(a, false, b, false, false, true, &builder)
	ClipBoundary(b, false, a, false, false, false, &builder)
	if !builder.AssemblePolygon(p, nil) {
		panic("Bad directed edges in InitToIntersection")
	}
}

func (p *Polygon) InitToUnion(a, b *Polygon) {
	p.InitToUnionSloppy(a, b, intersectionTolerance)
}

func (p *Polygon) InitToUnionSloppy(a, b *Polygon, vertexMergeRadius s1.Angle) {
	// We want the boundary of A clipped to the exterior of B,
	// plus the boundary of B clipped to the exterior of A,
	// plus one copy of any directed edges that are in both boundaries.
	options := DIRECTED_XOR()
	options.vertex_merge_radius = vertexMergeRadius
	builder := NewPolygonBuilder(options)
	ClipBoundary(a, false, b, false, true, true, &builder)
	ClipBoundary(b, false, a, false, true, false, &builder)
	if !builder.AssemblePolygon(p, nil) {
		panic("Bad directed edges in InitToUnion")
	}
}

func (p *Polygon) InitToDifference(a, b *Polygon) {
	p.InitToDifferenceSloppy(a, b, intersectionTolerance)
}

func (p *Polygon) InitToDifferenceSloppy(a, b *Polygon, vertexMergeRadius s1.Angle) {
	// We want the boundary of A clipped to the exterior of B,
	// plus the reversed boundary of B clipped to the interior of A,
	// plus one copy of any edge in A that is also a reverse edge in B.
	options := DIRECTED_XOR()
	options.vertex_merge_radius = vertexMergeRadius
	builder := NewPolygonBuilder(options)
	ClipBoundary(a, false, b, true, true, true, &builder)
	ClipBoundary(b, true, a, false, false, false, &builder)
	if !builder.AssemblePolygon(p, nil) {
		panic("Bad directed edges in InitToDifference")
	}
}

func (p *Polygon) InternalClipPolyline(invert bool, a *Polyline, out *[]*Polyline, mergeRadius s1.Angle) {
	// Clip the polyline A to the interior of this polygon.
	// The resulting polyline(s) will be appended to "out".
	// If invert is true, we clip A to the exterior of this polygon
	// instead. Vertices will be dropped such that adjacent vertices
	// will not be closer than "mergeRadius".
	//
	// We do the intersection/subtraction by walking the polyline edges.
	// For each edge, we compute all intersections with the polygon
	// boundary and sort them in increasing order of distance along that
	// edge. We then divide the intersection points into pairs, and
	// output a clipped polyline segment for each one.
	// We keep track of whether we're inside or outside of the polygon
	// at all times to decide which segments to output.
	var intersections IntersectionSet
	var vertices []Point
	index := NewPolygonIndex(p, false)
	n := a.NumVertices()
	inside := p.ContainsPoint(a.Vertex(0)) != invert
	for j := 0; j < n-1; j++ {
		a0 := a.Vertex(j)
		a1 := a.Vertex(j + 1)
		ClipEdge(a0, a1, index, true, &intersections)
		if inside {
			intersections = append(intersections, FloatPointPair{0, a0})
		}
		inside = (len(intersections) & 1) != 0
		if inside {
			intersections = append(intersections, FloatPointPair{1, a1})
		}
		sort.Sort(intersections)
		// At this point we have a sorted array of vertex pairs
		// representing the edge(s) obtained after clipping (a0,a1)
		// against the polygon.
		for k := 0; k < len(intersections); k += 2 {
			if intersections[k] == intersections[k+1] {
				continue
			}
			v0 := intersections[k].second
			v1 := intersections[k+1].second

			// If the gap from the previous vertex to this one is
			// large enough, start a new polyline.
			if len(vertices) > 0 && vertices[len(vertices)-1].Angle(v0.Vector).Radians() > mergeRadius.Radians() {
				*out = append(*out, PolylineFromPoints(vertices))
				vertices = []Point{}
			}
			// Append this segment to the current polyline,
			// ignoring any vertices that are too close to the
			// previous vertex.
			if len(vertices) == 0 {
				vertices = append(vertices, v0)
			}
			if vertices[len(vertices)-1].Angle(v1.Vector).Radians() > mergeRadius.Radians() {
				vertices = append(vertices, v1)
			}
		}
		intersections = IntersectionSet{}
	}
	if len(vertices) > 0 {
		*out = append(*out, PolylineFromPoints(vertices))
	}
}

func (p *Polygon) IntersectWithPolyline(a *Polyline, out *[]*Polyline) {
	p.IntersectWithPolylineSloppy(a, out, intersectionTolerance)
}

func (p *Polygon) IntersectWithPolylineSloppy(a *Polyline, out *[]*Polyline, vertexMergeRadius s1.Angle) {
	p.InternalClipPolyline(false, a, out, vertexMergeRadius)
}

func (p *Polygon) SubtractFromPolyline(a *Polyline, out *[]*Polyline) {
	p.SubtractFromPolylineSloppy(a, out, intersectionTolerance)
}

func (p *Polygon) SubtractFromPolylineSloppy(a *Polyline, out *[]*Polyline, mergeRadius s1.Angle) {
	p.InternalClipPolyline(true, a, out, mergeRadius)
}

func DestructiveUnion(polygons *[]*Polygon) *Polygon {
	return DestructiveUnionSloppy(polygons, intersectionTolerance)
}

func DestructiveUnionSloppy(polygons *[]*Polygon, vertexMergeRadius s1.Angle) *Polygon {
	// Create a priority queue of polygons in order of number of vertices.
	// Repeatedly union the two smallest polygons and add the result to the
	// queue until we have a single polygon to return.
	pq := make(IntPolygonQueue, len(*polygons))
	for i := 0; i < len(*polygons); i++ {
		pq[i] = &IntPolygonPair{
			size: (*polygons)[i].numVertices,
			poly: (*polygons)[i],
		}
	}
	heap.Init(&pq)
	*polygons = []*Polygon{}
	for pq.Len() > 1 {
		// Pop two simplest polygons from queue.
		a := heap.Pop(&pq).(*IntPolygonPair)
		b := heap.Pop(&pq).(*IntPolygonPair)
		// Union and add result back to queue.
		var c Polygon
		c.InitToUnionSloppy(a.poly, b.poly, vertexMergeRadius)
		heap.Push(&pq, &IntPolygonPair{
			size: a.size + b.size,
			poly: &c,
		})
	}
	if pq.Len() == 0 {
		return &Polygon{}
	}
	return heap.Pop(&pq).(*IntPolygonPair).poly
}

type IntPolygonPair struct {
	size int
	poly *Polygon
	idx  int
}

type IntPolygonQueue []*IntPolygonPair

func (p IntPolygonQueue) Len() int           { return len(p) }
func (p IntPolygonQueue) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }
func (p IntPolygonQueue) Less(i, j int) bool { return p[i].size < p[j].size }

func (p *IntPolygonQueue) Push(x interface{}) {
	n := len(*p)
	pair := x.(*IntPolygonPair)
	pair.idx = n
	*p = append(*p, pair)
}

func (p *IntPolygonQueue) Pop() interface{} {
	old := *p
	n := len(old)
	pair := old[n-1]
	pair.idx = -1
	*p = old[0 : n-1]
	return pair
}

func (a *Polygon) BoundaryApproxEquals(b *Polygon, maxError float64) bool {
	if a.NumLoops() != b.NumLoops() {
		return false
	}
	for i := 0; i < a.NumLoops(); i++ {
		aLoop := a.Loop(i)
		success := false
		for j := 0; j < a.NumLoops(); j++ {
			bLoop := b.Loop(j)
			if bLoop.Depth() == aLoop.Depth() && bLoop.BoundaryApproxEquals(aLoop, maxError) {
				success = true
				break
			}
		}
		if !success {
			return false
		}
	}
	return true
}

func (a *Polygon) BoundaryNear(b *Polygon, maxError float64) bool {
	if a.NumLoops() != b.NumLoops() {
		return false
	}
	for i := 0; i < a.NumLoops(); i++ {
		aLoop := a.Loop(i)
		success := false
		for j := 0; j < a.NumLoops(); j++ {
			bLoop := b.Loop(j)
			if bLoop.Depth() == aLoop.Depth() && bLoop.BoundaryNear(aLoop, maxError) {
				success = true
				break
			}
		}
		if !success {
			return false
		}
	}
	return true
}

type FloatPointPair struct {
	first  float64
	second Point
}

type IntersectionSet []FloatPointPair

func (x IntersectionSet) Len() int      { return len(x) }
func (x IntersectionSet) Swap(i, j int) { x[i], x[j] = x[j], x[i] }
func (x IntersectionSet) Less(i, j int) bool {
	if x[i].first < x[j].first {
		return true
	}
	if x[j].first < x[i].first {
		return false
	}
	if x[i].second.LessThan(x[j].second.Vector) {
		return true
	}
	return false
}

func AddIntersection(a0, a1, b0, b1 Point, add_shared_edges bool, crossing int, intersections *IntersectionSet) {
	if crossing > 0 {
		// There is a proper edge crossing.
		x := GetIntersection(a0, a1, b0, b1)
		t := GetDistanceFraction(x, a0, a1)
		*intersections = append(*intersections, FloatPointPair{t, x})
	} else if VertexCrossing(a0, a1, b0, b1) {
		t := 1.0
		if a0 == b0 || a0 == b1 {
			t = 0.0
		}
		if !add_shared_edges && a1 == b1 {
			t = 1.0
		}
		pair := FloatPointPair{t, a0}
		if t != 0 {
			pair = FloatPointPair{t, a1}
		}
		*intersections = append(*intersections, pair)
	}
}

func ClipEdge(a0, a1 Point, bIndex *PolygonIndex, add_shared_edges bool, intersections *IntersectionSet) {
	it := NewEdgeIndexIterator(bIndex)
	it.GetCandidates(a0, a1)
	crosser := NewEdgeCrosser(&a0, &a1, &a0)
	var from, to *Point
	for ; !it.Done(); it.Next() {
		prev := to
		from, to = bIndex.EdgeFromTo(it.Index())
		if prev != from {
			crosser.RestartAt(from)
		}
		crossing := crosser.RobustCrossing(to)
		if crossing < 0 {
			continue
		}
		AddIntersection(a0, a1, *from, *to, add_shared_edges, crossing, intersections)
	}
}

func ClipBoundary(a *Polygon, reverse_a bool, b *Polygon, reverse_b bool, invert_b, add_shared_edges bool, builder *PolygonBuilder) {
	bIndex := NewPolygonIndex(b, reverse_b)
	PredictAdditionalCalls(bIndex, a.numVertices)
	for _, aLoop := range a.loops {
		n := aLoop.NumVertices()
		dir := 1
		if aLoop.IsHole() != reverse_a {
			dir = -1
		}
		inside := b.ContainsPoint(*aLoop.vertex(0)) != invert_b
		j := n
		if dir > 0 {
			j = 0
		}
		for ; n > 0; n, j = n-1, j+dir {
			a0 := *aLoop.vertex(j)
			a1 := *aLoop.vertex(j + dir)
			var intersections IntersectionSet
			ClipEdge(a0, a1, bIndex, add_shared_edges, &intersections)
			if inside {
				intersections = append(intersections, FloatPointPair{0, a0})
			}
			inside = (len(intersections) & 1) != 0
			if inside {
				intersections = append(intersections, FloatPointPair{1, a1})
			}
			sort.Sort(intersections)
			for k := 0; k < len(intersections); k += 2 {
				if intersections[k] == intersections[k+1] {
					continue
				}
				builder.AddEdge(intersections[k].second, intersections[k+1].second)
			}
		}
	}
}

type LoopSequenceIndexer interface {
	EdgeFromTo(i int) (*Point, *Point)
}

type LoopSequenceIndex struct {
	*EdgeIndex
	index_to_loop       []int
	loop_to_first_index []int
	num_edges           int
	num_loops           int
}

func NewLoopSequenceIndex() LoopSequenceIndex {
	return LoopSequenceIndex{
		EdgeIndex:           NewEdgeIndex(),
		index_to_loop:       []int{},
		loop_to_first_index: []int{},
	}
}

func (idx LoopSequenceIndex) NumEdges() int { return idx.num_edges }

func (idx LoopSequenceIndex) DecodeIndex(i int) (loopIndex, vertexInLoop int) {
	loopIndex = idx.index_to_loop[i]
	vertexInLoop = i - idx.loop_to_first_index[loopIndex]
	return
}

type PolygonIndex struct {
	LoopSequenceIndex
	poly    *Polygon
	reverse bool
}

func (idx *PolygonIndex) EdgeFromTo(i int) (from *Point, to *Point) {
	loopIndex, vertexInLoop := idx.DecodeIndex(i)
	loop := idx.poly.Loop(loopIndex)
	var fromIdx, toIdx int
	if loop.IsHole() != idx.reverse {
		fromIdx = loop.NumVertices() - 1 - vertexInLoop
		toIdx = 2*loop.NumVertices() - 2 - vertexInLoop
	} else {
		fromIdx = vertexInLoop
		toIdx = vertexInLoop + 1
	}
	from = loop.vertex(fromIdx)
	to = loop.vertex(toIdx)
	return
}

func NewPolygonIndex(poly *Polygon, reverse bool) *PolygonIndex {
	p := &PolygonIndex{
		LoopSequenceIndex: NewLoopSequenceIndex(),
		poly:              poly,
		reverse:           reverse,
	}
	for i := 0; i < p.poly.NumLoops(); i++ {
		p.AddLoop(p.poly.Loop(i).NumVertices())
	}
	return p
}

func (p *PolygonIndex) AddLoop(numVerts int) {
	verts_so_far := p.LoopSequenceIndex.num_edges
	p.LoopSequenceIndex.loop_to_first_index = append(p.LoopSequenceIndex.loop_to_first_index, verts_so_far)
	for i := 0; i < numVerts; i++ {
		p.LoopSequenceIndex.index_to_loop = append(p.LoopSequenceIndex.index_to_loop,
			p.LoopSequenceIndex.num_loops)
		verts_so_far++
	}
	p.LoopSequenceIndex.num_edges += numVerts
	p.LoopSequenceIndex.num_loops++
}
