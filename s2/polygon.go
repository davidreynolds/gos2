package s2

import (
	"log"
	"sort"

	"github.com/davidreynolds/gos2/s1"
)

type LoopMap map[*Loop][]*Loop

type Polygon struct {
	loops       []*Loop
	bound       Rect
	ownsLoops   bool
	hasHoles    bool
	numVertices int
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

func (p Polygon) NumLoops() int    { return len(p.loops) }
func (p Polygon) Loop(k int) *Loop { return p.loops[k] }

func (p *Polygon) Init(loops *[]*Loop) {
	p.loops = make([]*Loop, len(*loops))
	if !AreLoopsValid(*loops) {
		log.Println("invalid loops")
	}
	copy(p.loops, *loops)
	for _, loop := range p.loops {
		p.numVertices += len(loop.vertices)
	}

	loopMap := LoopMap{}
	for _, loop := range p.loops {
		p.InsertLoop(loop, nil, loopMap)
	}

	// Reorder the loops in depth-first order.
	p.loops = []*Loop{}
	p.InitLoop(nil, -1, loopMap)

	for i := 0; i < len(p.loops); i++ {
		for j := 0; j < len(p.loops); j++ {
			if i == j {
				continue
			}
			if ContainsChild(p.loops[i], p.loops[j], loopMap) != p.loops[i].ContainsNested(p.loops[j]) {
				log.Println("FAILED EQ TEST")
			}
		}
	}

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

func NewPolygonFromLoops(loops *[]*Loop) *Polygon {
	p := &Polygon{
		bound:     EmptyRect(),
		ownsLoops: true,
		hasHoles:  false,
	}
	p.Init(loops)
	return p
}

func NewPolygonFromLoop(loop *Loop) *Polygon {
	p := &Polygon{
		loops:       []*Loop{},
		bound:       loop.Bound(),
		ownsLoops:   false,
		hasHoles:    false,
		numVertices: len(loop.vertices),
	}
	p.loops = append(p.loops, loop)
	return p
}

func (p Polygon) CapBound() Cap {
	return p.bound.CapBound()
}

func (p Polygon) RectBound() Rect {
	return p.bound
}

func (a Polygon) ContainsOrCrosses(b *Loop) int {
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

func (a Polygon) AnyLoopContains(b *Loop) bool {
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

func (a Polygon) ContainsPoint(p Point) bool {
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

func (p Polygon) ContainsCell(cell Cell) bool {
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

func (p Polygon) MayIntersect(cell Cell) bool {
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
	// No loop may contain the complement of another loop. (Handling this
	// case is significantly more complicated).
	//
	// Some of the children of the parent loop may now be children of the
	// new loop.
	for i := 0; i < len(loopMap[parent]); {
		child := loopMap[parent][i]
		if newLoop.ContainsNested(child) {
			loopMap[newLoop] = append(loopMap[newLoop], child)
			//			copy(loopMap[parent][i:], loopMap[parent][i+1:])
			//			loopMap[parent][len(loopMap[parent])-1] = nil
			//			loopMap[parent] = loopMap[parent][:len(loopMap[parent])-1]
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
		edges := map[PointPair]IntPair{}
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
	vertices := map[Point]bool{}
	for i := 0; i < p.NumLoops(); i++ {
		child := p.Loop(i)
		if child.depth == 0 {
			continue
		}
		parent := p.Loop(p.Parent(i))
		if parent != lastParent {
			vertices := map[Point]bool{}
			for j := 0; j < parent.NumVertices(); j++ {
				vertices[*parent.vertex(j)] = true
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

	options := DIRECTED_XOR()
	options.vertex_merge_radius = vertexMergeRadius
	builder := NewPolygonBuilder(options)
	ClipBoundary(a, false, b, false, false, true, &builder)
	ClipBoundary(b, false, a, false, false, false, &builder)
	if !builder.AssemblePolygon(p, nil) {
		log.Fatalf("Bad directed edges in InitToIntersection")
	}
}

type FloatPointPair struct {
	first  float64
	second Point
}

type byFloatPointPair []FloatPointPair

func (x byFloatPointPair) Len() int      { return len(x) }
func (x byFloatPointPair) Swap(i, j int) { x[i], x[j] = x[j], x[i] }
func (x byFloatPointPair) Less(i, j int) bool {
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

func AddIntersection(a0, a1, b0, b1 Point, add_shared_edges bool, crossing int,
	intersections *[]FloatPointPair) {
	if crossing > 0 {
		// There is a proper edge crossing.
		x := GetIntersection(a0, a1, b0, b1)
		t := GetDistanceFraction(x, a0, a1)
		*intersections = append(*intersections, FloatPointPair{t, x})
	}
}

func ClipEdge(a0, a1 Point, bIndex *PolygonIndex, add_shared_edges bool,
	intersections *[]FloatPointPair) {
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

func ClipBoundary(a *Polygon, reverse_a bool,
	b *Polygon, reverse_b bool,
	invert_b, add_shared_edges bool,
	builder *PolygonBuilder) {
	bIndex := NewPolygonIndex(b, reverse_b)

	intersections := []FloatPointPair{}
	for _, la := range a.loops {
		n := la.NumVertices()
		dir := 1
		if la.IsHole() != reverse_a {
			dir = -1
		}
		inside := b.ContainsPoint(*la.vertex(0)) != invert_b
		j := n
		if dir > 0 {
			j = 0
		}
		for ; n > 0; n, j = n-1, j+dir {
			a0 := *la.vertex(j)
			a1 := *la.vertex(j + dir)
			intersections = []FloatPointPair{}
			ClipEdge(a0, a1, &bIndex, add_shared_edges, &intersections)
			if inside {
				intersections = append(intersections, FloatPointPair{0, a0})
			}
			inside = (len(intersections) & 1) != 0
			sort.Sort(byFloatPointPair(intersections))
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
	EdgeIndex
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

// TODO: interface EdgeIndex because it's used everywhere.
type PolygonIndex struct {
	LoopSequenceIndex
	poly    *Polygon
	reverse bool
}

func (idx PolygonIndex) EdgeFromTo(i int) (from *Point, to *Point) {
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

func (idx PolygonIndex) edge_from(i int) *Point {
	from, _ := idx.EdgeFromTo(i)
	return from
}
func (idx PolygonIndex) edge_to(i int) *Point {
	_, to := idx.EdgeFromTo(i)
	return to
}

func NewPolygonIndex(poly *Polygon, reverse bool) PolygonIndex {
	p := PolygonIndex{LoopSequenceIndex: NewLoopSequenceIndex(), poly: poly, reverse: reverse}
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
