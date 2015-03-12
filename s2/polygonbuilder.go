package s2

import (
	"log"
	"sort"

	"github.com/davidreynolds/gos2/s1"
)

type Edge struct {
	v0, v1 Point
}

type PolygonBuilderOptions struct {
	xor_edges            bool
	undirected_edges     bool
	validate             bool
	vertex_merge_radius  s1.Angle
	edge_splice_fraction float64
}

func (op *PolygonBuilderOptions) SetUndirectedEdges(b bool)           { op.undirected_edges = b }
func (op *PolygonBuilderOptions) SetXorEdges(b bool)                  { op.xor_edges = b }
func (op *PolygonBuilderOptions) SetValidate(b bool)                  { op.validate = b }
func (op *PolygonBuilderOptions) SetVertexMergeRadius(angle s1.Angle) { op.vertex_merge_radius = angle }

// These are the options that should be used for assembling well-behaved
// input data into polygons. All edges should be directed such that
// "shells" and "holes" have opposite orientations (typically CCW shells and
// clockwise holes), unless it is known that shells and holes do not share
// any edges.
func DIRECTED_XOR() PolygonBuilderOptions {
	return PolygonBuilderOptions{
		xor_edges:            true,
		undirected_edges:     false,
		validate:             false,
		vertex_merge_radius:  0,
		edge_splice_fraction: 0.866,
	}
}

// These are the options that should be used for assembling polygons that do
// not follow the conventions above, e.g. where edge directions may vary
// within a single loop, or shells and holes are not oppositely oriented and
// may share edges.
func UNDIRECTED_XOR() PolygonBuilderOptions {
	return PolygonBuilderOptions{
		xor_edges:            true,
		undirected_edges:     true,
		validate:             false,
		vertex_merge_radius:  0,
		edge_splice_fraction: 0.866,
	}
}

// These are the options that should be used for assembling edges where the
// desired output is a collection of loops rather than a polygon, and edges
// may occur more than once. Edges are treated as undirected and are not
// XORed together.
func UNDIRECTED_UNION() PolygonBuilderOptions {
	return PolygonBuilderOptions{
		xor_edges:            false,
		undirected_edges:     true,
		validate:             false,
		vertex_merge_radius:  0,
		edge_splice_fraction: 0.866,
	}
}

type CellIDPoint struct {
	id    CellID
	point Point
}

type byIDPoint []CellIDPoint

func (x byIDPoint) Len() int           { return len(x) }
func (x byIDPoint) Less(i, j int) bool { return uint64(x[i].id) < uint64(x[j].id) }
func (x byIDPoint) Swap(i, j int)      { x[i], x[j] = x[j], x[i] }

type PointIndex struct {
	map_          []CellIDPoint
	vertex_radius float64
	edge_fraction float64
	level         int
}

func NewPointIndex(vertex_radius, edge_fraction float64) PointIndex {
	max_level := MinWidth.MaxLevel(2 * vertex_radius)
	idx := PointIndex{
		vertex_radius: vertex_radius,
		edge_fraction: edge_fraction,
		level:         min(max_level, maxLevel-1),
	}
	idx.map_ = []CellIDPoint{CellIDPoint{Sentinel(), Point{}}}
	return idx
}

func (idx PointIndex) FindNearbyPoint(v0, v1 Point, out *Point) bool {
	// Return a point whose distance from the edge (v0, v1) is less than
	// vertex_radius, and which is not equal to v0 or v1. The current
	// implementation returns the closest such point.
	//
	// Strategy: we compute a very cheap covering by approximating the edge
	// as two spherical caps, one around each endpoint, and then computing a
	// 4-cell covering using each one. We could improve the quality of the
	// covering by using some intermediate points along the edge as well.
	length := v0.Distance(v1).Radians()
	normal := v0.PointCross(v1)
	level := min(idx.level, MinWidth.MaxLevel(length))
	ids := []CellID{}
	cellIDFromPoint(v0).AppendVertexNeighbors(level, &ids)
	cellIDFromPoint(v1).AppendVertexNeighbors(level, &ids)

	// Sort the cell ids so that we can skip duplicates in the loop below.
	sort.Sort(byID(ids))

	bestDist := 2 * idx.vertex_radius
	for i := len(ids) - 1; i >= 0; i-- {
		// Skip duplicates
		if i > 0 && ids[i-1] == ids[i] {
			continue
		}

		maxId := ids[i].RangeMax()
		j := sort.Search(len(idx.map_), func(k int) bool {
			return uint64(idx.map_[k].id) >= uint64(ids[i].RangeMin())
		})
		for ; idx.map_[j].id <= maxId; j++ {
			p := idx.map_[j].point
			if p == v0 || p == v1 {
				continue
			}
			dist := p.DistanceToEdgeWithNormal(v0, v1, normal).Radians()
			if dist < bestDist {
				bestDist = dist
				*out = p
			}
		}
	}
	return bestDist < idx.edge_fraction*idx.vertex_radius
}

func (idx *PointIndex) Insert(p Point) {
	ids := []CellID{}
	cellIDFromPoint(p).AppendVertexNeighbors(idx.level, &ids)
	for i := len(ids) - 1; i >= 0; i-- {
		idx.map_ = append(idx.map_, CellIDPoint{ids[i], p})
	}
	sort.Sort(byIDPoint(idx.map_))
}

func (idx *PointIndex) Erase(p Point) {
	ids := []CellID{}
	cellIDFromPoint(p).AppendVertexNeighbors(idx.level, &ids)
	for i := len(ids) - 1; i >= 0; i-- {
		j := sort.Search(len(idx.map_), func(k int) bool {
			return uint64(idx.map_[k].id) >= uint64(ids[i])
		})
		for ; idx.map_[j].point != p; j++ {
			if ids[i] != idx.map_[j].id {
				log.Println("don't match:", ids[i], idx.map_[j])
			}
		}
		idx.map_ = append(idx.map_[:j], idx.map_[j+1:]...)
	}
}

func (idx PointIndex) QueryCap(axis Point) []Point {
	out := []Point{}
	id := cellIDFromPoint(axis).Parent(idx.level)
	i := sort.Search(len(idx.map_), func(k int) bool {
		return uint64(idx.map_[k].id) >= uint64(id)
	})
	for ; idx.map_[i].id == id; i++ {
		p := idx.map_[i].point
		dist := axis.Distance(p).Radians()
		if dist < idx.vertex_radius {
			out = append(out, p)
		}
	}
	return out
}

type VertexSet struct {
	vertices []Point
}

func NewVertexSet() *VertexSet {
	return &VertexSet{
		vertices: []Point{},
	}
}

func (vs *VertexSet) Len() int           { return len(vs.vertices) }
func (vs *VertexSet) Swap(i, j int)      { vs.vertices[i], vs.vertices[j] = vs.vertices[j], vs.vertices[i] }
func (vs *VertexSet) Less(i, j int) bool { return vs.vertices[i].LessThan(vs.vertices[j].Vector) }

func (vs *VertexSet) Insert(p Point) {
	vs.vertices = append(vs.vertices, p)
	sort.Sort(vs)
}

func (vs VertexSet) Find(p Point) int {
	idx := sort.Search(vs.Len(), func(k int) bool {
		return vs.vertices[k].GTE(p.Vector)
	})
	if idx < vs.Len() && vs.vertices[idx] == p {
		return idx
	}
	return vs.Len()
}

func (vs *VertexSet) Erase(i int) {
	if i >= 0 && i < vs.Len() {
		vs.vertices = append(vs.vertices[:i], vs.vertices[i+1:]...)
	}
}

type EdgeSet map[Point]*VertexSet

type PolygonBuilder struct {
	options             PolygonBuilderOptions
	edges               EdgeSet
	starting_vertices   []Point
	vertex_merge_radius s1.Angle
}

func NewPolygonBuilder(options PolygonBuilderOptions) PolygonBuilder {
	return PolygonBuilder{
		options:           options,
		edges:             EdgeSet{},
		starting_vertices: []Point{},
	}
}

func (b PolygonBuilder) HasEdge(v0, v1 Point) bool {
	if vset, ok := b.edges[v0]; ok {
		return vset.Find(v1) != vset.Len()
	}
	return false
}

func (b *PolygonBuilder) AddEdge(v0, v1 Point) bool {
	if v0 == v1 {
		return false
	}
	if b.options.xor_edges && b.HasEdge(v1, v0) {
		b.EraseEdge(v1, v0)
		return false
	}
	if _, ok := b.edges[v0]; !ok {
		b.edges[v0] = NewVertexSet()
		b.starting_vertices = append(b.starting_vertices, v0)
	}
	b.edges[v0].Insert(v1)
	if b.options.undirected_edges {
		if _, ok := b.edges[v1]; !ok {
			b.edges[v1] = NewVertexSet()
			b.starting_vertices = append(b.starting_vertices, v1)
		}
		b.edges[v1].Insert(v0)
	}
	return true
}

func (b *PolygonBuilder) AddLoop(loop *Loop) {
	sign := loop.Sign()
	for i := len(loop.vertices); i > 0; i-- {
		// Vertex indices need to be in the range [0, 2*num_vertices-1].
		b.AddEdge(*loop.vertex(i), *loop.vertex(i + sign))
	}
}

func (b *PolygonBuilder) AddPolygon(polygon *Polygon) {
	for _, loop := range polygon.loops {
		b.AddLoop(loop)
	}
}

func (b *PolygonBuilder) EraseEdge(v0, v1 Point) {
	if vset, ok := b.edges[v0]; ok {
		vset.Erase(vset.Find(v1))
		if vset.Len() == 0 {
			delete(b.edges, v0)
		}
	}

	if b.options.undirected_edges {
		if vset, ok := b.edges[v1]; ok {
			vset.Erase(vset.Find(v0))
			if vset.Len() == 0 {
				delete(b.edges, v1)
			}
		}
	}
}

type MergeMap map[Point]Point

// The overall strategy is to start from each vertex and grow a maximal
// cluster of mergeable vertices. In graph theoretic terms, we find the
// connected components of the undirected graph whose edges connect pairs of
// vertices that are separated by at most vertex_merge_radius.
//
// We then choose a single representative vertex for each cluster, and
// update "merge_map" appropriately. We choose an arbitrary existing vertex
// rather than computing the centroid of all the vertices to avoid creating new
// vertex pairs that need to be merged. (We guarantee that all vertex pairs
// are separated by at least the merge radius in the output.)
func (b *PolygonBuilder) BuildMergeMap(index *PointIndex) MergeMap {
	// First, we build the set of all the distinct vertices in the input.
	// We need to include the source and destination of every edge.
	vertices := make(map[Point]bool)
	merge_map := MergeMap{}

	for v0, vset := range b.edges {
		vertices[v0] = true
		for _, v1 := range vset.vertices {
			vertices[v1] = true
		}
	}

	// Build a spatial index containing all the distinct vertices
	for p := range vertices {
		index.Insert(p)
	}

	// Next we loop through all the vertices and attempt to grow a maximal
	// mergeable group starting from each vertex
	frontier := []Point{}
	for p := range vertices {
		// Skip any vertices that have already been merged with
		// another vertex
		if _, ok := merge_map[p]; ok {
			continue
		}

		// Grow a maximal mergeable component starting from "p", the
		// canonical representation of the mergeable group
		frontier = append(frontier, p)
		for len(frontier) > 0 {
			i := len(frontier) - 1
			mergeable := index.QueryCap(frontier[i])
			frontier = frontier[:i]
			for j := len(mergeable) - 1; j >= 0; j-- {
				v1 := mergeable[j]
				if v1 != p {
					// Erase from the index any vertices that will be merged.
					// This ensures that we won't try to merge the same vertex twice.
					index.Erase(v1)
					frontier = append(frontier, v1)
					merge_map[v1] = p
				}
			}
		}
	}
	return merge_map
}

func (b *PolygonBuilder) MoveVertices(merge_map MergeMap) {
	if len(merge_map) == 0 {
		return
	}
	edges := []Edge{}

	for v0, vset := range b.edges {
		for _, v1 := range vset.vertices {
			_, ok0 := merge_map[v0]
			_, ok1 := merge_map[v1]
			if ok0 || ok1 {
				if !b.options.undirected_edges || v0.LessThan(v1.Vector) {
					edges = append(edges, Edge{v0, v1})
				}
			}
		}
	}

	// Now erase all the old edges and add all the new edges. This will
	// automatically take care of any XORing that needs to be done, because
	// EraseEdge also erases the sibling of undirected edges.
	for _, e := range edges {
		v0 := e.v0
		v1 := e.v1
		b.EraseEdge(v0, v1)
		if new0, ok := merge_map[v0]; ok {
			v0 = new0
		}
		if new1, ok := merge_map[v1]; ok {
			v1 = new1
		}
		b.AddEdge(v0, v1)
	}
}

func (b *PolygonBuilder) SpliceEdges(index *PointIndex) {
	// We keep a stack of unprocessed edges. Initially all edges are
	// pushed onto the stack.
	edges := []Edge{}
	for v0, vset := range b.edges {
		for _, v1 := range vset.vertices {
			if !b.options.undirected_edges || v0.LessThan(v1.Vector) {
				edges = append(edges, Edge{v0, v1})
			}
		}
	}

	// For each edge, we check whether there are any nearby vertices that
	// should be spliced into it. If there are, we choose one such vertex,
	// split the edge into two pieces, and iterate on each piece.
	for len(edges) > 0 {
		i := len(edges) - 1
		v0 := edges[i].v0
		v1 := edges[i].v1
		edges = edges[:i]

		// If we are XORing, edges may be erased before we get to them.
		if b.options.xor_edges && !b.HasEdge(v0, v1) {
			continue
		}

		var vmid Point
		if !index.FindNearbyPoint(v0, v1, &vmid) {
			continue
		}

		b.EraseEdge(v0, v1)
		if b.AddEdge(v0, vmid) {
			edges = append(edges, Edge{v0, vmid})
		}
		if b.AddEdge(vmid, v1) {
			edges = append(edges, Edge{vmid, v1})
		}
	}
}

func (b *PolygonBuilder) EraseLoop(loop *Loop) {
	n := len(loop.vertices)
	for i, j := n-1, 0; j < n; i, j = j, j+1 {
		b.EraseEdge(loop.vertices[i], loop.vertices[j])
	}
}

func (b *PolygonBuilder) AssembleLoop(v0, v1 Point, unused_edges *[]Edge) *Loop {
	path := []Point{}            // The path so far.
	index := make(map[Point]int) // Maps a vertex to its index in "path"
	path = append(path, v0)
	path = append(path, v1)
	index[v1] = 1

	for len(path) >= 2 {
		v0 := path[len(path)-2]
		v1 := path[len(path)-1]
		var v2 Point
		v2_found := false

		if vset, ok := b.edges[v1]; ok {
			for _, v := range vset.vertices {
				if v == v0 {
					continue
				}
				if !v2_found || OrderedCCW(v0, v2, v, v1) {
					v2 = v
				}
				v2_found = true
			}
		}

		if !v2_found {
			// We've hit a dead end. Remove this edge and backtrack.
			*unused_edges = append(*unused_edges, Edge{v0, v1})
			b.EraseEdge(v0, v1)
			delete(index, v1)
			path = path[:len(path)-1]
		} else {
			if _, ok := index[v2]; !ok {
				// This is the first time we've visited this
				// vertex.
				index[v2] = len(path)
				path = append(path, v2)
			} else {
				// We've completed a loop. Throw away any
				// initial vertices that are not part of the
				// loop.
				path = path[index[v2]:]
				loop := NewLoopFromPath(path)
				if b.options.validate && !loop.IsValid() {
					b.RejectLoop(loop, unused_edges)
					b.EraseLoop(loop)
					return nil
				}
				if b.options.undirected_edges && !loop.IsNormalized() {
					return b.AssembleLoop(path[1], path[0], unused_edges)
				}
				return loop
			}
		}
	}
	return nil
}

func (b *PolygonBuilder) AssembleLoops(loops *[]*Loop, unused_edges *[]Edge) bool {
	if b.options.vertex_merge_radius.Radians() > 0 {
		index := NewPointIndex(b.options.vertex_merge_radius.Radians(),
			b.options.edge_splice_fraction)
		mergemap := b.BuildMergeMap(&index)
		b.MoveVertices(mergemap)
		if b.options.edge_splice_fraction > 0 {
			b.SpliceEdges(&index)
		}
	}

	var dummy_unused_edges []Edge
	if unused_edges == nil {
		unused_edges = &dummy_unused_edges
	}
	*unused_edges = []Edge{}

	// We repeatedly choose an edge and attempt to assemble a loop
	// starting from that edge. (This is always possible unless the
	// input includes extra edges that are not part of any loop.) To
	// maintain a consistent scanning order over b.edges between
	// different machine architectures (e.g. 'clovertown' vs 'opteron'),
	// we follow the order they were added to the builder.

	for i := 0; i < len(b.starting_vertices); {
		v0 := b.starting_vertices[i]
		if candidates, ok := b.edges[v0]; ok {
			v1 := candidates.vertices[0]
			loop := b.AssembleLoop(v0, v1, unused_edges)
			if loop != nil {
				*loops = append(*loops, loop)
				b.EraseLoop(loop)
			}
		} else {
			i++
		}
	}
	return len(*unused_edges) == 0
}

func (b *PolygonBuilder) AssemblePolygon(polygon *Polygon, unusedEdges *[]Edge) bool {
	var loops []*Loop
	success := b.AssembleLoops(&loops, unusedEdges)
	// If edges are undirected, then all loops are already CCW. Otherwise
	// we need to make sure the loops are normalized.
	if !b.options.undirected_edges {
		for _, loop := range loops {
			loop.Normalize()
		}
	}
	if b.options.validate && !AreLoopsValid(loops) {
		if unusedEdges != nil {
			for _, loop := range loops {
				b.RejectLoop(loop, unusedEdges)
			}
		}
		return false
	}
	polygon.Init(loops)
	return success
}

func (b *PolygonBuilder) RejectLoop(loop *Loop, unusedEdges *[]Edge) {
	n := len(loop.vertices)
	for i, j := n-1, 0; j < n; i, j = j, j+1 {
		*unusedEdges = append(*unusedEdges, Edge{loop.vertices[i], loop.vertices[j]})
	}
}
