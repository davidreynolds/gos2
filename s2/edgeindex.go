package s2

import (
	"math"
	"sort"
)

var (
	alwaysRecurseOnChildren bool
)

func init() {
	alwaysRecurseOnChildren = false
}

type CellEdge struct {
	cellId CellID
	edgeId int
}

type CellEdgeMultimap struct {
	items []CellEdge
}

func NewCellEdgeMultimap() CellEdgeMultimap { return CellEdgeMultimap{[]CellEdge{}} }

func (m *CellEdgeMultimap) Len() int      { return len(m.items) }
func (m *CellEdgeMultimap) Swap(i, j int) { m.items[i], m.items[j] = m.items[j], m.items[i] }
func (m *CellEdgeMultimap) Less(i, j int) bool {
	if uint64(m.items[i].cellId) < uint64(m.items[j].cellId) {
		return true
	}
	if uint64(m.items[j].cellId) < uint64(m.items[i].cellId) {
		return false
	}
	if m.items[i].edgeId < m.items[j].edgeId {
		return true
	}
	return false
}

func (m *CellEdgeMultimap) Insert(ce CellEdge) {
	m.items = append(m.items, ce)
	sort.Sort(m)
}

func (m CellEdgeMultimap) LowerBound(id CellID) int {
	idx := sort.Search(m.Len(), func(k int) bool {
		return uint64(m.items[k].cellId) >= uint64(id)
	})
	return idx
}

func (m CellEdgeMultimap) UpperBound(id CellID) int {
	// XXX: Is this a good way to get the upper bound?
	return m.LowerBound(id.Next())
}

type EdgeIndexer interface {
	QueryCount() int
	SetQueryCount(n int)
	IncrementQueryCount()
	MinLevelUsed() int
	SetMinLevelUsed(level int)
	NumEdges() int
	IndexComputed() bool
	SetIndexComputed(b bool)
	Reset()
	MappingInsert(id CellID, k int)
	Mapping() *CellEdgeMultimap
	edge_from(i int) *Point
	edge_to(i int) *Point
}

type EdgeIndex struct {
	mapping       CellEdgeMultimap
	minLevelUsed  int
	queryCount    int
	indexComputed bool
}

func NewEdgeIndex() EdgeIndex {
	return EdgeIndex{
		indexComputed: false,
		queryCount:    0,
		minLevelUsed:  maxLevel,
		mapping:       NewCellEdgeMultimap(),
	}
}

func (idx *EdgeIndex) Reset() {
	idx.SetMinLevelUsed(maxLevel)
	idx.SetIndexComputed(false)
	idx.SetQueryCount(0)
	idx.mapping = NewCellEdgeMultimap()
}

func (idx EdgeIndex) QueryCount() int                 { return idx.queryCount }
func (idx *EdgeIndex) SetQueryCount(n int)            { idx.queryCount = n }
func (idx *EdgeIndex) IncrementQueryCount()           { idx.queryCount++ }
func (idx EdgeIndex) MinLevelUsed() int               { return idx.minLevelUsed }
func (idx *EdgeIndex) SetMinLevelUsed(level int)      { idx.minLevelUsed = level }
func (idx EdgeIndex) IndexComputed() bool             { return idx.indexComputed }
func (idx *EdgeIndex) SetIndexComputed(b bool)        { idx.indexComputed = b }
func (idx *EdgeIndex) MappingInsert(id CellID, k int) { idx.mapping.Insert(CellEdge{id, k}) }
func (idx *EdgeIndex) Mapping() *CellEdgeMultimap     { return &idx.mapping }

// Appends to "candidate_crossings" all edge references which may cross the
// given edge. This is done by covering the edge and then finding all references
// of edges whose coverings overlap this covering. Parent cells are checked
// level by level. Child cells are checked all at once by taking advantage of
// the natural ordering of CellIDs.
func FindCandidateCrossings(idx EdgeIndexer, a, b Point, candidate_crossings *[]int) {
	cover, _ := GetCovering(idx, a, b, false)
	EdgesInParentCells(idx, cover, idx.Mapping(), idx.MinLevelUsed(), candidate_crossings)
	EdgesInChildrenCells(idx, a, b, &cover, idx.Mapping(), candidate_crossings)
	uniq := map[int]bool{}
	for _, c := range *candidate_crossings {
		uniq[c] = true
	}
	*candidate_crossings = []int{}
	for k, _ := range uniq {
		*candidate_crossings = append(*candidate_crossings, k)
	}
}

func EdgesInParentCells(idx EdgeIndexer, cover []CellID, mapping *CellEdgeMultimap, min_level_used int, candidate_crossings *[]int) {
	// Find all parent cells of covering cells.
	parent_cells := map[CellID]bool{}
	for _, cid := range cover {
		for parentLevel := cid.Level() - 1; parentLevel >= min_level_used; parentLevel-- {
			if _, ok := parent_cells[cid.Parent(parentLevel)]; !ok {
				parent_cells[cid.Parent(parentLevel)] = true
			} else {
				break
			}
		}
	}

	// Put parent cell edge references into result.
	for pi, _ := range parent_cells {
		start := mapping.LowerBound(pi)
		for i := start; i != mapping.Len() && mapping.items[i].cellId == pi; i++ {
			*candidate_crossings = append(*candidate_crossings, mapping.items[i].edgeId)
		}
	}
}

func EdgesInChildrenCells(idx EdgeIndexer, a, b Point, cover *[]CellID, mapping *CellEdgeMultimap, candidate_crossings *[]int) {
	num_cells := 0

	// Put all the edge references of (covering cells + descendant cells)
	// into result. This relies on the natural ordering of CellIDs.
	for len(*cover) > 0 {
		back := len(*cover) - 1
		cell := (*cover)[back]
		*cover = (*cover)[:back]
		num_cells++

		start := mapping.LowerBound(cell.RangeMin())
		end := mapping.UpperBound(cell.RangeMax())

		rewind := alwaysRecurseOnChildren
		num_edges := 0
		if !rewind {
			for i := start; i < mapping.Len() && i != end; i++ {
				*candidate_crossings = append(*candidate_crossings, mapping.items[i].edgeId)
				num_edges++
				if num_edges == 16 && !cell.IsLeaf() {
					rewind = true
					break
				}
			}
		}

		// If there were too many to insert, uninsert and recurse.
		if rewind {
			for i := 0; i < num_edges; i++ {
				back := len(*candidate_crossings) - 1
				*candidate_crossings = (*candidate_crossings)[:back]
			}

			i := mapping.LowerBound(cell)
			j := i
			for ; j != mapping.Len() && mapping.items[j].cellId == cell; j++ {
				*candidate_crossings = append(*candidate_crossings, mapping.items[j].edgeId)
			}

			// Recurse on the children -- hopefully some will be empty.
			if i != start || j != end {
				children := [4]Cell{}
				c := CellFromCellID(cell)
				c.Subdivide(&children)
				for _, child := range children {
					if EdgeIntersectsCellBoundary(a, b, child) {
						*cover = append(*cover, child.Id())
					}
				}
			}
		}
	}
}

func ComputeIndex(idx EdgeIndexer) {
	for i := 0; i < idx.NumEdges(); i++ {
		cover, level := GetCovering(idx, *idx.edge_from(i), *idx.edge_to(i), true)
		idx.SetMinLevelUsed(min(idx.MinLevelUsed(), level))
		for _, cid := range cover {
			idx.MappingInsert(cid, i)
		}
	}
	idx.SetIndexComputed(true)
}

// Returns the smallest cell containing all four points, or Sentinel if they
// are not all on the same face. The points don't need to be normalized.
func containingCell4(pa, pb, pc, pd Point) CellID {
	a := cellIDFromPoint(pa)
	b := cellIDFromPoint(pb)
	c := cellIDFromPoint(pc)
	d := cellIDFromPoint(pd)
	if a.Face() != b.Face() || a.Face() != c.Face() || a.Face() != d.Face() {
		return Sentinel()
	}

	for a != b || a != c || a != d {
		a = a.immediateParent()
		b = b.immediateParent()
		c = c.immediateParent()
		d = d.immediateParent()
	}
	return a
}

// Returns the smallest cell containing both points, or Sentinel if they
// are not all on the same face. The points don't need to be normalized.
func containingCell2(pa, pb Point) CellID {
	a := cellIDFromPoint(pa)
	b := cellIDFromPoint(pb)
	if a.Face() != b.Face() {
		return Sentinel()
	}
	for a != b {
		a = a.immediateParent()
		b = b.immediateParent()
	}
	return a
}

func GetCovering(idx EdgeIndexer, a, b Point, thicken_edge bool) ([]CellID, int) {
	covering := []CellID{}
	// Thicken the edge in all directions by roughly 1% of the edge length
	// when thicken_edge is true.
	thickening := 0.01

	// Selects the ideal S2 level at which to cover the edge, this will be
	// the level whose S2 cells have a width roughly commensurate to the
	// length of the edge. We multiply the edge length by 2*thickening to
	// guarantee the thickening is honored when doing the covering-by-cap
	// trick.
	edge_length := a.Angle(b.Vector).Radians()
	ideal_level := MinWidth.MaxLevel(edge_length * (1 + 2*thickening))
	var containing_cell CellID
	if !thicken_edge {
		containing_cell = containingCell2(a, b)
	} else {
		if ideal_level == maxLevel {
			// If the edge is tiny, instabilities are more likely,
			// so we want to limit the number of operations.
			// We pretend we are in a cell much larger so as to
			// trigger the 'needs covering' case, so we won't try to
			// thicken the edge.
			containing_cell = CellID(0xFFF0).Parent(3)
		} else {
			pq := b.Sub(a.Vector).Mul(thickening)
			ortho := pq.Cross(a.Vector).Normalize().Mul(edge_length * thickening)
			p := a.Sub(pq)
			q := b.Add(pq)
			containing_cell = containingCell4(Point{p.Sub(ortho)}, Point{p.Add(ortho)},
				Point{q.Sub(ortho)}, Point{q.Add(ortho)})
		}
	}

	// Best case: edge is fully contained in a cell that's not too big.
	if containing_cell != Sentinel() && containing_cell.Level() >= ideal_level-2 {
		covering = append(covering, containing_cell)
		return covering, containing_cell.Level()
	}

	if ideal_level == 0 {
		// Edge is very long, maybe even longer than a face width, so
		// the trick below doesn't work. For now, we will add the whole
		// S2 sphere.
		for cid := CellIDBegin(0); cid != CellIDEnd(0); cid = cid.Next() {
			covering = append(covering, cid)
		}
		return covering, 0
	}

	// Cover the edge by a cap centered on the edge midpoint, then cover
	// the cap by four big-enough cells around the cell vertex closest to
	// the cap center.
	middle := Point{a.Add(b.Vector).Mul(0.5).Normalize()}
	actual_level := min(ideal_level, maxLevel-1)
	cellIDFromPoint(middle).AppendVertexNeighbors(actual_level, &covering)
	return covering, actual_level
}

func PredictAdditionalCalls(idx EdgeIndexer, n int) {
	if !idx.IndexComputed() {
		if idx.NumEdges() > 100 && (idx.QueryCount()+n) > 30 {
			ComputeIndex(idx)
		}
	}
}

type EdgeIndexerIterator interface {
	GetCandidates(a, b Point)
	Index() int
	Next()
	Done() bool
}

// An iterator on loops/edges that may cross a query edge (a,b).
type EdgeIndexIterator struct {
	current                     int
	num_edges                   int
	is_brute_force              bool
	index                       EdgeIndexer
	candidates                  []int
	current_index_in_candidates int
}

func NewEdgeIndexIterator(idx EdgeIndexer) *EdgeIndexIterator {
	return &EdgeIndexIterator{current: 0, is_brute_force: false, index: idx}
}

func (it *EdgeIndexIterator) GetCandidates(a, b Point) {
	PredictAdditionalCalls(it.index, 1)
	it.is_brute_force = !it.index.IndexComputed()
	if it.is_brute_force {
		it.index.IncrementQueryCount()
		it.current = 0
		it.num_edges = it.index.NumEdges()
	} else {
		it.candidates = []int{}
		FindCandidateCrossings(it.index, a, b, &it.candidates)
		it.current_index_in_candidates = 0
		if len(it.candidates) != 0 {
			it.current = it.candidates[0]
		}
	}
}

// Index of the current loop in the iteration.
func (it EdgeIndexIterator) Index() int { return it.current }

// Iterate to the next available candidate.
func (it *EdgeIndexIterator) Next() {
	if it.is_brute_force {
		it.current++
	} else {
		it.current_index_in_candidates++
		if it.current_index_in_candidates < len(it.candidates) {
			it.current = it.candidates[it.current_index_in_candidates]
		}
	}
}

// True if there are no more candidates.
func (it EdgeIndexIterator) Done() bool {
	if it.is_brute_force {
		return it.current >= it.num_edges
	} else {
		return it.current_index_in_candidates >= len(it.candidates)
	}
}

func LenientCrossing(a, b, c, d Point) bool {
	maxDetError := 1e-14
	acb := a.Cross(c.Vector).Dot(b.Vector)
	bda := b.Cross(d.Vector).Dot(a.Vector)
	if math.Abs(acb) < maxDetError || math.Abs(bda) < maxDetError {
		return true
	}
	if acb*bda < 0 {
		return false
	}
	cbd := c.Cross(b.Vector).Dot(d.Vector)
	dac := d.Cross(a.Vector).Dot(c.Vector)
	if math.Abs(cbd) < maxDetError || math.Abs(dac) < maxDetError {
		return true
	}
	return (acb*cbd >= 0) && (acb*dac >= 0)
}

func EdgeIntersectsCellBoundary(a, b Point, cell Cell) bool {
	vertices := make([]Point, 4)
	for i := 0; i < 4; i++ {
		vertices[i] = cell.Vertex(i)
	}
	for i := 0; i < 4; i++ {
		from := vertices[i]
		to := vertices[(i+1)%4]
		if LenientCrossing(a, b, from, to) {
			return true
		}
	}
	return false
}
