package s2

import "container/heap"

const (
	// By default, the covering uses at most 8 cells at any level. This
	// gives a reasonable tradeoff between the number of cells used and the
	// accuracy of the approximation.
	defaultMaxCells = 8
)

type Candidate struct {
	cell       Cell
	isTerminal bool // Cell should not be expanded further
	children   []*Candidate
}

type CandidateEntry struct {
	candidate *Candidate
	priority  int
}

type PriorityQueue []*CandidateEntry

func (pq PriorityQueue) Len() int            { return len(pq) }
func (pq PriorityQueue) Less(i, j int) bool  { return pq[i].priority > pq[j].priority }
func (pq PriorityQueue) Swap(i, j int)       { pq[i], pq[j] = pq[j], pq[i] }
func (pq *PriorityQueue) Push(x interface{}) { *pq = append(*pq, x.(*CandidateEntry)) }

func (pq *PriorityQueue) Pop() interface{} {
	old := *pq
	n := len(*pq)
	entry := old[n-1]
	*pq = old[0 : n-1]
	return entry
}

type RegionCoverer struct {
	minLevel                 int
	maxLevel                 int
	levelMod                 int
	maxCells                 int
	result                   []CellID
	interiorCovering         bool
	candidatesCreatedCounter int
	region                   *Region
	pq                       *PriorityQueue
}

func NewRegionCoverer() RegionCoverer {
	r := RegionCoverer{
		minLevel: 0,
		maxLevel: maxLevel,
		levelMod: 1,
		maxCells: defaultMaxCells,
		result:   []CellID{},
		pq:       &PriorityQueue{},
	}
	heap.Init(r.pq)
	return r
}

// Returns the log base 2 of the maximum number of children of a candidate.
func (r RegionCoverer) maxChildrenShift() int  { return 2 * r.levelMod }
func (r *RegionCoverer) SetMinLevel(level int) { r.minLevel = max(0, min(maxLevel, level)) }
func (r *RegionCoverer) SetMaxLevel(level int) { r.maxLevel = max(0, min(maxLevel, level)) }
func (r *RegionCoverer) SetLevelMod(mod int)   { r.levelMod = max(1, min(3, mod)) }
func (r *RegionCoverer) SetMaxCells(cells int) { r.maxCells = cells }

func (r *RegionCoverer) NewCandidate(cell Cell) *Candidate {
	if !(*r.region).MayIntersect(cell) {
		return nil
	}
	isTerminal := false
	cellLevel := int(cell.level)
	if cellLevel >= r.minLevel {
		if r.interiorCovering {
			if (*r.region).ContainsCell(cell) {
				isTerminal = true
			} else if cellLevel+r.levelMod > r.maxLevel {
				return nil
			}
		} else {
			if cellLevel+r.levelMod > r.maxLevel || (*r.region).ContainsCell(cell) {
				isTerminal = true
			}
		}
	}

	candidate := &Candidate{cell: cell, isTerminal: isTerminal}
	if !isTerminal {
		candidate.children = make([]*Candidate, 0, 8<<uint(r.maxChildrenShift()))
	}
	r.candidatesCreatedCounter++
	return candidate
}

func (r *RegionCoverer) ExpandChildren(candidate *Candidate, cell Cell, numLevels int) int {
	numLevels--
	childCells := make([]Cell, 0, 4)
	cell.Subdivide(&childCells)
	numTerminals := 0
	for _, c := range childCells {
		if numLevels > 0 {
			if (*r.region).MayIntersect(c) {
				numTerminals += r.ExpandChildren(candidate, c, numLevels)
			}
			continue
		}
		child := r.NewCandidate(c)
		if child != nil {
			candidate.children = append(candidate.children, child)
			if child.isTerminal {
				numTerminals++
			}
		}
	}
	return numTerminals
}

func (r *RegionCoverer) DeleteCandidate(candidate **Candidate, deleteChildren bool) {
	if deleteChildren {
		for _, c := range (*candidate).children {
			r.DeleteCandidate(&c, true)
		}
	}
	*candidate = nil
}

func (r *RegionCoverer) AddCandidate(candidate *Candidate) {
	if candidate == nil {
		return
	}
	if candidate.isTerminal {
		r.result = append(r.result, candidate.cell.id)
		r.DeleteCandidate(&candidate, true)
		return
	}

	// Expand one level at a time until we hit minLevel to ensure that we
	// don't skip over it.
	var numLevels int
	level := int(candidate.cell.level)
	if level < r.minLevel {
		numLevels = 1
	} else {
		numLevels = r.levelMod
	}
	numTerminals := r.ExpandChildren(candidate, candidate.cell, numLevels)
	shift := uint(r.maxChildrenShift())
	numChildren := len(candidate.children)

	if numChildren == 0 {
		r.DeleteCandidate(&candidate, false)
	} else if !r.interiorCovering &&
		numTerminals == 1<<shift &&
		level >= r.minLevel {
		// Optimization: add the parent cell rather than all of its
		// children. We can't do this for interior coverings, since
		// the children just intersect the region, but may not be
		// contained by it - we need to subdivide them further.
		candidate.isTerminal = true
		r.AddCandidate(candidate)
	} else {
		// We negate the priority so that smaller absolute priorities
		// are returned first. The heuristic is designed to refine the
		// largest cells first, since those are where we have the
		// largest potential gain. Among cells at the same level, we
		// prefer the cells with the smallest number of intersecting
		// children. Finally, we prefer cells that have the smallest
		// number of children that cannot be refined any further.
		priority := -((((level << shift) + numChildren) << shift) + numTerminals)
		heap.Push(r.pq, &CandidateEntry{candidate, priority})
		//		fmt.Printf("Push: %v (%v)\n", candidate.cell.id, priority)
	}
}

func (r *RegionCoverer) InitialCandidates() {
	// Optimization: if at least 4 cells are desired (the normal case),
	// start with a 4-cell covering of the region's bounding cap. This
	// lets us skip quite a few levels of refinement when the region to
	// be covered is relatively small.
	if r.maxCells >= 4 {
		// Find the maximum level such that the bounding cap contains
		// at most one cell vertex at that level.
		s2cap := (*r.region).CapBound()
		level := min(MinWidth.MaxLevel(2*float64(s2cap.Radius())), min(r.maxLevel, maxLevel-1))
		if r.levelMod > 1 && level > r.minLevel {
			level -= (level - r.minLevel) % r.levelMod
		}
		// We don't bother trying to optimize the level == 0 case,
		// since more than four face cells may be required.
		if level > 0 {
			// Find the leaf cell containing the cap axis, and
			// determine which subcell of the parent cell contains
			// it.
			base := make([]CellID, 0, 4)
			id := cellIDFromPoint(s2cap.center)
			id.AppendVertexNeighbors(level, &base)
			for _, cid := range base {
				r.AddCandidate(r.NewCandidate(CellFromCellID(cid)))
			}
			return
		}
	}
	// Default: start with all six cube faces.
	for face := 0; face < 6; face++ {
		r.AddCandidate(r.NewCandidate(faceCells[face]))
	}
}

func (r *RegionCoverer) CoveringInternal(region Region) {
	// Strategy: Start with the 6 faces of the cube. Discard any that do
	// not intersect the shape. Then repeatedly choose the largest cell
	// that intersects the shape and subdivide it.
	//
	// r.result contains the cells that will be part of the output, while
	// r.pq contains cells that we may still need to subdivide further.
	// Cells that are entirely contained within the region are immediately
	// added to the output, while cells that do not intersect the region
	// are immediately discarded. Therefore r.pq only contains cells that
	// partially intersect the region. Candidates are prioritized first
	// according to cell size (larger cells first), then by the number of
	// intersecting children they have (fewest children first), and then
	// by the number of fully contained children (fewest children first).
	r.region = &region
	r.candidatesCreatedCounter = 0

	r.InitialCandidates()
	for r.pq.Len() > 0 && (!r.interiorCovering || len(r.result) < r.maxCells) {
		candidate := heap.Pop(r.pq).(*CandidateEntry).candidate
		numChildren := len(candidate.children)
		count := 0
		if !r.interiorCovering {
			count = r.pq.Len()
		}
		if int(candidate.cell.level) < r.minLevel ||
			numChildren == 1 ||
			len(r.result)+numChildren+count <= r.maxCells {
			// Expand this candidate into its children.
			for _, child := range candidate.children {
				r.AddCandidate(child)
			}
			r.DeleteCandidate(&candidate, false)
		} else if r.interiorCovering {
			r.DeleteCandidate(&candidate, true)
		} else {
			candidate.isTerminal = true
			r.AddCandidate(candidate)
		}
	}
	for r.pq.Len() > 0 {
		candidate := heap.Pop(r.pq).(*CandidateEntry).candidate
		r.DeleteCandidate(&candidate, true)
	}
	r.region = nil
}

func (r *RegionCoverer) Covering(region Region) []CellID {
	// Rather than just returning the raw list of cell ids generated by
	// GetCoveringInternal(), we construct a cell union and then
	// denormalize it. This has the effect of replacing four child cells
	// with their parent whenever this does not violate the covering
	// parameters specified (minLevel, levelMod, etc). This strategy
	// significantly reduces the number of cells returned in many cases,
	// and it is cheap compared to computing the covering in the first
	// place.
	cu := r.CellUnionCovering(region)
	return cu.Denormalize(r.minLevel, r.levelMod)
}

func (r *RegionCoverer) CellUnionCovering(region Region) CellUnion {
	covering := CellUnion{}
	r.interiorCovering = false
	r.CoveringInternal(region)
	covering.InitSwap(&r.result)
	return covering
}

func FloodFill(region Region, start CellID) []CellID {
	all := map[CellID]bool{}
	frontier := []CellID{}
	output := []CellID{}
	all[start] = true
	frontier = append(frontier, start)
	for len(frontier) > 0 {
		back := len(frontier) - 1
		id := frontier[back]
		frontier = frontier[:back]
		if !region.MayIntersect(CellFromCellID(id)) {
			continue
		}
		output = append(output, id)
		neighbors := id.EdgeNeighbors()
		for _, nbr := range neighbors {
			if _, ok := all[nbr]; !ok {
				all[nbr] = true
				frontier = append(frontier, nbr)
			}
		}
	}
	return output
}

func SimpleCovering(region Region, start Point, level int) []CellID {
	return FloodFill(region, cellIDFromPoint(start).Parent(level))
}

var (
	faceCells [6]Cell
)

func init() {
	for face := 0; face < 6; face++ {
		faceCells[face] = CellFromCellID(CellIDFromFacePosLevel(face, 0, 0))
	}
}
