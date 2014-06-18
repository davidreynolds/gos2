package s2

import (
	"sort"
)

// A CellUnion is a collection of CellIDs.
//
// It is normalized if it is sorted, and does not contain redundancy.
// Specifically, it may not contain the same CellID twice, nor a CellID that is contained by another,
// nor the four sibling CellIDs that are children of a single higher level CellID.
type CellUnion []CellID

func (cu *CellUnion) Init(cellIds []CellID) {
	cu.InitRaw(cellIds)
	cu.Normalize()
}

func (cu *CellUnion) InitRaw(cellIds []CellID) {
	*cu = cellIds
}

func (cu *CellUnion) InitSwap(cellIds *[]CellID) {
	cu.InitRawSwap(cellIds)
	cu.Normalize()
}

func (cu *CellUnion) InitRawSwap(cellIds *[]CellID) {
	*cu = *cellIds
	*cellIds = []CellID{}
}

// Normalize normalizes the CellUnion.
func (cu *CellUnion) Normalize() {
	sort.Sort(byID(*cu))

	output := make([]CellID, 0, len(*cu)) // the list of accepted cells
	// Loop invariant: output is a sorted list of cells with no redundancy.
	for _, ci := range *cu {
		// The first two passes here either ignore this new candidate,
		// or remove previously accepted cells that are covered by this candidate.

		// Ignore this cell if it is contained by the previous one.
		// We only need to check the last accepted cell. The ordering of the
		// cells implies containment (but not the converse), and output has no redundancy,
		// so if this candidate is not contained by the last accepted cell
		// then it cannot be contained by any previously accepted cell.
		if len(output) > 0 && output[len(output)-1].Contains(ci) {
			continue
		}

		// Discard any previously accepted cells contained by this one.
		// This could be any contiguous trailing subsequence, but it can't be
		// a discontiguous subsequence because of the containment property of
		// sorted S2 cells mentioned above.
		j := len(output) - 1 // last index to keep
		for j >= 0 {
			if !ci.Contains(output[j]) {
				break
			}
			j--
		}
		output = output[:j+1]

		// See if the last three cells plus this one can be collapsed.
		// We loop because collapsing three accepted cells and adding a higher level cell
		// could cascade into previously accepted cells.
		for len(output) >= 3 {
			fin := output[len(output)-3:]

			// fast XOR test; a necessary but not sufficient condition
			if fin[0]^fin[1]^fin[2]^ci != 0 {
				break
			}

			// more expensive test; exact.
			// Compute the two bit mask for the encoded child position,
			// then see if they all agree.
			mask := CellID(ci.lsb() << 1)
			mask = ^(mask + mask<<1)
			should := ci & mask
			if (fin[0]&mask != should) || (fin[1]&mask != should) || (fin[2]&mask != should) || ci.isFace() {
				break
			}

			output = output[:len(output)-3]
			ci = ci.immediateParent() // checked !ci.isFace above
		}
		output = append(output, ci)
	}
	*cu = output
}

func (cu CellUnion) Denormalize(minLevel, levelMod int) []CellID {
	output := make([]CellID, 0, len(cu))
	for _, id := range cu {
		level := id.Level()
		newLevel := max(minLevel, level)
		if levelMod > 1 {
			// Round up so that (newLevel - minLevel) is a multiple
			// of levelMod. (Note that maxLevel is a multiple
			// of 1, 2, and 3.)
			newLevel += (maxLevel - (newLevel - minLevel)) % levelMod
			newLevel = min(maxLevel, newLevel)
		}
		if newLevel == level {
			output = append(output, id)
		} else {
			end := id.ChildEndAtLevel(newLevel)
			for id = id.ChildBeginAtLevel(newLevel); id != end; id = id.Next() {
				output = append(output, id)
			}
		}
	}
	return output
}

func (cu CellUnion) ContainsCellID(id CellID) bool {
	// This function requires that Normalize has been called first.
	//
	// This is an exact test. Each cell occupies a linear span of the S2
	// space-filling curve, and the cell id is simply the position at the
	// center of this span. The cell union ids are sorted in increasing
	// order along the space-filling curve. So we simply find the pair of
	// cell ids that surround the given cell id (using binary search).
	// There is containment if and only if one of these two cell ids
	// contains this cell.
	idx := sort.Search(len(cu), func(i int) bool { return cu[i] >= id })
	if idx < len(cu) && cu[idx].RangeMin() <= id {
		return true
	}
	return idx > 0 && cu[idx-1].RangeMax() >= id
}

func (cu CellUnion) IntersectsCellID(id CellID) bool {
	// This function requires that Normalize has been called first.
	// This is an exact test; see ContainsCellID above.
	idx := sort.Search(len(cu), func(i int) bool { return cu[i] >= id })
	if idx < len(cu) && cu[idx].RangeMin() <= id.RangeMax() {
		return true
	}
	return idx > 0 && cu[idx-1].RangeMax() >= id.RangeMin()
}

type byID []CellID

func (cu byID) Len() int           { return len(cu) }
func (cu byID) Less(i, j int) bool { return uint64(cu[i]) < uint64(cu[j]) }
func (cu byID) Swap(i, j int)      { cu[i], cu[j] = cu[j], cu[i] }
