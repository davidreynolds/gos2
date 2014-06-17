package s2

type Region interface {
	MayIntersect(cell Cell) bool
	ContainsCell(cell Cell) bool
	ContainsPoint(p Point) bool
	CapBound() Cap
}
