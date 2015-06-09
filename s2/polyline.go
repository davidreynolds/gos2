package s2

type Polyline struct {
	vertices []Point
}

func PolylineFromPoints(vertices []Point) *Polyline {
	pl := &Polyline{}
	pl.initPoints(vertices)
	return pl
}

func (pl *Polyline) initPoints(vertices []Point) {
	pl.vertices = make([]Point, len(vertices))
	copy(pl.vertices, vertices)
}

func (pl Polyline) NumVertices() int   { return len(pl.vertices) }
func (pl Polyline) Vertex(k int) Point { return pl.vertices[k] }

func (pl *Polyline) GetRectBound() *Rect {
	rb := NewRectBounder()
	for _, pt := range pl.vertices {
		rb.AddPoint(&pt)
	}
	bound := rb.Bound()
	return &bound
}

func (pl *Polyline) MayIntersect(cell Cell) bool {
	if pl.NumVertices() == 0 {
		return false
	}

	for _, pt := range pl.vertices {
		if cell.ContainsPoint(pt) {
			return true
		}
	}

	cellVertices := []Point{
		cell.Vertex(0),
		cell.Vertex(1),
		cell.Vertex(2),
		cell.Vertex(3),
	}

	for j := 0; j < 4; j++ {
		crosser := NewEdgeCrosser(&cellVertices[j], &cellVertices[(j+1)&3], &pl.vertices[0])
		for _, vetex := range pl.vertices[1:] {
			if crosser.RobustCrossing(&vetex) >= 0 {
				return true
			}
		}
	}

	return false
}

func (pl *Polyline) ContainsCell(cell Cell) bool { return false }
func (pl *Polyline) ContainsPoint(p Point) bool  { return false }

func (pl *Polyline) CapBound() Cap {
	return pl.GetRectBound().CapBound()
}
