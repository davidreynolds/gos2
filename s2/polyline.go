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
