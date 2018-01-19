// ContainsPoint returns true if the loop contains the point.
func (l *Loop) ContainsPoint(p Point) bool {
// Empty and full loops don't need a special case, but invalid loops with
// zero vertices do, so we might as well handle them all at once.
if len(l.vertices) < 3 {
return l.originInside
}

// For small loops, and during initial construction, it is faster to just
// check all the crossing.
const maxBruteForceVertices = 32
if len(l.vertices) < maxBruteForceVertices || l.index == nil {
return l.bruteForceContainsPoint(p)
}

// Otherwise, look up the point in the index.
it := l.index.Iterator()
if !it.LocatePoint(p) {
return false
}
return l.iteratorContainsPoint(it, p)
}

// ContainsCell reports whether the given Cell is contained by this Loop.
func (l *Loop) ContainsCell(target Cell) bool {
it := l.index.Iterator()
relation := it.LocateCellID(target.ID())

// If "target" is disjoint from all index cells, it is not contained.
// Similarly, if "target" is subdivided into one or more index cells then it
// is not contained, since index cells are subdivided only if they (nearly)
// intersect a sufficient number of edges.  (But note that if "target" itself
// is an index cell then it may be contained, since it could be a cell with
// no edges in the loop interior.)
if relation != Indexed {
return false
}

// Otherwise check if any edges intersect "target".
if l.boundaryApproxIntersects(it, target) {
return false
}

// Otherwise check if the loop contains the center of "target".
return l.iteratorContainsPoint(it, target.Center())
}

// IntersectsCell reports whether this Loop intersects the given cell.
func (l *Loop) IntersectsCell(target Cell) bool {
it := l.index.Iterator()
relation := it.LocateCellID(target.ID())

// If target does not overlap any index cell, there is no intersection.
if relation == Disjoint {
return false
}
// If target is subdivided into one or more index cells, there is an
// intersection to within the ShapeIndex error bound (see Contains).
if relation == Subdivided {
return true
}
// If target is an index cell, there is an intersection because index cells
// are created only if they have at least one edge or they are entirely
// contained by the loop.
if it.CellID() == target.id {
return true
}
// Otherwise check if any edges intersect target.
if l.boundaryApproxIntersects(it, target) {
return true
}
// Otherwise check if the loop contains the center of target.
return l.iteratorContainsPoint(it, target.Center())
}

// CellUnionBound computes a covering of the Loop.
func (l *Loop) CellUnionBound() []CellID {
return l.CapBound().CellUnionBound()
}

// boundaryApproxIntersects reports if the loop's boundary intersects target.
// It may also return true when the loop boundary does not intersect target but
// some edge comes within the worst-case error tolerance.
//
// This requires that it.Locate(target) returned Indexed.
func (l *Loop) boundaryApproxIntersects(it *ShapeIndexIterator, target Cell) bool {
aClipped := it.IndexCell().findByShapeID(0)

// If there are no edges, there is no intersection.
if len(aClipped.edges) == 0 {
return false
}

// We can save some work if target is the index cell itself.
if it.CellID() == target.ID() {
return true
}

// Otherwise check whether any of the edges intersect target.
maxError := (faceClipErrorUVCoord + intersectsRectErrorUVDist)
bound := target.BoundUV().ExpandedByMargin(maxError)
for _, ai := range aClipped.edges {
v0, v1, ok := ClipToPaddedFace(l.Vertex(ai), l.Vertex(ai+1), target.Face(), maxError)
if ok && edgeIntersectsRect(v0, v1, bound) {
return true
}
}
return false
}

// iteratorContainsPoint reports if the iterator that is positioned at the ShapeIndexCell
// that may contain p, contains the point p.
func (l *Loop) iteratorContainsPoint(it *ShapeIndexIterator, p Point) bool {
// Test containment by drawing a line segment from the cell center to the
// given point and counting edge crossings.
aClipped := it.IndexCell().findByShapeID(0)
inside := aClipped.containsCenter
if len(aClipped.edges) > 0 {
center := it.Center()
crosser := NewEdgeCrosser(center, p)
aiPrev := -2
for _, ai := range aClipped.edges {
if ai != aiPrev+1 {
crosser.RestartAt(l.Vertex(ai))
}
aiPrev = ai
inside = inside != crosser.EdgeOrVertexChainCrossing(l.Vertex(ai+1))
}
}
return inside
}

// RegularLoop creates a loop with the given number of vertices, all
// located on a circle of the specified radius around the given center.
func RegularLoop(center Point, radius s1.Angle, numVertices int) *Loop {
return RegularLoopForFrame(getFrame(center), radius, numVertices)
}

// RegularLoopForFrame creates a loop centered around the z-axis of the given
// coordinate frame, with the first vertex in the direction of the positive x-axis.
func RegularLoopForFrame(frame matrix3x3, radius s1.Angle, numVertices int) *Loop {
return LoopFromPoints(regularPointsForFrame(frame, radius, numVertices))
}

// CanonicalFirstVertex returns a first index and a direction (either +1 or -1)
// such that the vertex sequence (first, first+dir, ..., first+(n-1)*dir) does
// not change when the loop vertex order is rotated or inverted. This allows the
// loop vertices to be traversed in a canonical order. The return values are
// chosen such that (first, ..., first+n*dir) are in the range [0, 2*n-1] as
// expected by the Vertex method.
func (l *Loop) CanonicalFirstVertex() (firstIdx, direction int) {
firstIdx = 0
n := len(l.vertices)
for i := 1; i < n; i++ {
if l.Vertex(i).Cmp(l.Vertex(firstIdx).Vector) == -1 {
firstIdx = i
}
}

// 0 <= firstIdx <= n-1, so (firstIdx+n*dir) <= 2*n-1.
if l.Vertex(firstIdx+1).Cmp(l.Vertex(firstIdx+n-1).Vector) == -1 {
return firstIdx, 1
}

// n <= firstIdx <= 2*n-1, so (firstIdx+n*dir) >= 0.
firstIdx += n
return firstIdx, -1
}

// TurningAngle returns the sum of the turning angles at each vertex. The return
// value is positive if the loop is counter-clockwise, negative if the loop is
// clockwise, and zero if the loop is a great circle. Degenerate and
// nearly-degenerate loops are handled consistently with Sign. So for example,
// if a loop has zero area (i.e., it is a very small CCW loop) then the turning
// angle will always be negative.
//
// This quantity is also called the "geodesic curvature" of the loop.
func (l *Loop) TurningAngle() float64 {
// For empty and full loops, we return the limit value as the loop area
// approaches 0 or 4*Pi respectively.
if l.isEmptyOrFull() {
if l.ContainsOrigin() {
return -2 * math.Pi
}
return 2 * math.Pi
}

// Don't crash even if the loop is not well-defined.
if len(l.vertices) < 3 {
return 0
}

// To ensure that we get the same result when the vertex order is rotated,
// and that the result is negated when the vertex order is reversed, we need
// to add up the individual turn angles in a consistent order. (In general,
// adding up a set of numbers in a different order can change the sum due to
// rounding errors.)
//
// Furthermore, if we just accumulate an ordinary sum then the worst-case
// error is quadratic in the number of vertices. (This can happen with
// spiral shapes, where the partial sum of the turning angles can be linear
// in the number of vertices.) To avoid this we use the Kahan summation
// algorithm (http://en.wikipedia.org/wiki/Kahan_summation_algorithm).
n := len(l.vertices)
i, dir := l.CanonicalFirstVertex()
sum := TurnAngle(l.Vertex((i+n-dir)%n), l.Vertex(i), l.Vertex((i+dir)%n))

compensation := s1.Angle(0)
for n-1 > 0 {
i += dir
angle := TurnAngle(l.Vertex(i-dir), l.Vertex(i), l.Vertex(i+dir))
oldSum := sum
angle += compensation
sum += angle
compensation = (oldSum - sum) + angle
n--
}
return float64(dir) * float64(sum+compensation)
}

// turningAngleMaxError return the maximum error in TurningAngle. The value is not
// constant; it depends on the loop.
func (l *Loop) turningAngleMaxError() float64 {
// The maximum error can be bounded as follows:
//   2.24 * dblEpsilon    for RobustCrossProd(b, a)
//   2.24 * dblEpsilon    for RobustCrossProd(c, b)
//   3.25 * dblEpsilon    for Angle()
//   2.00 * dblEpsilon    for each addition in the Kahan summation
//   ------------------
//   9.73 * dblEpsilon
maxErrorPerVertex := 9.73 * dblEpsilon
return maxErrorPerVertex * float64(len(l.vertices))
}

// IsNormalized reports whether the loop area is at most 2*pi. Degenerate loops are
// handled consistently with Sign, i.e., if a loop can be
// expressed as the union of degenerate or nearly-degenerate CCW triangles,
// then it will always be considered normalized.
func (l *Loop) IsNormalized() bool {
// Optimization: if the longitude span is less than 180 degrees, then the
// loop covers less than half the sphere and is therefore normalized.
if l.bound.Lng.Length() < math.Pi {
return true
}

// We allow some error so that hemispheres are always considered normalized.
// TODO(roberts): This is no longer required by the Polygon implementation,
// so alternatively we could create the invariant that a loop is normalized
// if and only if its complement is not normalized.
return l.TurningAngle() >= -l.turningAngleMaxError()
}

// Normalize inverts the loop if necessary so that the area enclosed by the loop
// is at most 2*pi.
func (l *Loop) Normalize() {
if !l.IsNormalized() {
l.Invert()
}
}

// Invert reverses the order of the loop vertices, effectively complementing the
// region represented by the loop. For example, the loop ABCD (with edges
// AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
// Notice that the last edge is the same in both cases except that its
// direction has been reversed.
func (l *Loop) Invert() {
l.index.Reset()
if l.isEmptyOrFull() {
if l.IsFull() {
l.vertices[0] = emptyLoopPoint
} else {
l.vertices[0] = fullLoopPoint
}
} else {
// For non-special loops, reverse the slice of vertices.
for i := len(l.vertices)/2 - 1; i >= 0; i-- {
opp := len(l.vertices) - 1 - i
l.vertices[i], l.vertices[opp] = l.vertices[opp], l.vertices[i]
}
}

// originInside must be set correctly before building the ShapeIndex.
l.originInside = l.originInside != true
if l.bound.Lat.Lo > -math.Pi/2 && l.bound.Lat.Hi < math.Pi/2 {
// The complement of this loop contains both poles.
l.bound = FullRect()
l.subregionBound = l.bound
} else {
l.initBound()
}
l.index.Add(l)
}

// ContainsNested reports whether the given loops is contained within this loop.
// This function does not test for edge intersections. The two loops must meet
// all of the Polygon requirements; for example this implies that their
// boundaries may not cross or have any shared edges (although they may have
// shared vertices).
func (l *Loop) ContainsNested(other *Loop) bool {
if !l.subregionBound.Contains(other.bound) {
return false
}

// Special cases to handle either loop being empty or full.  Also bail out
// when B has no vertices to avoid heap overflow on the vertex(1) call
// below.  (This method is called during polygon initialization before the
// client has an opportunity to call IsValid().)
if l.isEmptyOrFull() || other.NumVertices() < 2 {
return l.IsFull() || other.IsEmpty()
}

// We are given that A and B do not share any edges, and that either one
// loop contains the other or they do not intersect.
m, ok := l.findVertex(other.Vertex(1))
if !ok {
// Since other.vertex(1) is not shared, we can check whether A contains it.
return l.ContainsPoint(other.Vertex(1))
}

// Check whether the edge order around other.Vertex(1) is compatible with
// A containing B.
return WedgeContains(l.Vertex(m-1), l.Vertex(m), l.Vertex(m+1), other.Vertex(0), other.Vertex(2))
}

// surfaceIntegralFloat64 computes the oriented surface integral of some quantity f(x)
// over the loop interior, given a function f(A,B,C) that returns the
// corresponding integral over the spherical triangle ABC. Here "oriented
// surface integral" means:
//
// (1) f(A,B,C) must be the integral of f if ABC is counterclockwise,
//     and the integral of -f if ABC is clockwise.
//
// (2) The result of this function is *either* the integral of f over the
//     loop interior, or the integral of (-f) over the loop exterior.
//
// Note that there are at least two common situations where it easy to work
// around property (2) above:
//
//  - If the integral of f over the entire sphere is zero, then it doesn't
//    matter which case is returned because they are always equal.
//
//  - If f is non-negative, then it is easy to detect when the integral over
//    the loop exterior has been returned, and the integral over the loop
//    interior can be obtained by adding the integral of f over the entire
//    unit sphere (a constant) to the result.
//
// Any changes to this method may need corresponding changes to surfaceIntegralPoint as well.
func (l *Loop) surfaceIntegralFloat64(f func(a, b, c Point) float64) float64 {
// We sum f over a collection T of oriented triangles, possibly
// overlapping. Let the sign of a triangle be +1 if it is CCW and -1
// otherwise, and let the sign of a point x be the sum of the signs of the
// triangles containing x. Then the collection of triangles T is chosen
// such that either:
//
//  (1) Each point in the loop interior has sign +1, and sign 0 otherwise; or
//  (2) Each point in the loop exterior has sign -1, and sign 0 otherwise.
//
// The triangles basically consist of a fan from vertex 0 to every loop
// edge that does not include vertex 0. These triangles will always satisfy
// either (1) or (2). However, what makes this a bit tricky is that
// spherical edges become numerically unstable as their length approaches
// 180 degrees. Of course there is not much we can do if the loop itself
// contains such edges, but we would like to make sure that all the triangle
// edges under our control (i.e., the non-loop edges) are stable. For
// example, consider a loop around the equator consisting of four equally
// spaced points. This is a well-defined loop, but we cannot just split it
// into two triangles by connecting vertex 0 to vertex 2.
//
// We handle this type of situation by moving the origin of the triangle fan
// whenever we are about to create an unstable edge. We choose a new
// location for the origin such that all relevant edges are stable. We also
// create extra triangles with the appropriate orientation so that the sum
// of the triangle signs is still correct at every point.

// The maximum length of an edge for it to be considered numerically stable.
// The exact value is fairly arbitrary since it depends on the stability of
// the function f. The value below is quite conservative but could be
// reduced further if desired.
const maxLength = math.Pi - 1e-5

var sum float64
origin := l.Vertex(0)
for i := 1; i+1 < len(l.vertices); i++ {
// Let V_i be vertex(i), let O be the current origin, and let length(A,B)
// be the length of edge (A,B). At the start of each loop iteration, the
// "leading edge" of the triangle fan is (O,V_i), and we want to extend
// the triangle fan so that the leading edge is (O,V_i+1).
//
// Invariants:
//  1. length(O,V_i) < maxLength for all (i > 1).
//  2. Either O == V_0, or O is approximately perpendicular to V_0.
//  3. "sum" is the oriented integral of f over the area defined by
//     (O, V_0, V_1, ..., V_i).
if l.Vertex(i+1).Angle(origin.Vector) > maxLength {
// We are about to create an unstable edge, so choose a new origin O'
// for the triangle fan.
oldOrigin := origin
if origin == l.Vertex(0) {
// The following point is well-separated from V_i and V_0 (and
// therefore V_i+1 as well).
origin = Point{l.Vertex(0).PointCross(l.Vertex(i)).Normalize()}
} else if l.Vertex(i).Angle(l.Vertex(0).Vector) < maxLength {
// All edges of the triangle (O, V_0, V_i) are stable, so we can
// revert to using V_0 as the origin.
origin = l.Vertex(0)
} else {
// (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0 are
// perpendicular. Therefore V_0.CrossProd(O) is approximately
// perpendicular to all of {O, V_0, V_i, V_i+1}, and we can choose
// this point O' as the new origin.
origin = Point{l.Vertex(0).Cross(oldOrigin.Vector)}

// Advance the edge (V_0,O) to (V_0,O').
sum += f(l.Vertex(0), oldOrigin, origin)
}
// Advance the edge (O,V_i) to (O',V_i).
sum += f(oldOrigin, l.Vertex(i), origin)
}
// Advance the edge (O,V_i) to (O,V_i+1).
sum += f(origin, l.Vertex(i), l.Vertex(i+1))
}
// If the origin is not V_0, we need to sum one more triangle.
if origin != l.Vertex(0) {
// Advance the edge (O,V_n-1) to (O,V_0).
sum += f(origin, l.Vertex(len(l.vertices)-1), l.Vertex(0))
}
return sum
}

// surfaceIntegralPoint mirrors the surfaceIntegralFloat64 method but over Points;
// see that method for commentary. The C++ version uses a templated method.
// Any changes to this method may need corresponding changes to surfaceIntegralFloat64 as well.
func (l *Loop) surfaceIntegralPoint(f func(a, b, c Point) Point) Point {
const maxLength = math.Pi - 1e-5
var sum r3.Vector

origin := l.Vertex(0)
for i := 1; i+1 < len(l.vertices); i++ {
if l.Vertex(i+1).Angle(origin.Vector) > maxLength {
oldOrigin := origin
if origin == l.Vertex(0) {
origin = Point{l.Vertex(0).PointCross(l.Vertex(i)).Normalize()}
} else if l.Vertex(i).Angle(l.Vertex(0).Vector) < maxLength {
origin = l.Vertex(0)
} else {
origin = Point{l.Vertex(0).Cross(oldOrigin.Vector)}
sum = sum.Add(f(l.Vertex(0), oldOrigin, origin).Vector)
}
sum = sum.Add(f(oldOrigin, l.Vertex(i), origin).Vector)
}
sum = sum.Add(f(origin, l.Vertex(i), l.Vertex(i+1)).Vector)
}
if origin != l.Vertex(0) {
sum = sum.Add(f(origin, l.Vertex(len(l.vertices)-1), l.Vertex(0)).Vector)
}
return Point{sum}
}

// Area returns the area of the loop interior, i.e. the region on the left side of
// the loop. The return value is between 0 and 4*pi. (Note that the return
// value is not affected by whether this loop is a "hole" or a "shell".)
func (l *Loop) Area() float64 {
// It is suprisingly difficult to compute the area of a loop robustly. The
// main issues are (1) whether degenerate loops are considered to be CCW or
// not (i.e., whether their area is close to 0 or 4*pi), and (2) computing
// the areas of small loops with good relative accuracy.
//
// With respect to degeneracies, we would like Area to be consistent
// with ContainsPoint in that loops that contain many points
// should have large areas, and loops that contain few points should have
// small areas. For example, if a degenerate triangle is considered CCW
// according to s2predicates Sign, then it will contain very few points and
// its area should be approximately zero. On the other hand if it is
// considered clockwise, then it will contain virtually all points and so
// its area should be approximately 4*pi.
//
// More precisely, let U be the set of Points for which IsUnitLength
// is true, let P(U) be the projection of those points onto the mathematical
// unit sphere, and let V(P(U)) be the Voronoi diagram of the projected
// points. Then for every loop x, we would like Area to approximately
// equal the sum of the areas of the Voronoi regions of the points p for
// which x.ContainsPoint(p) is true.
//
// The second issue is that we want to compute the area of small loops
// accurately. This requires having good relative precision rather than
// good absolute precision. For example, if the area of a loop is 1e-12 and
// the error is 1e-15, then the area only has 3 digits of accuracy. (For
// reference, 1e-12 is about 40 square meters on the surface of the earth.)
// We would like to have good relative accuracy even for small loops.
//
// To achieve these goals, we combine two different methods of computing the
// area. This first method is based on the Gauss-Bonnet theorem, which says
// that the area enclosed by the loop equals 2*pi minus the total geodesic
// curvature of the loop (i.e., the sum of the "turning angles" at all the
// loop vertices). The big advantage of this method is that as long as we
// use Sign to compute the turning angle at each vertex, then
// degeneracies are always handled correctly. In other words, if a
// degenerate loop is CCW according to the symbolic perturbations used by
// Sign, then its turning angle will be approximately 2*pi.
//
// The disadvantage of the Gauss-Bonnet method is that its absolute error is
// about 2e-15 times the number of vertices (see turningAngleMaxError).
// So, it cannot compute the area of small loops accurately.
//
// The second method is based on splitting the loop into triangles and
// summing the area of each triangle. To avoid the difficulty and expense
// of decomposing the loop into a union of non-overlapping triangles,
// instead we compute a signed sum over triangles that may overlap (see the
// comments for surfaceIntegral). The advantage of this method
// is that the area of each triangle can be computed with much better
// relative accuracy (using l'Huilier's theorem). The disadvantage is that
// the result is a signed area: CCW loops may yield a small positive value,
// while CW loops may yield a small negative value (which is converted to a
// positive area by adding 4*pi). This means that small errors in computing
// the signed area may translate into a very large error in the result (if
// the sign of the sum is incorrect).
//
// So, our strategy is to combine these two methods as follows. First we
// compute the area using the "signed sum over triangles" approach (since it
// is generally more accurate). We also estimate the maximum error in this
// result. If the signed area is too close to zero (i.e., zero is within
// the error bounds), then we double-check the sign of the result using the
// Gauss-Bonnet method. (In fact we just call IsNormalized, which is
// based on this method.) If the two methods disagree, we return either 0
// or 4*pi based on the result of IsNormalized. Otherwise we return the
// area that we computed originally.
if l.isEmptyOrFull() {
if l.ContainsOrigin() {
return 4 * math.Pi
}
return 0
}
area := l.surfaceIntegralFloat64(SignedArea)

// TODO(roberts): This error estimate is very approximate. There are two
// issues: (1) SignedArea needs some improvements to ensure that its error
// is actually never higher than GirardArea, and (2) although the number of
// triangles in the sum is typically N-2, in theory it could be as high as
// 2*N for pathological inputs. But in other respects this error bound is
// very conservative since it assumes that the maximum error is achieved on
// every triangle.
maxError := l.turningAngleMaxError()

// The signed area should be between approximately -4*pi and 4*pi.
if area < 0 {
// We have computed the negative of the area of the loop exterior.
area += 4 * math.Pi
}

if area > 4*math.Pi {
area = 4 * math.Pi
}
if area < 0 {
area = 0
}

// If the area is close enough to zero or 4*pi so that the loop orientation
// is ambiguous, then we compute the loop orientation explicitly.
if area < maxError && !l.IsNormalized() {
return 4 * math.Pi
} else if area > (4*math.Pi-maxError) && l.IsNormalized() {
return 0
}

return area
}

// Centroid returns the true centroid of the loop multiplied by the area of the
// loop. The result is not unit length, so you may want to normalize it. Also
// note that in general, the centroid may not be contained by the loop.
//
// We prescale by the loop area for two reasons: (1) it is cheaper to
// compute this way, and (2) it makes it easier to compute the centroid of
// more complicated shapes (by splitting them into disjoint regions and
// adding their centroids).
//
// Note that the return value is not affected by whether this loop is a
// "hole" or a "shell".
func (l *Loop) Centroid() Point {
// surfaceIntegralPoint() returns either the integral of position over loop
// interior, or the negative of the integral of position over the loop
// exterior. But these two values are the same (!), because the integral of
// position over the entire sphere is (0, 0, 0).
return l.surfaceIntegralPoint(TrueCentroid)
}

// Encode encodes the Loop.
func (l Loop) Encode(w io.Writer) error {
e := &encoder{w: w}
l.encode(e)
return e.err
}

func (l Loop) encode(e *encoder) {
e.writeInt8(encodingVersion)
e.writeUint32(uint32(len(l.vertices)))
for _, v := range l.vertices {
e.writeFloat64(v.X)
e.writeFloat64(v.Y)
e.writeFloat64(v.Z)
}

e.writeBool(l.originInside)
e.writeInt32(int32(l.depth))

// Encode the bound.
l.bound.encode(e)
}

// Decode decodes a loop.
func (l *Loop) Decode(r io.Reader) error {
*l = Loop{}
d := &decoder{r: asByteReader(r)}
l.decode(d)
return d.err
}

func (l *Loop) decode(d *decoder) {
version := int8(d.readUint8())
if d.err != nil {
return
}
if version != encodingVersion {
d.err = fmt.Errorf("cannot decode version %d", version)
return
}

// Empty loops are explicitly allowed here: a newly created loop has zero vertices
// and such loops encode and decode properly.
nvertices := d.readUint32()
if nvertices > maxEncodedVertices {
if d.err == nil {
d.err = fmt.Errorf("too many vertices (%d; max is %d)", nvertices, maxEncodedVertices)

}
return
}
l.vertices = make([]Point, nvertices)
for i := range l.vertices {
l.vertices[i].X = d.readFloat64()
l.vertices[i].Y = d.readFloat64()
l.vertices[i].Z = d.readFloat64()
}
l.originInside = d.readBool()
l.depth = int(d.readUint32())
l.bound.decode(d)
l.subregionBound = ExpandForSubregions(l.bound)

l.index = NewShapeIndex()
l.index.Add(l)
}

// Bitmasks to read from properties.
const (
originInside = 1 << iota
boundEncoded
)

func (l *Loop) xyzFaceSiTiVertices() []xyzFaceSiTi {
ret := make([]xyzFaceSiTi, len(l.vertices))
for i, v := range l.vertices {
ret[i].xyz = v
ret[i].face, ret[i].si, ret[i].ti, ret[i].level = xyzToFaceSiTi(v)
}
return ret
}

func (l *Loop) encodeCompressed(e *encoder, snapLevel int, vertices []xyzFaceSiTi) {
if len(l.vertices) != len(vertices) {
panic("encodeCompressed: vertices must be the same length as l.vertices")
}
if len(vertices) > maxEncodedVertices {
if e.err == nil {
e.err = fmt.Errorf("too many vertices (%d; max is %d)", len(vertices), maxEncodedVertices)
}
return
}
e.writeUvarint(uint64(len(vertices)))
encodePointsCompressed(e, vertices, snapLevel)

props := l.compressedEncodingProperties()
e.writeUvarint(props)
e.writeUvarint(uint64(l.depth))
if props&boundEncoded != 0 {
l.bound.encode(e)
}
}

func (l *Loop) compressedEncodingProperties() uint64 {
var properties uint64
if l.originInside {
properties |= originInside
}

// Write whether there is a bound so we can change the threshold later.
// Recomputing the bound multiplies the decode time taken per vertex
// by a factor of about 3.5.  Without recomputing the bound, decode
// takes approximately 125 ns / vertex.  A loop with 63 vertices
// encoded without the bound will take ~30us to decode, which is
// acceptable.  At ~3.5 bytes / vertex without the bound, adding
// the bound will increase the size by <15%, which is also acceptable.
const minVerticesForBound = 64
if len(l.vertices) >= minVerticesForBound {
properties |= boundEncoded
}

return properties
}

func (l *Loop) decodeCompressed(d *decoder, snapLevel int) {
nvertices := d.readUvarint()
if d.err != nil {
return
}
if nvertices > maxEncodedVertices {
d.err = fmt.Errorf("too many vertices (%d; max is %d)", nvertices, maxEncodedVertices)
return
}
l.vertices = make([]Point, nvertices)
decodePointsCompressed(d, snapLevel, l.vertices)
properties := d.readUvarint()

// Make sure values are valid before using.
if d.err != nil {
return
}

l.originInside = (properties & originInside) != 0

l.depth = int(d.readUvarint())

if (properties & boundEncoded) != 0 {
l.bound.decode(d)
if d.err != nil {
return
}
l.subregionBound = ExpandForSubregions(l.bound)
} else {
l.initBound()
}

l.index = NewShapeIndex()
l.index.Add(l)
}
