func TestLoopAreaAndCentroid(t *testing.T) {
var p Point

if got, want := EmptyLoop().Area(), 0.0; got != want {
t.Errorf("EmptyLoop.Area() = %v, want %v", got, want)
}
if got, want := FullLoop().Area(), 4*math.Pi; got != want {
t.Errorf("FullLoop.Area() = %v, want %v", got, want)
}
if got := EmptyLoop().Centroid(); !p.ApproxEqual(got) {
t.Errorf("EmptyLoop.Centroid() = %v, want %v", got, p)
}
if got := FullLoop().Centroid(); !p.ApproxEqual(got) {
t.Errorf("FullLoop.Centroid() = %v, want %v", got, p)
}

if got, want := northHemi.Area(), 2*math.Pi; !float64Eq(got, want) {
t.Errorf("northHemi.Area() = %v, want %v", got, want)
}

eastHemiArea := eastHemi.Area()
if eastHemiArea < 2*math.Pi-1e-12 || eastHemiArea > 2*math.Pi+1e-12 {
t.Errorf("eastHemi.Area() = %v, want between [%v, %v]", eastHemiArea, 2*math.Pi-1e-12, 2*math.Pi+1e-12)
}

// Construct spherical caps of random height, and approximate their boundary
// with closely spaces vertices. Then check that the area and centroid are
// correct.
for i := 0; i < 50; i++ {
// Choose a coordinate frame for the spherical cap.
f := randomFrame()
x := f.col(0)
y := f.col(1)
z := f.col(2)

// Given two points at latitude phi and whose longitudes differ by dtheta,
// the geodesic between the two points has a maximum latitude of
// atan(tan(phi) / cos(dtheta/2)). This can be derived by positioning
// the two points at (-dtheta/2, phi) and (dtheta/2, phi).
//
// We want to position the vertices close enough together so that their
// maximum distance from the boundary of the spherical cap is maxDist.
// Thus we want abs(atan(tan(phi) / cos(dtheta/2)) - phi) <= maxDist.
const maxDist = 1e-6
height := 2 * randomFloat64()
phi := math.Asin(1.0 - height)
maxDtheta := 2 * math.Acos(math.Tan(math.Abs(phi))/math.Tan(math.Abs(phi)+maxDist))
maxDtheta = math.Min(math.Pi, maxDtheta)

var vertices []Point
for theta := 0.0; theta < 2*math.Pi; theta += randomFloat64() * maxDtheta {
vertices = append(vertices,
Point{x.Mul(math.Cos(theta) * math.Cos(phi)).Add(
y.Mul(math.Sin(theta) * math.Cos(phi))).Add(
z.Mul(math.Sin(phi)))},
)
}

loop := LoopFromPoints(vertices)
area := loop.Area()
centroid := loop.Centroid()
expectedArea := 2 * math.Pi * height
if delta, want := math.Abs(area-expectedArea), 2*math.Pi*maxDist; delta > want {
t.Errorf("%v.Area() = %v, want to be within %v of %v, got %v", loop, area, want, expectedArea, delta)
}
expectedCentroid := z.Mul(expectedArea * (1 - 0.5*height))
if delta, want := centroid.Sub(expectedCentroid).Norm(), 2*maxDist; delta > want {
t.Errorf("%v.Centroid() = %v, want to be within %v of %v, got %v", loop, centroid, want, expectedCentroid, delta)
}
}
}

// TODO(roberts): Test that Area() has an accuracy significantly better
// than 1e-15 on loops whose area is small.

func TestLoopAreaConsistentWithTurningAngle(t *testing.T) {
// Check that the area computed using GetArea() is consistent with the
// turning angle of the loop computed using GetTurnAngle().  According to
// the Gauss-Bonnet theorem, the area of the loop should be equal to 2*Pi
// minus its turning angle.
for x, loop := range allLoops {
area := loop.Area()
gaussArea := 2*math.Pi - loop.TurningAngle()

// TODO(roberts): The error bound below is much larger than it should be.
if math.Abs(area-gaussArea) > 1e-9 {
t.Errorf("%d. %v.Area() = %v want %v", x, loop, area, gaussArea)
}
}
}

func TestLoopGetAreaConsistentWithSign(t *testing.T) {
// TODO(roberts): Uncomment when Loop has IsValid
/*
// Test that Area() returns an area near 0 for degenerate loops that
// contain almost no points, and an area near 4*pi for degenerate loops that
// contain almost all points.
const maxVertices = 6

for i := 0; i < 50; i++ {
numVertices := 3 + randomUniformInt(maxVertices-3+1)
// Repeatedly choose N vertices that are exactly on the equator until we
// find some that form a valid loop.
var loop = Loop{}
for !loop.IsValid() {
var vertices []Point
for i := 0; i < numVertices; i++ {
// We limit longitude to the range [0, 90] to ensure that the loop is
// degenerate (as opposed to following the entire equator).
vertices = append(vertices,
PointFromLatLng(LatLng{0, s1.Angle(randomFloat64()) * math.Pi / 2}))
}
loop.vertices = vertices
break
}

ccw := loop.IsNormalized()
want := 0.0
if !ccw {
want = 4 * math.Pi
}

// TODO(roberts): The error bound below is much larger than it should be.
if got := loop.Area(); !float64Near(got, want, 1e-8) {
t.Errorf("%v.Area() = %v, want %v", loop, got, want)
}
p := PointFromCoords(0, 0, 1)
if got := loop.ContainsPoint(p); got != !ccw {
t.Errorf("%v.ContainsPoint(%v) = %v, want %v", got, p, !ccw)
}
}
*/
}

func TestLoopNormalizedCompatibleWithContains(t *testing.T) {
p := parsePoint("40:40")

tests := []*Loop{
lineTriangle,
skinnyChevron,
}

// Checks that if a loop is normalized, it doesn't contain a
// point outside of it, and vice versa.
for _, loop := range tests {
flip := cloneLoop(loop)

flip.Invert()
if norm, contains := loop.IsNormalized(), loop.ContainsPoint(p); norm == contains {
t.Errorf("loop.IsNormalized() = %b == loop.ContainsPoint(%v) = %v, want !=", norm, p, contains)
}
if norm, contains := flip.IsNormalized(), flip.ContainsPoint(p); norm == contains {
t.Errorf("flip.IsNormalized() = %b == flip.ContainsPoint(%v) = %v, want !=", norm, p, contains)
}
if a, b := loop.IsNormalized(), flip.IsNormalized(); a == b {
t.Errorf("a loop and it's invert can not both be normalized")
}
flip.Normalize()
if flip.ContainsPoint(p) {
t.Errorf("%v.ContainsPoint(%v) = true, want false", flip, p)
}
}
}

func TestLoopIsValidDetectsInvalidLoops(t *testing.T) {
tests := []struct {
msg    string
points []Point
}{
// Not enough vertices. Note that all single-vertex loops are valid; they
// are interpreted as being either "empty" or "full".
{
msg:    "loop has no vertices",
points: parsePoints(""),
},
{
msg:    "loop has too few vertices",
points: parsePoints("20:20, 21:21"),
},
// degenerate edge checks happen in validation before duplicate vertices.
{
msg:    "loop has degenerate first edge",
points: parsePoints("20:20, 20:20, 20:21"),
},
{
msg:    "loop has degenerate third edge",
points: parsePoints("20:20, 20:21, 20:20"),
},
// TODO(roberts): Uncomment these cases when FindAnyCrossings is in.
/*
{
msg:    "loop has duplicate points",
points: parsePoints("20:20, 21:21, 21:20, 20:20, 20:21"),
},
{
msg:    "loop has crossing edges",
points: parsePoints("20:20, 21:21, 21:20.5, 21:20, 20:21"),
},
*/
{
// Ensure points are not normalized.
msg: "loop with non-normalized vertices",
points: []Point{
Point{r3.Vector{2, 0, 0}},
Point{r3.Vector{0, 1, 0}},
Point{r3.Vector{0, 0, 1}},
},
},
{
// Adjacent antipodal vertices
msg: "loop with antipodal points",
points: []Point{
Point{r3.Vector{1, 0, 0}},
Point{r3.Vector{-1, 0, 0}},
Point{r3.Vector{0, 0, 1}},
},
},
}

for _, test := range tests {
loop := LoopFromPoints(test.points)
err := loop.findValidationError()
if err == nil {
t.Errorf("%s. %v.findValidationError() = nil, want err to be non-nil", test.msg, loop)
}
// The C++ tests also tests that the returned error message string contains
// a specific set of text. That part of the test is skipped here.
}
}

const (
// TODO(roberts): Convert these into changeable flags or parameters.
// A loop with a 10km radius and 4096 vertices has an edge length of 15 meters.
defaultRadiusKm   = 10.0
numLoopSamples    = 16
numQueriesPerLoop = 100
)

func BenchmarkLoopContainsPoint(b *testing.B) {
// Benchmark ContainsPoint() on regular loops. The query points for a loop are
// chosen so that they all lie in the loop's bounding rectangle (to avoid the
// quick-rejection code path).

// C++ ranges from 4 -> 256k by powers of 2 for number of vertices for benchmarking.
vertices := 4
for n := 1; n <= 17; n++ {
b.Run(fmt.Sprintf("%d", vertices),
func(b *testing.B) {
b.StopTimer()
loops := make([]*Loop, numLoopSamples)
for i := 0; i < numLoopSamples; i++ {
loops[i] = RegularLoop(randomPoint(), kmToAngle(10.0), vertices)
}

queries := make([][]Point, numLoopSamples)

for i, loop := range loops {
queries[i] = make([]Point, numQueriesPerLoop)
for j := 0; j < numQueriesPerLoop; j++ {
queries[i][j] = samplePointFromRect(loop.RectBound())
}
}

b.StartTimer()
for i := 0; i < b.N; i++ {
loops[i%numLoopSamples].ContainsPoint(queries[i%numLoopSamples][i%numQueriesPerLoop])
}
})
vertices *= 2
}
}
