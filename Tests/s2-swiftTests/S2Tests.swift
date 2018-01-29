//
//  S2Tests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


let	epsilon = 1e-14

// float64Eq reports whether the two values are within the default epsilon.
func float64Eq(_ x: Double, _ y: Double) -> Bool {
  return float64Near(x, y, epsilon: epsilon)
}

// float64Near reports whether the two values are within the given epsilon.
func float64Near(_ x: Double, _ y: Double, epsilon: Double) -> Bool {
	return abs(x - y) <= epsilon
}

// TODO: Add in flag to allow specifying the random seed for repeatable tests.

// kmToAngle converts a distance on the Earth's surface to an angle.
func kmToAngle(km: Double) -> Double {
	// The Earth's mean radius in kilometers (according to NASA).
	let earthRadiusKm = 6371.01
  // ???
	return km / earthRadiusKm
}

func llDegrees(_ lat: Double, _ lng: Double) -> LatLng {
  return LatLng(lat: lat * toRadians, lng: lng * toRadians)
}

func p(_ x: Double, _ y: Double, _ z: Double) -> S2Point {
  return S2Point(x: x, y: y, z: z)
}

func p(_ lat: Double, _ lng: Double) -> S2Point {
  return llDegrees(lat, lng).toPoint()
}

// randomBits returns a 64-bit random unsigned integer whose lowest "num" are random, and
// whose other bits are zero.
func randomBits(num: UInt32) -> UInt64 {
	// Make sure the request is for not more than 63 bits.
	let num = min(num, 63)
  let r = randomUInt64()
	return r & ((UInt64(1) << UInt64(num)) - 1)
}

// Return a uniformly distributed 64-bit unsigned integer.
func randomUInt64() -> UInt64 {
	return UInt64(arc4random()) | (UInt64(arc4random()) << 32)
}

// Return a uniformly distributed 32-bit unsigned integer.
func randomUInt32() -> UInt32 {
  return UInt32(randomBits(num: 32))
}

// Return a uniformly distributed 32-bit unsigned integer.
func randomInt(hi: Int) -> Int {
  return Int(UInt32(randomBits(num: 32)))
}

// randomFloat64 returns a uniformly distributed value in the range [0,1).
// Note that the values returned are all multiples of 2**-53, which means that
// not all possible values in this range are returned.
func randomFloat64() -> Double {
  let randomFloatBits = UInt32(53)
	return ldexp(Double(randomBits(num: randomFloatBits)), -Int(randomFloatBits))
}

// randomUniformInt returns a uniformly distributed integer in the range [0,n).
// NOTE: This is replicated here to stay in sync with how the C++ code generates
// uniform randoms. (instead of using Go's math/rand package directly).
func randomUniformInt(n: Int) -> Int {
	return Int(randomFloat64() * Double(n))
}

// randomUniformFloat64 returns a uniformly distributed value in the range [min, max).
func randomUniformFloat64(min: Double, max: Double) -> Double {
	return min + randomFloat64() * (max - min)
}

// oneIn returns true with a probability of 1/n.
func oneIn(n: Int) -> Bool {
	return randomUniformInt(n: n) == 0
}

// randomPoint returns a random unit-length vector.
func randomPoint() -> S2Point {
  let x = randomUniformFloat64(min: -1, max: 1)
  let y = randomUniformFloat64(min: -1, max: 1)
  let z = randomUniformFloat64(min: -1, max: 1)
	return S2Point(x: x, y: y, z: z)
}

// randomFrame returns a right-handed coordinate frame (three orthonormal vectors) for
// a randomly generated point.
func randomFrame() -> R3Matrix {
	return randomFrameAtPoint(z: randomPoint())
}

// randomFrameAtPoint returns a right-handed coordinate frame using the given
// point as the z-axis. The x- and y-axes are computed such that (x,y,z) is a
// right-handed coordinate frame (three orthonormal vectors).
func randomFrameAtPoint(z: S2Point) -> R3Matrix {
  let x = z.v.cross(randomPoint().v).s2
	let y = z.v.cross(x.v).s2
  //
  return R3Matrix(
    r1: [x.x, y.x, z.x],
    r2: [x.y, y.y, z.y],
    r3: [x.z, y.z, z.z])
}

// randomCellIDForLevel returns a random CellID at the given level.
// The distribution is uniform over the space of cell ids, but only
// approximately uniform over the surface of the sphere.
func randomCellIdForLevel(level: Int) -> CellId {
	let face = randomUniformInt(n: CellId.numFaces)
	let pos = randomUInt64() & UInt64((1 << CellId.posBits)-1)
	return CellId(face: face, pos: pos, level: level)
}

// randomCellID returns a random CellID at a randomly chosen
// level. The distribution is uniform over the space of cell ids,
// but only approximately uniform over the surface of the sphere.
func randomCellId() -> CellId {
	return randomCellIdForLevel(level: randomUniformInt(n: CellId.maxLevel + 1))
}

// parsePoint returns an Point from the latitude-longitude coordinate in degrees
// in the given string, or the origin if the string was invalid.
// e.g., "-20:150"
func parsePoint(_ s: String) -> S2Point {
	let p = parsePoints(s)
	if p.count > 0 {
		return p[0]
	}
  return S2Point(x: 0, y: 0, z: 0)
}

// parseRect returns the minimal bounding Rect that contains the one or more
// latitude-longitude coordinates in degrees in the given string.
// Examples of input:
//   "-20:150"                     // one point
//   "-20:150, -20:151, -19:150"   // three points
func parseRect(_ s: String) -> S2Rect {
	var rect = S2Rect.empty
	let lls = parseLatLngs(s)
//	if lls.count > 0 {
//		rect = S2Rect((latLng: lls[0])
//	}
//
	for ll in lls {
		rect = rect.add(ll)
	}
	return rect
}

// parseLatLngs splits up a string of lat:lng points and returns the list of parsed
// entries.
func parseLatLngs(_ s: String) -> [LatLng] {
	let pieces = s.components(separatedBy: ",")
	var lls = [LatLng]()
	for piece in pieces {
    // get a trimmed non-empty string
    let piece = piece.trimmingCharacters(in: NSCharacterSet.whitespaces)
		if piece == "" {
			continue
		}
    let p = piece.components(separatedBy: ":")
		if p.count != 2 {
			fatalError("invalid input string for parseLatLngs")
		}
		guard let lat = Double(p[0]) else {
			fatalError("invalid float in parseLatLngs")
		}
		guard let lng = Double(p[1]) else {
			fatalError("invalid float in parseLatLngs")
		}

		lls.append(llDegrees(lat, lng))
	}
	return lls
}

// parsePoints takes a string of lat:lng points and returns the set of Points it defines.
func parsePoints(_ s: String) -> [S2Point] {
	let lls = parseLatLngs(s)
	var points = [S2Point]()
	for ll in lls {
		points.append(ll.toPoint())
	}
	return points
}

// skewedInt returns a number in the range [0,2^max_log-1] with bias towards smaller numbers.
func skewedInt(maxLog: Int) -> Int {
	let base = Int32(arc4random() & 0x7fffffff) % (Int32(maxLog + 1))
	return Int(randomBits(num: 31)) & Int((1 << base) - 1)
}

// randomCap returns a cap with a random axis such that the log of its area is
// uniformly distributed between the logs of the two given values. The log of
// the cap angle is also approximately uniformly distributed.
func randomCap(minArea: Double, maxArea: Double) -> S2Cap {
  let capArea = maxArea * pow(minArea/maxArea, randomFloat64())
  return S2Cap(center: randomPoint(), area: capArea)
}

func random(lo: Double, hi: Double) -> Double {
  let r = randomFloat64()
  return lo + r * (hi - lo)
}

func randomRect() -> S2Rect {
  let lat0 = random(lo: -.pi / 2.0, hi: .pi / 2.0)
  let lat1 = random(lo: -.pi / 2.0, hi: .pi / 2.0)
  let lng0 = random(lo: -.pi, hi: .pi)
  let lng1 = random(lo: -.pi, hi: .pi)
  let latInterval = R1Interval(p0: lat0, p1: lat1)
  let lngInterval = S1Interval(p0: lng0, p1: lng1)
  return S2Rect(lat: latInterval, lng: lngInterval)
}

func randomPolyline() -> S2Polyline {
  let n = randomUInt32() % 10
  var latLngs: [LatLng] = []
  for _ in 0..<n {
    let lat0 = random(lo: -.pi / 2.0, hi: .pi / 2.0)
    let lng0 = random(lo: -.pi, hi: .pi)
    let latLng = LatLng(lat: lat0, lng: lng0)
    latLngs.append(latLng)
  }
  return S2Polyline(latLngs: latLngs)
}

// pointsApproxEquals reports whether the two points are within the given distance
// of each other. This is the same as Point.ApproxEquals but permits specifying
// the epsilon.
func pointsApproxEquals(a: S2Point, b:S2Point, epsilon: Double) -> Bool {
	return Double(a.angle(b)) <= epsilon
}

let	rectErrorLat = 10 * Cell.dblEpsilon
let rectErrorLng = Cell.dblEpsilon

// r2PointsApproxEqual reports whether the two points are within the given epsilon.
func r2PointsApproxEquals(a: S2Point, b: S2Point, epsilon: Double) -> Bool {
	return float64Near(a.x, b.x, epsilon: epsilon) && float64Near(a.y, b.y, epsilon: epsilon)
}

// rectsApproxEqual reports whether the two rect are within the given tolerances
// at each corner from each other. The tolerances are specific to each axis.
func rectsApproxEqual(a: S2Rect, b: S2Rect, tolLat: Double, tolLng: Double) -> Bool {
	return abs(a.lat.lo-b.lat.lo) < tolLat &&
		abs(a.lat.hi-b.lat.hi) < tolLat &&
		abs(a.lng.lo-b.lng.lo) < tolLng &&
		abs(a.lng.hi-b.lng.hi) < tolLng
}

// matricesApproxEqual reports whether all cells in both matrices are equal within
// the default floating point epsilon.
func matricesApproxEqual(m1: R3Matrix, m2: R3Matrix) -> Bool {
	return float64Eq(m1[0, 0], m2[0, 0]) &&
		float64Eq(m1[0, 1], m2[0, 1]) &&
		float64Eq(m1[0, 2], m2[0, 2]) &&

		float64Eq(m1[1, 0], m2[1, 0]) &&
		float64Eq(m1[1, 1], m2[1, 1]) &&
		float64Eq(m1[1, 2], m2[1, 2]) &&

		float64Eq(m1[2, 0], m2[2, 0]) &&
		float64Eq(m1[2, 1], m2[2, 1]) &&
		float64Eq(m1[2, 2], m2[2, 2])
}

// samplePointFromRect returns a point chosen uniformly at random (with respect
// to area on the sphere) from the given rectangle.
func samplePointFromRect(rect: S2Rect) -> S2Point {
	// First choose a latitude uniformly with respect to area on the sphere.
  let sinLo = sin(rect.lat.lo)
	let sinHi = sin(rect.lat.hi)
	let lat = asin(randomUniformFloat64(min: sinLo, max: sinHi))

	// Now choose longitude uniformly within the given range.
	let lng = rect.lng.lo + randomFloat64()*rect.lng.length

  return LatLng(lat: lat, lng: lng).toPoint()
}

// samplePointFromCap returns a point chosen uniformly at random (with respect
// to area) from the given cap.
func samplePointFromCap(c: S2Cap) -> S2Point {
	// We consider the cap axis to be the "z" axis. We choose two other axes to
	// complete the coordinate frame.
	let m = S2Point.getFrame(c.center)

	// The surface area of a spherical cap is directly proportional to its
	// height. First we choose a random height, and then we choose a random
	// point along the circle at that height.
	let h = randomFloat64() * c.height
	let theta = 2 * .pi * randomFloat64()
	let r = sqrt(h * (2 - h))

	// The result should already be very close to unit-length, but we might as
	// well make it accurate as possible.
  let p = S2Point(x: cos(theta) * r, y: sin(theta) * r, z: 1 - h)
  return R3Matrix.fromFrame(m, point: p)
}

// perturbATowardsB returns a point that has been shifted some distance towards the
// second point based on a random number.
func perturbATowardsB(a: S2Point, b: S2Point) -> S2Point {
	let choice = randomFloat64()
	if choice < 0.1 {
		return a
	}
	if choice < 0.3 {
		// Return a point that is exactly proportional to A and that still
		// satisfies IsUnitlength.
		while true {
      let b3 = (randomFloat64() - 0.5) * Cell.dblEpsilon
      let b2 = 2 - a.v.norm + 5 * b3
      let b = a.v.mul(b2)
			if !b.approxEquals(a.v) && b.isUnit {
        return b.s2
			}
		}
	}
	if choice < 0.5 {
		// Return a point such that the distance squared to A will underflow.
    return interpolateAtDistance(ax: 1e-300, a: a, b: b)
	}
	// Otherwise return a point whose distance from A is near dblEpsilon such
	// that the log of the pdf is uniformly distributed.
	let distance = Cell.dblEpsilon * 1e-5 * pow(1e6, randomFloat64())
  return interpolateAtDistance(ax: distance, a: a, b: b)
}

// perturbedCornerOrMidpoint returns a Point from a line segment whose endpoints are
// difficult to handle correctly. Given two adjacent cube vertices P and Q,
// it returns either an edge midpoint, face midpoint, or corner vertex that is
// in the plane of PQ and that has been perturbed slightly. It also sometimes
// returns a random point from anywhere on the sphere.
func perturbedCornerOrMidpoint(p: S2Point, q: S2Point) -> S2Point {
	var a = p.v.mul(Double(randomUniformInt(n: 3) - 1)).add(q.v.mul(Double(randomUniformInt(n: 3) - 1)))
	if oneIn(n: 10) {
		// This perturbation often has no effect except on coordinates that are
		// zero, in which case the perturbed value is so small that operations on
		// it often result in underflow.
		a = a.add(randomPoint().v.mul(pow(1e-300, randomFloat64())))
	} else if oneIn(n: 2) {
		// For coordinates near 1 (say > 0.5), this perturbation yields values
		// that are only a few representable values away from the initial value.
		a = a.add(randomPoint().v.mul(4 * Cell.dblEpsilon))
	} else {
		// A perturbation whose magnitude is in the range [1e-25, 1e-10].
		a = a.add(randomPoint().v.mul(1e-10 * pow(1e-15, randomFloat64())))
	}

	if a.norm2 < Cell.dblEpsilon {
		// If a.norm2 is denormalized, normalized() loses too much precision.
		return perturbedCornerOrMidpoint(p: p, q: q)
	}
  return a.s2
}

/// A Shape representing an arbitrary set of edges. It
/// is used for testing, but it can also be useful if you have, say, a
/// collection of polylines and don't care about memory efficiency (since
/// this type would store most of the vertices twice).
struct EdgeVectorShape {
  var edges: [Edge]
}

extension EdgeVectorShape: Shape {

  /// Creates an edgeVectorShape of length 1 from the given points.
  init(a: S2Point, b: S2Point) {
    edges = [Edge(v0: a, v1: b)]
  }

  mutating func add(a: S2Point, b: S2Point) {
    edges.append(Edge(v0: a, v1: b))
  }
  
  func numEdges() -> Int { return edges.count }
  func edge(_ id: Int) -> Edge { return edges[id] }
  func hasInterior() -> Bool { return false }
  func referencePoint() -> ReferencePoint { return ReferencePoint(origin: true, contained: false) }
  func numChains() -> Int { return edges.count }
  func chain(_ chainId: Int) -> Chain { return Chain(start: chainId, length: 1) }
  func chainEdge(chainId: Int, offset: Int) -> Edge { return edges[chainId] }
  func chainPosition(_ edgeId: Int) -> ChainPosition { return ChainPosition(chainId: edgeId, offset: 0) }
  func dimension() -> ShapeDimension { return .polylineGeometry }
  func containsOrigin() -> Bool { return false }

}

/// Constructs a polygon from the sequence of loops in the input
/// string. Loops are automatically normalized by inverting them if necessary
/// so that they enclose at most half of the unit sphere. (Historically this was
/// once a requirement of polygon loops. It also hides the problem that if the
/// user thinks of the coordinates as X:Y rather than LAT:LNG, it yields a loop
/// with the opposite orientation.)
///
/// Loops are semicolon separated in the input string with each loop using the
/// same format as above.
///
/// Examples of the input format:
///     "10:20, 90:0, 20:30"                                  // one loop
///     "10:20, 90:0, 20:30; 5.5:6.5, -90:-180, -15.2:20.3"   // two loops
///     ""       // the empty polygon (consisting of no loops)
///     "empty"  // the empty polygon (consisting of no loops)
///     "full"   // the full polygon (consisting of one full loop)
func makePolygon(_ s: String, normalize: Bool) -> S2Polygon {
  let strs = (s == "empty") ? [] : s.components(separatedBy: ";")
  var loops: [S2Loop] = []
  for str in strs {
    if str == "" {
      continue
    }
    let loop = makeLoop(str.trimmingCharacters(in: .whitespaces))
    if normalize && !loop.isFull {
      // TODO(roberts): Uncomment once Normalize is implemented.
//      loop.normalize()
    }
    loops.append(loop)
  }
  return S2Polygon(loops: loops)
}

/// Constructs a Polyline from the given string.
func makePolyline(_ s: String) -> S2Polyline {
  let points = parsePoints(s)
  return S2Polyline(points: points)
}

// concentricLoopsPolygon constructs a polygon with the specified center as a
// number of concentric loops and vertices per loop.
func concentricLoopsPolygon(center: S2Point, numLoops: Int, verticesPerLoop: Int) -> S2Polygon {
  var loops: [S2Loop] = []
  for li in 0..<numLoops {
    let radius = S1Angle(0.005 * Double(li + 1) / Double(numLoops))
    loops.append(S2Loop.regularLoop(center: center, radius: radius, numVertices: verticesPerLoop))
  }
  return S2Polygon(loops: loops)
}

