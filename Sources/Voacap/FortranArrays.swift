import Foundation

// One-based, column-major arrays that mirror the Fortran originals, so the
// translated routines can keep VOACAP's own indices.

struct A1 {
    var d: [Float]
    init(_ n: Int, _ value: Float = 0) { d = Array(repeating: value, count: n) }
    init(_ values: [Float]) { d = values }
    @inline(__always) subscript(i: Int) -> Float {
        get { d[i - 1] }
        set { d[i - 1] = newValue }
    }
    mutating func fill(_ value: Float) { for i in d.indices { d[i] = value } }
}

struct A2 {
    var d: [Float]
    let n1: Int
    init(_ n1: Int, _ n2: Int, _ value: Float = 0) {
        self.n1 = n1
        d = Array(repeating: value, count: n1 * n2)
    }
    @inline(__always) subscript(i: Int, j: Int) -> Float {
        get { d[(j - 1) * n1 + i - 1] }
        set { d[(j - 1) * n1 + i - 1] = newValue }
    }
    mutating func fill(_ value: Float) { for i in d.indices { d[i] = value } }
}

struct A3 {
    var d: [Float]
    let n1: Int
    let n2: Int
    init(_ n1: Int, _ n2: Int, _ n3: Int, _ value: Float = 0) {
        self.n1 = n1
        self.n2 = n2
        d = Array(repeating: value, count: n1 * n2 * n3)
    }
    @inline(__always) subscript(i: Int, j: Int, k: Int) -> Float {
        get { d[((k - 1) * n2 + j - 1) * n1 + i - 1] }
        set { d[((k - 1) * n2 + j - 1) * n1 + i - 1] = newValue }
    }
    mutating func fill(_ value: Float) { for i in d.indices { d[i] = value } }
}

struct I1 {
    var d: [Int]
    init(_ n: Int, _ value: Int = 0) { d = Array(repeating: value, count: n) }
    init(_ values: [Int]) { d = values }
    @inline(__always) subscript(i: Int) -> Int {
        get { d[i - 1] }
        set { d[i - 1] = newValue }
    }
    mutating func fill(_ value: Int) { for i in d.indices { d[i] = value } }
}

struct I2 {
    var d: [Int]
    let n1: Int
    init(_ n1: Int, _ n2: Int, _ value: Int = 0) {
        self.n1 = n1
        d = Array(repeating: value, count: n1 * n2)
    }
    @inline(__always) subscript(i: Int, j: Int) -> Int {
        get { d[(j - 1) * n1 + i - 1] }
        set { d[(j - 1) * n1 + i - 1] = newValue }
    }
    mutating func fill(_ value: Int) { for i in d.indices { d[i] = value } }
}

struct I3 {
    var d: [Int]
    let n1: Int
    let n2: Int
    init(_ n1: Int, _ n2: Int, _ n3: Int, _ value: Int = 0) {
        self.n1 = n1
        self.n2 = n2
        d = Array(repeating: value, count: n1 * n2 * n3)
    }
    @inline(__always) subscript(i: Int, j: Int, k: Int) -> Int {
        get { d[((k - 1) * n2 + j - 1) * n1 + i - 1] }
        set { d[((k - 1) * n2 + j - 1) * n1 + i - 1] = newValue }
    }
    mutating func fill(_ value: Int) { for i in d.indices { d[i] = value } }
}

// Fortran intrinsics with the same rounding behaviour.

/// SIGN(a, b): the magnitude of `a` with the sign of `b`.
@inline(__always) func fsign(_ a: Float, _ b: Float) -> Float { b >= 0 ? abs(a) : -abs(a) }
/// INT(): truncation toward zero.
@inline(__always) func fint(_ x: Float) -> Int { Int(x) }
/// NINT(): round half away from zero.
@inline(__always) func nint(_ x: Float) -> Int { Int(x.rounded(.toNearestOrAwayFromZero)) }
@inline(__always) func alog(_ x: Float) -> Float { logf(x) }
@inline(__always) func alog10(_ x: Float) -> Float { log10f(x) }
@inline(__always) func fsin(_ x: Float) -> Float { sinf(x) }
@inline(__always) func fcos(_ x: Float) -> Float { cosf(x) }
@inline(__always) func ftan(_ x: Float) -> Float { tanf(x) }
@inline(__always) func fasin(_ x: Float) -> Float { asinf(x) }
@inline(__always) func facos(_ x: Float) -> Float { acosf(x) }
@inline(__always) func fatan(_ x: Float) -> Float { atanf(x) }
@inline(__always) func fatan2(_ y: Float, _ x: Float) -> Float { atan2f(y, x) }
@inline(__always) func fexp(_ x: Float) -> Float { expf(x) }
@inline(__always) func fsqrt(_ x: Float) -> Float { x.squareRoot() }
/// AMOD(a, p): remainder with the sign of `a`.
@inline(__always) func amod(_ a: Float, _ p: Float) -> Float { fmodf(a, p) }
