import Foundation

// VOACAP's own horizontal dipole (antenna type 23, the ITSA-1 model in
// IONGAIN) and the sine/cosine integrals it needs (CISI).

extension VoacapAntenna {
    /// CISI: cosine and sine integrals.
    private static func cisi(_ x: Float) -> (ci: Float, si: Float) {
        if x - 1 < 0 {
            let sq = x * x
            var ci = 0.5772156649 + alog(x)
            var term = -1 * sq / 4
            var g: Float = 4
            repeat {
                ci = ci + term
                term = -1 * term * sq * (g - 2) / ((g - 1) * g * g)
                g = g + 2
            } while abs(term) - 0.00005 > 0
            var si: Float = 0
            term = x
            g = 3
            repeat {
                si = si + term
                term = -1 * term * sq * (g - 2) / ((g - 1) * g * g)
                g = g + 2
            } while abs(term) - 0.00005 > 0
            return (ci, si)
        }
        let x2 = x * x
        let t = fcos(x)
        let s = fsin(x)
        var g = ((((x2 + 48.196927) * x2 + 482.485984) * x2 + 1114.978885) * x2 + 449.690326) * x2
        g = ((((x2 + 42.252855) * x2 + 302.757865) * x2 + 352.018498) * x2 + 21.821899) / g
        var f = ((((x2 + 40.021433) * x2 + 322.624911) * x2 + 570.236280) * x2 + 157.105423) * x
        f = ((((x2 + 38.027264) * x2 + 265.187033) * x2 + 335.677320) * x2 + 38.102495) / f
        return (f * s - g * t, 1.5708 - f * t - g * s)
    }

    /// IONGAIN2 for the horizontal dipole. `length` and `height` are metres
    /// when positive and wavelengths when negative, as in VOACAP's antenna files.
    private static func dipoleGain(offAzimuth beta: Float, sigma: Float, er: Float, length el: Float,
                                   height h: Float, extraGain: Float, delta: Float, fmc: Float) -> Float {
        let floor: Float = -30
        if delta <= 0 { return floor }
        var x = 18000 * sigma / fmc
        let t = fcos(delta)
        let q = fsin(delta)
        let r = q * q
        let s = r * r
        let ert = er - t * t
        let rho = fsqrt(ert * ert + x * x)
        let rho12 = fsqrt(rho)
        let alpha = -fatan(x / ert)
        let u = er * er + x * x
        let v = fsqrt(u)
        let asxv = fasin(x / v)
        var cv = fsqrt(rho * rho + u * u * s - 2 * rho * u * r * fcos(alpha + 2 * asxv))
            / (rho + u * r + 2 * rho12 * v * q * fcos(alpha * 0.5 + asxv))
        let a = 2 * rho12 * q * v * fsin(alpha * 0.5 + asxv)
        let wave = 299.7925 / fmc
        let b = rho - u * r
        var psiv: Float
        if b < 0 {
            psiv = fatan(a / b) + 3.141593
        } else if b == 0 {
            psiv = a < 0 ? -1.570796 : (a == 0 ? 0 : 1.570796)
        } else {
            psiv = fatan(a / b)
        }
        let reta = beta * 0.017453293
        let psih = fatan(2 * rho12 * q * fsin(alpha * 0.5) / (rho - r))
        let sb = fsin(reta)
        let cb = fcos(reta)
        x = 1
        let ch = fsqrt(rho * rho + s - 2 * rho * r * fcos(alpha)) / (rho + r + 2 * rho12 * q * fcos(alpha * 0.5))
        if abs(cv - ch) <= 0.001 { cv = ch }
        let fac: Float = el < 0 ? 3.1416 * abs(el) : 3.1416 * el / wave
        let hwave: Float = h > 0 ? 6.2832 * h / wave : 6.2832 * abs(h)
        let four = cisi(fac * 4)
        let two = cisi(fac * 2)
        let g = 1 / ftan(fac)
        let ss = sb * sb
        let c = t * t
        let cinKl = 0.577 + alog(2 * fac) - two.ci
        let cin2Kl = 0.577 + alog(4 * fac) - four.ci
        let xR = 30 * ((1 - g * g) * cin2Kl + 4 * g * g * cinKl + 2 * g * (four.si - 2 * two.si))
        let xKv = cv * cv + 1 - 2 * cv * fcos(psiv - 2 * hwave * q)
        let xKh = ch * ch + 1 - 2 * ch * fcos(psih - 2 * hwave * q)
        let xt1 = fcos(fac * sb * t) - fcos(fac)
        let xb1 = 1 - ss * c
        let xt2 = xKv * ss * r + xKh * cb * cb
        var rain = 120 / xR * (xt1 / xb1) * (xt1 / xb1) * xt2
        rain = max(rain, 0.00001)
        rain = 10 * alog10(rain)
        rain = rain + extraGain
        if rain <= floor { rain = floor }
        return rain
    }

    /// VOACAP's type 23 horizontal dipole over real ground.
    ///
    /// - Parameters:
    ///   - offAzimuth: degrees between the path bearing and the dipole's broadside direction.
    ///   - length: metres when positive, wavelengths when negative (−0.5 is a half-wave dipole).
    ///   - height: metres when positive, wavelengths when negative.
    public static func horizontalDipole(offAzimuth: Float, length: Float = -0.5, height: Float = -0.25,
                                 dielectric: Float = 4, conductivity: Float = 0.001,
                                 extraGain: Float = 0) -> VoacapAntenna {
        var antenna = VoacapAntenna.isotropic(gain: 0)
        for mhz in 2...30 {
            let frequency = Float(mhz)
            // VOACAP blends toward the overhead gain so 90° is the same at every azimuth.
            let overhead = dipoleGain(offAzimuth: 0, sigma: conductivity, er: dielectric, length: length,
                                      height: height, extraGain: extraGain, delta: 1.570796, fmc: frequency)
            for elevation in 0...90 {
                let delta = Float(elevation) * 0.01745329
                var rain = dipoleGain(offAzimuth: offAzimuth, sigma: conductivity, er: dielectric, length: length,
                                      height: height, extraGain: extraGain, delta: delta, fmc: frequency)
                let sine = fsin(delta)
                let gmorph = sine * sine * sine * sine
                rain = overhead * gmorph + rain * (1 - gmorph)
                antenna.gain[elevation * 30 + mhz - 1] = rain
            }
        }
        return antenna
    }
}
