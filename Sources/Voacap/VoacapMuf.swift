import Foundation

// Electron density profile, virtual heights and the circuit MUF:
// LECDEN, GETHP, BENDY, PEN, GENION, F2DIS and CURMUF.

extension VoacapEngine {
    /// The piecewise-linear table search GETHP uses in both directions.
    @inline(__always)
    private func lookup(_ x: Float, _ xs: UnsafePointer<Float>, _ ys: UnsafePointer<Float>) -> Float {
        // xs and ys are 50-element columns, zero-based here.
        for ih in 0..<49 {
            let a = xs[ih], b = xs[ih + 1]
            if a == b {
                if x == a { return ys[ih] }
                continue
            } else if a < b {
                if a - x <= 0 && !(x - b < 0) { continue }
            } else {
                if a - x < 0 { continue }
                if x - b <= 0 { continue }
            }
            return ys[ih] + (x - a) * (ys[ih + 1] - ys[ih]) / (b - a)
        }
        return ys[49]
    }

    /// GETHP: true and virtual height of reflection at a vertical frequency,
    /// by Gaussian integration through the electron density profile.
    func gethp(_ fxx: Float) -> (hp: Float, ht: Float) {
        let fr = fxx * fxx
        return htr.d.withUnsafeBufferPointer { htrBuffer in
            fnsq.d.withUnsafeBufferPointer { fnsqBuffer in
                let h = htrBuffer.baseAddress!, f = fnsqBuffer.baseAddress!
                if fr - f[0] <= 0 { return (h[0], h[0]) }
                var ht = lookup(fr, f, h)
                var hp: Float = 0
                ht = ht - h[0]
                let hrmz = ht * twdiv * xnpl
                for ig in 1...nt {
                    var zg = 1.0 - xt[ig]
                    var zi = ht * (1 - twdiv * zg) + h[0]
                    var ysq = zi - h[0] <= 0 ? f[0] : lookup(zi, h, f)
                    ysq = ysq / fr
                    ysq = min(ysq, 0.9999)
                    let ymup = 1 / fsqrt(1 - ysq)
                    zg = 1 + xt[ig]
                    zi = ht * (1 - twdiv * zg) + h[0]
                    ysq = zi - h[0] <= 0 ? f[0] : lookup(zi, h, f)
                    ysq = ysq / fr
                    ysq = min(ysq, 0.9999)
                    let zmup = 1 / fsqrt(1 - ysq)
                    hp = hp + wt[ig] * (ymup + zmup)
                }
                return (h[0] + hrmz * hp, h[0] + ht)
            }
        }
    }

    /// BENDY: bending in a parabolic layer.
    func bendy(_ i: Int, _ k: Int, _ f: Float) -> Float {
        var x = f / fi[i, k]
        x = min(x, 0.999)
        return 0.5 * x * yi[i, k] * alog((1 + x) / (1 - x))
    }

    /// PEN: retardation in a parabolic layer.
    func pen(_ i: Int, _ k: Int, _ f: Float) -> Float {
        var x = f / fi[i, k]
        x = max(x, 1.001)
        return yi[i, k] * alog((1 + x) / (x - 1)) * x
    }

    /// LECDEN: build the electron density profile for control point `k`.
    func lecden(_ k: Int) {
        let xlow: Float = 0.8516
        let hz = hi[1, k] - yi[1, k]
        let xup = 0.98 * fi[1, k] / fi[3, k]
        let hlow = hz + yi[1, k] * (1 + fsqrt(1 - xlow * xlow))
        let hte = hi[1, k] + yi[1, k]
        var hb2 = hi[3, k] - yi[3, k]
        let fce = fi[1, k] * fi[1, k]
        let fc2 = fi[3, k] * fi[3, k]
        let hup = hb2 + yi[3, k] * (1 - fsqrt(1 - xup * xup))
        let fup = xup * xup * fc2
        let flow = xlow * xlow * fce
        var asp: Float = 0
        if hup - hlow > 0 { asp = (fup - flow) / (hup - hlow) }
        var fc1: Float = 0, hb1: Float = 0, ht1: Float = 0, s1: Float = 0
        var lin = 0
        if fi[2, k] > 0 {
            fc1 = fi[2, k] * fi[2, k]
            hb1 = hi[2, k] - yi[2, k]
            ht1 = hi[2, k] + yi[2, k]
            // Height of F2 at the F1 critical frequency.
            let htw = hb2 + yi[3, k] * (1 - fsqrt(1 - fc1 / fc2))
            var parabolic = false
            if htw > hi[2, k] + 0.001 {
                parabolic = true
            } else {
                var ys = htw - hb1
                ys = max(1, ys)
                s1 = fc1 / ys
                hi[2, k] = htw
                yi[2, k] = ys
                if hb2 - hb1 < 0 {
                    yi[2, k] = hi[2, k] - hb2
                    parabolic = true
                } else {
                    lin = 1
                    fsecv[k] = fsqrt(flow)
                    var denom = 1 - (fi[1, k] / fi[2, k]) * (fi[1, k] / fi[2, k])
                    denom = max(0.17, denom)
                    let yb = (hi[2, k] - hi[1, k]) / denom
                    if ys - yb >= 0 {
                        ys = yb
                        yi[2, k] = ys
                        s1 = fc1 / ys
                        hb1 = hi[1, k]
                    }
                    ht1 = htw
                }
            }
            if parabolic {
                lin = 0
                yi[2, k] = min(yi[2, k], hi[2, k] - hi[1, k] + 1)
                fsecv[k] = -1
            }
        }
        // D region tail.
        let hd: Float = 70
        let xtail: Float = 0.85
        let hex = hi[1, k] - xtail * yi[1, k]
        let fnx = 1 - xtail * xtail
        let alp = 2 * (hi[1, k] - hex) / (fnx * yi[1, k] * yi[1, k])
        let fsq = fnx * fexp(-alp * (hex - hd))
        htr[1, 1] = hd
        htr[5, 1] = hex
        var hdif = (htr[5, 1] - htr[1, 1]) * 0.25
        hdif = max(0, hdif)
        htr[4, 1] = htr[5, 1] - min(1.0, hdif)
        htr[2, 1] = htr[1, 1] + hdif
        htr[3, 1] = (htr[2, 1] + htr[4, 1]) * 0.5
        htr[11, 1] = hi[1, k]
        hdif = (htr[11, 1] - htr[5, 1]) / 6
        for ih in 6...10 { htr[ih, 1] = htr[ih - 1, 1] + hdif }
        htr[17, 1] = hi[1, k] + yi[1, k]
        hdif = (htr[17, 1] - htr[11, 1]) / 6
        for ih in 12...16 { htr[ih, 1] = htr[ih - 1, 1] + hdif }
        htr[11, 1] = 0.5 * (htr[10, 1] + htr[12, 1])
        htr[50, 1] = hi[3, k]
        var withF1 = false
        if fi[2, k] > 0 {
            hb2 = hi[3, k] - yi[3, k]
            hb1 = hi[2, k] - yi[2, k] + 0.00001
            withF1 = hb2 - hb1 > 0
        }
        if !withF1 {
            htr[18, 1] = hi[3, k] - yi[3, k]
            hdif = (htr[50, 1] - htr[18, 1]) / 32
            for ih in 19...49 { htr[ih, 1] = htr[ih - 1, 1] + hdif }
        } else {
            htr[18, 1] = hi[2, k] - yi[2, k]
            htr[18, 1] = max(htr[18, 1], htr[17, 1] + 1)
            htr[28, 1] = hi[2, k]
            hdif = (htr[28, 1] - htr[18, 1]) / 10
            for ih in 19...27 { htr[ih, 1] = htr[ih - 1, 1] + hdif }
            hdif = (htr[50, 1] - htr[28, 1]) / 22
            for ih in 29...49 { htr[ih, 1] = htr[ih - 1, 1] + hdif }
        }
        hb1 = max(hi[1, k], hi[2, k] - yi[2, k])
        for ih in 1...50 {
            var fnd: Float = 0, fne: Float = 0, fn1: Float = 0, fn2: Float = 0, fnval: Float = 0
            let h = htr[ih, 1]
            if h - hlow > 0 && h - hup < 0 { fnval = fup + asp * (h - hup) }
            if h - hex >= 0 {
                if h - hte <= 0 {
                    let z = (h - hi[1, k]) / yi[1, k]
                    fne = fce * (1 - z * z)
                }
            } else {
                fnd = fce * fsq * fexp(alp * (h - hd))
            }
            if fi[2, k] > 0 && hb1 - h <= 0 && h - ht1 <= 0 {
                if lin <= 0 {
                    let z = (h - hi[2, k]) / yi[2, k]
                    fn1 = fc1 * (1.0 - z * z)
                } else {
                    fn1 = s1 * (h - (hi[2, k] - yi[2, k]))
                }
            }
            if hb2 - h <= 0 {
                let z = (h - hi[3, k]) / yi[3, k]
                fn2 = fc2 * (1.0 - z * z)
            }
            fnsq[ih, 1] = max(fnd, fne, fnval, fn1, fn2)
        }
    }

    /// GENION: the vertical ionogram (frequency, true and virtual height) for control point `k`.
    func genion(_ k: Int) {
        let xtail: Float = 0.85
        let fex = fi[1, k] * fsqrt(1 - xtail * xtail)
        fvert[1, k] = 0.01
        fvert[4, k] = fex
        var fdif = (fvert[4, k] - fvert[1, k]) / 3
        fvert[2, k] = fvert[1, k] + fdif
        fvert[3, k] = fvert[2, k] + fdif
        fvert[9, k] = 0.957 * fi[1, k]
        fvert[10, k] = 0.99 * fi[1, k]
        fdif = (fvert[9, k] - fvert[4, k]) / 5
        fvert[5, k] = fvert[4, k] + fdif
        fvert[6, k] = fvert[5, k] + fdif
        fvert[7, k] = fvert[6, k] + fdif
        fvert[8, k] = fvert[7, k] + fdif
        fvert[11, k] = 1.05 * fi[1, k]
        fvert[30, k] = 0.99 * fi[3, k]
        fvert[29, k] = 0.98 * fi[3, k]
        fvert[28, k] = 0.96 * fi[3, k]
        fvert[27, k] = 0.92 * fi[3, k]
        if fi[2, k] <= 0 {
            fdif = (fvert[27, k] - fvert[11, k]) / 16
            for i in 12...26 { fvert[i, k] = fvert[i - 1, k] + fdif }
        } else {
            fvert[20, k] = 0.99 * fi[2, k]
            fdif = (fvert[20, k] - fvert[11, k]) / 9
            for i in 12...19 { fvert[i, k] = fvert[i - 1, k] + fdif }
            fvert[21, k] = 1.01 * fi[2, k]
            fdif = (fvert[27, k] - fvert[21, k]) / 6
            for i in 22...26 { fvert[i, k] = fvert[i - 1, k] + fdif }
        }
        // IEDP is always negative here: integrate the full profile.
        for i in 1...30 {
            let heights = gethp(fvert[i, k])
            hprim[i, k] = heights.hp
            htrue[i, k] = heights.ht
        }
    }

    /// F2DIS: standard deviation of the F2 MUF from the F2 M(3000) decile tables.
    func f2dis(fmuf: Float, ssn: Float, clat: Float, freq: Float, abb: Float) -> Float {
        let cl = clat * 57.29577
        var icc = fint(abb / 4.0 + 1.55)
        if icc > 6 { icc = 1 }
        let i = freq <= fmuf ? 8 : 0
        var icl = fint(9.5 - abs(cl) / 10)
        if icl < 1 { icl = 1 }
        if icl > 8 { icl = 8 }
        let iz = icl + i
        var j = 2
        if ssn <= 50 { j = 1 }
        if ssn > 100 { j = 3 }
        if cl <= 0 { j = j + 3 }
        var sig = abs(fmuf - f2d[iz, j, icc] * fmuf) / 1.28
        if sig <= 0.001 { sig = 0.001 }
        return sig
    }

    /// The hop geometry CURMUF repeats for each layer: the take-off angle
    /// and the sine of the incidence angle for the fewest hops that clear
    /// the minimum angle.
    private struct HopGeometry {
        var xhops: Float
        var cpsi: Float
        var spsi: Float
    }

    private func hopGeometry(hp: Float) -> HopGeometry {
        let dele = max(0.0, amind)
        let del = dele * d2r
        let phe = fasin(rz * fcos(del) / (rz + hp))
        let nhops = fint(0.5 * gcdkm / ((pio2 - del - phe) * rz))
        var psi = gcdkm / (2 * rz)
        let xhops = Float(nhops + 1)
        psi = psi / xhops
        return HopGeometry(xhops: xhops, cpsi: fcos(psi), spsi: fsin(psi))
    }

    /// CURMUF: MUF of each layer and of the circuit (sporadic E kept separate).
    func curmuf() {
        let fqdel: Float = 0.1
        let nfq = 4
        let ae: Float = 0.5, a1: Float = 0.5, a2: Float = 0.5
        let dzf: Float = 2000, bex: Float = 9.5
        let kt = 1
        let kr = kfx < 2 ? 1 : (kfx == 2 ? 2 : 3)
        // Select the controlling sample area.
        var ks: Int
        func byE() -> Int { fi[1, kt] - fi[1, kr] <= 0 ? kt : kr }
        if kr < 2 {
            ks = kt
        } else if kr == 2 {
            ks = byE()
        } else {
            if abs(fi[3, kt] - fi[3, kr]) - 0.01 <= 0 {
                ks = byE()
            } else if fi[3, kt] - fi[3, kr] < 0 {
                ks = kt
            } else {
                ks = kr
            }
        }
        lecden(ks)
        // Tangent frequencies.
        let xte = 1.0 / fsqrt(1.0 + ae * yi[1, ks] / hi[1, ks])
        let fxe = xte * fi[1, ks]
        var fx1: Float = 0
        if fi[2, ks] > 0 {
            let xt1 = 1.0 / fsqrt(1.0 + a1 * yi[2, ks] / hi[2, ks])
            fx1 = xt1 * fi[2, ks]
        }
        var xt2 = 1.0 / fsqrt(1.0 + a2 * yi[3, ks] / hi[3, ks])
        var fx2 = xt2 * fi[3, ks]
        if gcdkm - dzf < 0 {
            let a = -1 + 1 / xt2
            let beta = 1 + a * fexp(-bex * gcdkm / dzf)
            fx2 = beta * fx2
            xt2 = xt2 * beta
        }
        _ = xt2

        // E layer.
        let e = gethp(fxe)
        let hpe = e.hp, hte = e.ht
        var geometry = hopGeometry(hp: hpe)
        var psi = gcdkm / (2 * rz) / geometry.xhops
        var tanp = geometry.spsi / (1 - geometry.cpsi + hpe / rz)
        var phe = fatan(tanp)
        var del = pio2 - phe - psi
        var cdel = fcos(del)
        var sphe = rz * cdel / (rz + hte)
        var secp = 1 / fsqrt(1 - sphe * sphe)
        emuf[it] = fxe * secp
        let dele = del * r2d
        sigl[1] = max(0.01, 0.1 * emuf[it])
        sigu[1] = sigl[1]
        delmuf[1] = dele
        hpmuf[1] = hpe
        htmuf[1] = hte
        fvmuf[1] = fxe
        afmuf[1] = 0
        nhopmf[1] = Int(geometry.xhops)
        yfot[1] = emuf[it] - 1.28 * sigl[1]
        yhpf[1] = emuf[it] + 1.28 * sigu[1]
        ymuf[1] = emuf[it]

        // F2 layer.
        let f2 = gethp(fx2)
        let hp2 = f2.hp, ht2 = f2.ht
        geometry = hopGeometry(hp: hp2)
        psi = gcdkm / ((2 * rz) * geometry.xhops)
        tanp = geometry.spsi / (1 - geometry.cpsi + hp2 / rz)
        phe = fatan(tanp)
        del = pio2 - phe - psi
        cdel = fcos(del)
        sphe = rz * cdel / (rz + ht2)
        secp = 1 / fsqrt(1 - sphe * sphe)
        var fob1 = fx2 * secp
        var ntry = 0
        var fob2 = fob1
        var xhp = (hp2 - ht2) / rz
        var hpx2: Float = hp2
        while true {
            fob1 = fob2
            // Correction to Martyn's theorem.
            let xmut = sphe * sphe
            let xfsq = fob1 * fob1 / (fi[3, ks] * fi[3, ks])
            let sph = xfsq * xmut * xhp * (ht2 + 2 * (rz + ht2) * xhp)
            hpx2 = hp2 + sph
            tanp = geometry.spsi / (1 - geometry.cpsi + hpx2 / rz)
            phe = fatan(tanp)
            del = pio2 - phe - psi
            cdel = fcos(del)
            sphe = rz * cdel / (rz + ht2)
            secp = 1 / fsqrt(1 - sphe * sphe)
            fob2 = fx2 * secp
            ntry = ntry + 1
            if abs(fob2 - fob1) - fqdel <= 0 { break }
            if nfq - ntry < 0 { break }
        }
        f2muf[it] = fob2
        let del2 = del * r2d
        freq = 0.9 * f2muf[it]
        var sig = f2dis(fmuf: f2muf[it], ssn: ssn, clat: clat[ks], freq: freq, abb: clck[ks])
        sigl[3] = max(0.01, sig)
        freq = 1.1 * f2muf[it]
        sig = f2dis(fmuf: f2muf[it], ssn: ssn, clat: clat[ks], freq: freq, abb: clck[ks])
        sigu[3] = max(0.01, sig)
        delmuf[3] = del2
        hpmuf[3] = hpx2
        htmuf[3] = ht2
        fvmuf[3] = fx2
        afmuf[3] = 0
        nhopmf[3] = Int(geometry.xhops)
        yfot[3] = f2muf[it] - 1.28 * sigl[3]
        ymuf[3] = f2muf[it]
        yhpf[3] = f2muf[it] + 1.28 * sigu[3]

        if fi[2, ks] <= 0 {
            // No F1 layer: copy the E layer.
            f1muf[it] = emuf[it]
            sigl[2] = sigl[1]
            sigu[2] = sigu[1]
            delmuf[2] = delmuf[1]
            hpmuf[2] = hpmuf[1]
            htmuf[2] = htmuf[1]
            fvmuf[2] = fvmuf[1]
            afmuf[2] = afmuf[1]
            nhopmf[2] = nhopmf[1]
            yfot[2] = yfot[1]
            ymuf[2] = ymuf[1]
            yhpf[2] = yhpf[1]
        } else {
            let f1 = gethp(fx1)
            let hp1 = f1.hp, ht1 = f1.ht
            geometry = hopGeometry(hp: hp1)
            psi = gcdkm / (2 * rz) / geometry.xhops
            tanp = geometry.spsi / (1 - geometry.cpsi + hp1 / rz)
            phe = fatan(tanp)
            del = pio2 - phe - psi
            cdel = fcos(del)
            sphe = rz * cdel / (rz + ht1)
            secp = 1 / fsqrt(1 - sphe * sphe)
            fob1 = fx1 * secp
            ntry = 0
            fob2 = fob1
            xhp = (hp1 - ht1) / rz
            var hpy2 = hp1
            while true {
                fob1 = fob2
                let xfsq = fob1 * fob1 / (fi[2, ks] * fi[2, ks])
                let xmut = sphe * sphe
                let sph = xfsq * xmut * xhp * (ht1 + 2 * (rz + ht1) * xhp)
                hpy2 = hp1 + sph
                tanp = geometry.spsi / (1 - geometry.cpsi + hpy2 / rz)
                phe = fatan(tanp)
                del = pio2 - phe - psi
                cdel = fcos(del)
                sphe = rz * cdel / (rz + ht1)
                secp = 1 / fsqrt(1 - sphe * sphe)
                fob2 = fx1 * secp
                ntry = ntry + 1
                if abs(fob2 - fob1) - fqdel <= 0 { break }
                if nfq - ntry < 0 { break }
            }
            f1muf[it] = fob2
            let del1 = del * r2d
            sigl[2] = max(0.01, 0.1 * f1muf[it])
            sigu[2] = sigl[2]
            delmuf[2] = del1
            hpmuf[2] = hpy2
            htmuf[2] = ht1
            fvmuf[2] = fx1
            afmuf[2] = 0
            nhopmf[2] = Int(geometry.xhops)
            yfot[2] = f1muf[it] - 1.28 * sigl[2]
            ymuf[2] = f1muf[it]
            yhpf[2] = f1muf[it] + 1.28 * sigu[2]
        }

        // Sporadic E MUF, for a reflection probability of one half.
        let dels = max(0.0, amind)
        esmuf[it] = 1000
        var ksx = 1
        var tdelx: Float = 0, secpx: Float = 0, hop: Float = 0
        for k in 1...km {
            if fs[2, k] <= 0 { continue }
            del = dels * d2r
            phe = fasin(rz * fcos(del) / (rz + hs[k]))
            let nhops = fint(0.5 * gcdkm / ((pio2 - del - phe) * rz))
            hop = Float(nhops + 1)
            psi = 0.5 * gcd / hop
            let tdel = (fcos(psi) - rz / (rz + hs[k])) / fsin(psi)
            cdel = 1 / fsqrt(1 + tdel * tdel)
            sphe = rz * cdel / (rz + hs[k])
            secp = 1 / fsqrt(1 - sphe * sphe)
            let esd = fs[2, k] * secp
            if esmuf[it] - esd > 0 {
                esmuf[it] = esd
                ksx = k
                tdelx = tdel
                secpx = secp
            }
        }
        if esmuf[it] - 1000 >= 0 {
            delmuf[4] = 0; hpmuf[4] = 0; htmuf[4] = 0; fvmuf[4] = 0; afmuf[4] = 0
            nhopmf[4] = 0
            yfot[4] = 0; ymuf[4] = 0; yhpf[4] = 0; sigl[4] = 0; sigu[4] = 0
            esmuf[it] = 0
        } else {
            delmuf[4] = fatan(tdelx) * r2d
            hpmuf[4] = hs[ksx]
            htmuf[4] = hpmuf[4]
            fvmuf[4] = fs[2, ksx]
            afmuf[4] = 0
            nhopmf[4] = Int(hop)
            ymuf[4] = esmuf[it]
            let xes = fi[1, ks] - fs[1, ksx]
            let zdif: Float = 0.1
            let fzero: Float = 0.1
            sigu[4] = (fs[3, ksx] - fs[2, ksx]) / 1.28
            let sigb = (fs[2, ksx] - fs[1, ksx]) / 1.28
            let sigz = (fs[2, ksx] - fzero) / 3.1
            if xes + zdif < 0 {
                sigl[4] = sigb
            } else if xes - zdif < 0 {
                sigl[4] = sigb + ((xes + zdif) / (2 * zdif)) * (sigu[4] - sigb)
            } else {
                sigl[4] = sigu[4]
            }
            sigl[4] = min(sigl[4], sigz)
            fs[1, ksx] = fs[2, ksx] - 1.28 * sigl[4]
            sigl[4] = sigl[4] * secpx
            sigu[4] = sigu[4] * secpx
            yfot[4] = esmuf[it] - 1.28 * sigl[4]
            yhpf[4] = esmuf[it] + 1.28 * sigu[4]
        }

        // Circuit MUF (sporadic E not included).
        fot[it] = max(yfot[1], yfot[2], yfot[3])
        allmuf[it] = max(emuf[it], f1muf[it], f2muf[it])
        hpf[it] = max(yhpf[1], yhpf[2], yhpf[3])
        if emuf[it] - allmuf[it] >= 0 {
            angmuf[it] = delmuf[1]
            modmuf = 1
        } else if f1muf[it] - allmuf[it] >= 0 {
            angmuf[it] = delmuf[2]
            modmuf = 2
        } else {
            angmuf[it] = delmuf[3]
            modmuf = 3
        }
    }
}
