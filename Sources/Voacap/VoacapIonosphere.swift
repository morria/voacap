import Foundation

// Coefficient loading, path geometry and the ionospheric parameters at the
// path's sample areas: REDMAP, GEOM, MAGVAR, MAGFIN, NOISY, GEOTIM, CNGTIM,
// SANG, VIRTIM, TIMVAR, EF1VAR, VERSY, F2VAR, ESIND and IONSET.

extension VoacapEngine {
    /// REDMAP: load one month's coefficients and interpolate them to the sunspot number.
    func redmap(ssn: Float, month: Int) {
        self.ssn = ssn
        self.month = month
        let m = coefficients.months[month - 1]
        for i in 0..<60 { ikim.d[i] = Int(m.ikim[i]) }
        dud.d = m.dud
        fam.d = m.fam
        sys.d = m.sys
        f2d.d = m.f2d
        perr.d = m.perr
        // P(29,16,8): six noise maps, then the land-mass map, then the height ratio map.
        for i in 0..<2784 { p.d[i] = m.fakp[i] }
        for i in 0..<464 { p.d[2784 + i] = m.fakmap[i] }
        for i in 0..<12 { abp.d[i] = m.fakabp[i] }
        for i in 0..<6 { abp.d[12 + i] = m.abmap[i] }
        anew.d = Array(m.fone[0..<3])
        bnew.d = Array(m.fone[3..<6])
        achi.d = Array(m.fone[6..<8])
        bchi.d = Array(m.fone[8..<10])

        let ssn100 = 100 - ssn
        for i in 0..<988 { cofion.d[977 + i] = (m.xf2cof[i] * ssn100 + m.xf2cof[988 + i] * ssn) / 100 }
        let ssn150 = 150 - ssn
        let ssn10 = ssn - 10
        for i in 0..<427 { cofion.d[275 + i] = (m.xesmcf[i] * ssn150 + m.xesmcf[427 + i] * ssn10) / 140 }
        let ssn125 = 125 - ssn
        let ssn25 = ssn - 25
        for i in 0..<464 { p.d[3248 + i] = (m.xpmap[i] * ssn125 + m.xpmap[464 + i] * ssn25) / 100 }
        for i in 0..<275 {
            cofion.d[i] = (m.xeslcf[i] * ssn150 + m.xeslcf[275 + i] * ssn10) / 140
            cofion.d[702 + i] = (m.xesucf[i] * ssn150 + m.xesucf[275 + i] * ssn10) / 140
        }
        for i in 0..<441 { cofion.d[1965 + i] = (m.xfm3cf[i] * ssn100 + m.xfm3cf[441 + i] * ssn) / 100 }
        for i in 0..<198 { cofion.d[2406 + i] = (m.xercof[i] * ssn150 + m.xercof[198 + i] * ssn10) / 140 }
    }

    /// GEOM: great-circle distance, bearings and the sample areas along the path.
    func geom() {
        let epslon: Float = 1.0e-7
        let ddx: Float = 0.03
        let glt = gMagneticLat * d2r
        let glg = gMagneticLon * d2r
        if abs(tlatd - rlatd) <= ddx && abs(tlongd - rlongd) <= ddx {
            // Move the receiver away from the transmitter a little.
            var dd = ddx
            if rlatd < 0 { dd = -ddx }
            rlatd = tlatd - dd
        }
        tlat = tlatd * d2r
        tlong = tlongd * d2r
        rlat = rlatd * d2r
        rlong = rlongd * d2r
        if abs(rlatd) > 89.9 { rlong = 0 }
        if abs(tlatd) > 89.9 { tlong = 0 }
        dlong = tlong - rlong
        if abs(dlong) > pi { dlong = dlong - fsign(pi2, dlong) }
        if npsl == 1 { dlong = dlong - fsign(pi2, dlong) }
        var qcos = fsin(tlat) * fsin(rlat) + fcos(tlat) * fcos(rlat) * fcos(dlong)
        if abs(qcos) > 1 { qcos = fsign(1, qcos) }
        gcd = facos(qcos)
        if gcd < 0.000001 { gcd = 0.000001 }
        if npsl == 1 { gcd = pi2 - gcd }
        gcdkm = gcd * rz
        if fcos(tlat) - epslon <= 0 {
            btr = tlat <= 0 ? 0 : pi
        } else {
            qcos = (fsin(rlat) - fsin(tlat) * fcos(gcd)) / (fcos(tlat) * fsin(gcd))
            if abs(qcos) > 1 { qcos = fsign(1, qcos) }
            btr = facos(qcos)
        }
        if dlong > 0 { btr = pi2 - btr }
        var brt: Float
        if fcos(rlat) - epslon <= 0 {
            brt = rlat <= 0 ? 0 : pi
        } else {
            qcos = (fsin(tlat) - fsin(rlat) * fcos(gcd)) / (fcos(rlat) * fsin(gcd))
            if abs(qcos) > 1 { qcos = fsign(1, qcos) }
            brt = facos(qcos)
        }
        if dlong < 0 { brt = pi2 - brt }
        btrd = btr * r2d
        brtd = brt * r2d
        // Sample areas.
        if gcdkm <= 2000.01 {
            km = 1
            rd[1] = gcd / 2
        } else if gcdkm <= 4000 {
            km = 3
            rd[1] = 1000 / rz
            rd[2] = gcd / 2
            rd[3] = gcd - rd[1]
        } else {
            km = 5
            rd[1] = 1000 / rz
            rd[2] = rd[1] + rd[1]
            rd[3] = gcd / 2
            rd[4] = gcd - rd[2]
            rd[5] = gcd - rd[1]
        }
        for i in 1...km {
            let drf = rd[i]
            let ctlat = fcos(tlat)
            var rflt: Float
            var rflg: Float
            if ctlat - epslon < 0 {
                rflt = tlat - fsign(drf, tlat)
                if abs(rflt) > pio2 { rflt = pio2 * fsign(1, rflt) }
                rflg = rlong
            } else {
                qcos = fcos(drf) * fsin(tlat) + fsin(drf) * fcos(tlat) * fcos(btr)
                if abs(qcos) > 1 { qcos = fsign(1, qcos) }
                rflt = pio2 - facos(qcos)
                if fcos(rflt) - epslon <= 0 {
                    rflg = tlong
                } else {
                    qcos = (fcos(drf) - fsin(rflt) * fsin(tlat)) / (fcos(rflt) * fcos(tlat))
                    if abs(qcos) > 1 { qcos = fsign(1, qcos) }
                    rflg = facos(qcos)
                    if drf >= pi { rflg = pi2 - rflg }
                    rflg = tlong - fsign(rflg, dlong)
                    if abs(rflg) > pi { rflg = rflg - fsign(pi2, rflg) }
                }
            }
            qcos = fsin(glt) * fsin(rflt) + fcos(glt) * fcos(rflt) * fcos(rflg - glg)
            if abs(qcos) > 1 { qcos = fsign(1, qcos) }
            let gat = facos(qcos)
            glat[i] = pio2 - gat
            clat[i] = rflt
            clong[i] = rflg
            let clg = magvar(cenlat: rflt, cenlg: rflg, ii: i)
            let rfltd = rflt * r2d
            var clgd = clg * r2d
            if clgd < 0 { clgd = clgd + 360 }
            let wld = noisy(7, rfltd, clgd)
            if wld >= 0 {
                sigpat[i] = 0.001
                epspat[i] = 4
            } else {
                sigpat[i] = 5
                epspat[i] = 80
            }
        }
        alatd = abs(clat[1]) * r2d
        if km != 1 { alatd = abs((clat[1] + clat[2] + clat[3]) / 3.0) * r2d }
    }

    /// MAGVAR: gyrofrequency and modified magnetic dip at a sample area.
    /// Returns the east longitude.
    func magvar(cenlat: Float, cenlg: Float, ii: Int) -> Float {
        let c360 = 360 * d2r
        let clg = cenlg < 0 ? c360 + cenlg : cenlg
        let une = magfin(lat: cenlat, lon: clg, height: 300000)
        let hhm = fsqrt(une.0 * une.0 + une.1 * une.1 + une.2 * une.2)
        gyz[ii] = 2.8 * hhm
        let tmp = une.1 * une.1 + une.2 * une.2
        var gob = fcos(cenlat)
        if gob <= 0.000001 { gob = 0.000001 }
        gmdip[ii] = fatan(fatan(-une.0 / fsqrt(tmp)) / fsqrt(gob))
        return clg
    }

    private static let magCT: [Float] = [
        0.0, 0.0, 0.33333333, 0.26666667, 0.25714286, 0.25396825, 0.25252525,
        0, 0, 0, 0.2, 0.22857142, 0.23809523, 0.24242424,
        0, 0, 0, 0, 0.14285714, 0.19047619, 0.21212121,
        0, 0, 0, 0, 0, 0.11111111, 0.16161616,
        0, 0, 0, 0, 0, 0, 0.09090909,
        0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0]
    private static let magG: [Float] = [
        0.0, 0.304112, 0.024035, -0.031518, -0.041794, 0.016256, -0.019523,
        0, 0.021474, -0.051253, 0.062130, -0.045298, -0.034407, -0.004853,
        0, 0, -0.013381, -0.024898, -0.021795, -0.019447, 0.003212,
        0, 0, 0, -0.006496, 0.007008, -0.000608, 0.021413,
        0, 0, 0, 0, -0.002044, 0.002775, 0.001051,
        0, 0, 0, 0, 0, 0.000697, 0.000227,
        0, 0, 0, 0, 0, 0, 0.001115]
    private static let magH: [Float] = [
        0, 0, 0, 0, 0, 0, 0,
        0, -0.057989, 0.033124, 0.014870, -0.011825, -0.000796, -0.005758,
        0, 0, -0.001579, -0.004075, 0.010006, -0.002, -0.008735,
        0, 0, 0, 0.00021, 0.00043, 0.004597, -0.003406,
        0, 0, 0, 0, 0.001385, 0.002421, -0.000118,
        0, 0, 0, 0, 0, -0.001218, -0.001116,
        0, 0, 0, 0, 0, 0, -0.000325]

    /// MAGFIN: the magnetic field vector (up, north, east) from the 1960 Jensen and Cain model.
    func magfin(lat: Float, lon: Float, height: Float) -> (Float, Float, Float) {
        var pp = A2(7, 7), dp = A2(7, 7)
        pp[1, 1] = 1
        var cp = A1(7), sp = A1(7), aor = A1(7)
        cp[1] = 1
        @inline(__always) func ct(_ n: Int, _ m: Int) -> Float { Self.magCT[(m - 1) * 7 + n - 1] }
        @inline(__always) func g(_ n: Int, _ m: Int) -> Float { Self.magG[(m - 1) * 7 + n - 1] }
        @inline(__always) func h(_ n: Int, _ m: Int) -> Float { Self.magH[(m - 1) * 7 + n - 1] }
        let rdLimit: Float = 1.56905124
        let hc: Float = 6371200
        var p2 = lon
        var p1 = lat
        if p1 > rdLimit {
            p1 = rdLimit
            p2 = 0
        } else if p1 < -rdLimit {
            p1 = -rdLimit
            p2 = 0
        }
        let phi = p2
        let ar = hc / (hc + height)
        let c = fsin(p1)
        let s = fsqrt(cp[1] - c * c)
        sp[2] = fsin(phi)
        cp[2] = fcos(phi)
        aor[1] = ar * ar
        aor[2] = aor[1] * ar
        for m in 3...7 {
            sp[m] = sp[2] * cp[m - 1] + cp[2] * sp[m - 1]
            cp[m] = cp[2] * cp[m - 1] - sp[2] * sp[m - 1]
            aor[m] = ar * aor[m - 1]
        }
        var bv: Float = 0, bn: Float = 0, bphi: Float = 0
        for n in 2...7 {
            let fn = Float(n)
            var sumr: Float = 0, sumt: Float = 0, sump: Float = 0
            for m in 1...n {
                if n == m {
                    pp[n, n] = s * pp[n - 1, n - 1]
                    dp[n, n] = s * dp[n - 1, n - 1] + c * pp[n - 1, n - 1]
                } else if n == 2 {
                    pp[2, 1] = c
                    dp[2, 1] = -s
                } else {
                    pp[n, m] = c * pp[n - 1, m] - ct(n, m) * pp[n - 2, m]
                    dp[n, m] = c * dp[n - 1, m] - s * pp[n - 1, m] - ct(n, m) * dp[n - 2, m]
                }
                let fm = Float(m - 1)
                let ts = g(n, m) * cp[m] + h(n, m) * sp[m]
                sumr = sumr + pp[n, m] * ts
                sumt = sumt + dp[n, m] * ts
                sump = sump + fm * pp[n, m] * (-g(n, m) * sp[m] + h(n, m) * cp[m])
            }
            bv = bv + aor[n] * fn * sumr
            bn = bn - aor[n] * sumt
            bphi = bphi - aor[n] * sump
        }
        return (-bv, bn, -bphi / s)
    }

    /// NOISY: evaluate one of the two-variable Fourier maps (1–6 atmospheric
    /// noise, 7 land mass, 8 F2 height-to-semithickness ratio).
    func noisy(_ kj: Int, _ xla: Float, _ ceg: Float) -> Float {
        let alf = abp[1, kj]
        let bet = abp[2, kj]
        let lm = kj == 8 ? 15 : 29
        let ln = kj == 8 ? 10 : 15
        var sx = A1(15), zz = A1(29)
        var q = 0.0087266466 * ceg
        var c1 = fcos(q)
        var s1 = fsin(q)
        sx[1] = s1
        var cx = c1
        for k in 2...ln {
            let tx = sx[k - 1]
            sx[k] = tx * c1 + cx * s1
            cx = cx * c1 - tx * s1
        }
        for j in 1...lm {
            var r: Float = 0
            for k in 1...ln { r = r + sx[k] * p[j, k, kj] }
            zz[j] = r + p[j, 16, kj]
        }
        q = 0.01745329252 * (xla + 90)
        s1 = fsin(q)
        c1 = fcos(q)
        sx[1] = s1
        cx = c1
        var r: Float = 0
        for k in 1...lm {
            r = r + sx[1] * zz[k]
            let ss = sx[1] * c1 + cx * s1
            cx = cx * c1 - sx[1] * s1
            sx[1] = ss
        }
        return r + alf + bet * q
    }

    /// CNGTIM: convert between universal and local mean time at a longitude.
    func cngtim(_ time: inout Float, _ xtheta: Float, _ isw: Int) {
        if xtheta == 0 { return }
        let theta = xtheta < 0 ? 360 + xtheta : xtheta
        let fint1 = Float(fint(theta / 180))
        let fint2 = Float(fint(theta / 360))
        time = time + Float(isw) * (theta / 15.0 - fint1 * 24 + fint2 * 24)
        if time < 0 { time = time + 24 }
        time = time - fint2 * 24
        while true {
            if time < 0 {
                time = time + 24
            } else if time - 24 > 0 {
                time = time - 24
            } else {
                break
            }
        }
    }

    /// GEOTIM: local times at the terminals and sample areas for a UT hour.
    func geotim(_ jt: Int) {
        var ckc = Float(jt)
        gmt = ckc
        cngtim(&ckc, tlong * r2d, 1)
        let temp = ckc
        it = fint(gmt)
        if it <= 0 { it = 24 }
        utime[it] = gmt
        xlmt[it] = temp
        ckc = gmt
        cngtim(&ckc, rlong * r2d, 1)
        gmtr = ckc
        for ii in 1...km {
            ckc = gmt
            cngtim(&ckc, clong[ii] * r2d, 1)
            clck[ii] = ckc
        }
    }

    /// SANG: limit the elevation-angle scan by path length.
    func sang() {
        let nangx = [40, 34, 29, 24, 19, 14, 12, 9]
        var id = fint(gcdkm / 2000 + 1)
        id = min(id, 8)
        var iang = nangx[id - 1]
        if amind > 0 && ang[iang] - 10 - amind < 0 {
            for iaIndex in 1...40 {
                iang = iaIndex
                if ang[iang] - 8.5 - amind >= 0 { break }
            }
        }
        nang = iang
    }

    /// VIRTIM: collapse the coefficient maps to the current universal time.
    func virtim() {
        var c = A1(8), s = A1(8)
        let time = (15 * gmt - 180) * d2r
        c[1] = fcos(time)
        s[1] = fsin(time)
        for jb in 2...8 {
            c[jb] = c[1] * c[jb - 1] - s[1] * s[jb - 1]
            s[jb] = c[1] * s[jb - 1] + s[1] * c[jb - 1]
        }
        for iz in 1...6 {
            let j = ikim[10, iz]
            let i = ikim[9, iz] + 1
            for jb in 1...i {
                let isubb = ia[iz] + (jb - 1) * ib[iz]
                let isuba = icab[iz] + jb - 1
                ab[isuba] = cofion[isubb]
                if j >= 1 {
                    for ka in 1...j {
                        let isubc = isubb + 2 * ka - 1
                        ab[isuba] = ab[isuba] + s[ka] * cofion[isubc] + c[ka] * cofion[isubc + 1]
                    }
                }
            }
        }
    }

    /// TIMVAR: solar zenith angle, E and F1 layers, and the F2 height ratio.
    func timvar() {
        let ssl = (180 - 15 * gmt) * d2r
        for ii in 1...km {
            let cenlg = clong[ii]
            let cenlat = clat[ii]
            let cendog = abs(cenlat - d2r * sunLat(1, month))
            var ssp = sunLat(2, month) * d2r
            let cencat = abs(cenlat - ssp)
            if cendog - cencat > 0 { ssp = sunLat(1, month) * d2r }
            var clock = gmt + cenlg / 0.261799387
            if 24 - clock < 0 {
                clock = clock - 24
            } else if clock <= 0 {
                clock = clock + 24
            }
            clck[ii] = clock
            var z = cenlg - ssl
            var cycen = fsin(cenlat) * fsin(ssp) + fcos(cenlat) * fcos(ssp) * fcos(z)
            cycen = facos(cycen)
            ef1var(cycen: abs(cycen), ii: ii)
            let gm = abs(glat[ii]) * r2d - 45
            z = cycen * r2d
            if clock - 12.0 < 0 { z = -z }
            z = z + 180.0
            var ym = noisy(8, gm, z)
            if ym - 2.0 < 0 { ym = 2.0 }
            rat[ii] = ym
            abiy[ii] = -0.04 + fexp(-2.937 + 0.8445 * fi[1, ii])
        }
    }

    /// EF1VAR: E and F1 layer critical frequency, height and semithickness.
    func ef1var(cycen: Float, ii: Int) {
        let betae: Float = 5.5, betaf1: Float = 4.0
        versy(6, 6, ii)
        if gamma[6] - 0.36 < 0 { gamma[6] = 0.36 * fsqrt(1 + 0.0098 * ssn) }
        fi[1, ii] = gamma[6] * psc[1]
        hi[1, ii] = 110
        yi[1, ii] = hi[1, ii] / betae
        let z = cycen * r2d
        let cosz = fcos(cycen)
        let cosdi = fcos(gmdip[ii])
        let zmax = achi[1] + bchi[1] * ssn + (achi[2] + bchi[2] * ssn) * cosdi
        zenmax[ii] = zmax
        zenang[ii] = z
        if z <= zmax {
            fi[2, ii] = (anew[1] + bnew[1] * ssn) + (anew[2] + bnew[2] * ssn) * cosz
                + (anew[3] + bnew[3] * ssn) * cosz * cosz
            fi[2, ii] = fi[2, ii] * psc[2]
            hi[2, ii] = 165 + 0.6428 * z
            yi[2, ii] = hi[2, ii] / betaf1
        } else {
            fi[2, ii] = 0
            hi[2, ii] = 0
            yi[2, ii] = 0
        }
    }

    /// VERSY: evaluate maps `kik` through `lik` at a sample area
    /// (1–3 sporadic E, 4 foF2, 5 M(3000)F2, 6 foE).
    func versy(_ kik: Int, _ lik: Int, _ ii: Int) {
        var g = [Double](repeating: 0, count: 80)
        let c360 = 360 * d2r
        let cenlg = clong[ii]
        let cenlat = clat[ii]
        let clg = cenlg < 0 ? c360 + cenlg : cenlg
        let gob = fcos(cenlat)
        for iz in kik...lik {
            let x = Double(iz != 6 ? gmdip[ii] : cenlat)
            let y = clg
            let i = ikim[9, iz] + 1
            let k = ikim[1, iz]
            let sx = sin(x)
            g[1] = 1
            g[2] = sx
            if k >= 2 {
                for ka in 2...k { g[ka + 1] = sx * g[ka] }
            }
            var kdif = ikim[2, iz] - k
            if kdif != 0 {
                var jg = 1
                var cx = Double(gob)
                var t = Double(y)
                while true {
                    let kk = ikim[jg, iz] + 4
                    g[kk - 2] = cx * cos(t)
                    g[kk - 1] = cx * sin(t)
                    let lo = ikim[jg + 1, iz]
                    if kdif != 2 && lo >= kk {
                        var ka = kk
                        while ka <= lo {
                            g[ka] = sx * g[ka - 2]
                            g[ka + 1] = sx * g[ka - 1]
                            ka += 2
                        }
                    }
                    if jg == 8 { break }
                    kdif = ikim[jg + 2, iz] - lo
                    if kdif == 0 { break }
                    cx = cx * Double(gob)
                    jg = jg + 1
                    t = Double(Float(jg) * y)
                }
            }
            var value = Float(g[1] * Double(ab[icab[iz]]))
            if i >= 2 {
                for jb in 2...i {
                    value = Float(Double(value) + Double(ab[icab[iz] + jb - 1]) * g[jb])
                }
            }
            gamma[iz] = value
        }
    }

    /// F2VAR: F2 layer critical frequency, height and semithickness.
    func f2var() {
        let xf1: Float = 1.1, delz: Float = 2
        for ii in 1...km {
            versy(4, 5, ii)
            f2m3[ii] = gamma[5]
            hpf2[ii] = 1490 / gamma[5] - 176
            fi[3, ii] = (gamma[4] + 0.5 * gyz[ii]) * psc[3]
            let zmax = zenmax[ii]
            let z = zenang[ii]
            fi[2, ii] = min(fi[2, ii], fi[3, ii] - 0.2)
            let fc = 0.834 * fi[3, ii]
            let ec = fi[1, ii]
            var fcec = fc / ec
            fcec = max(fcec, 1.1)
            var ret = fcec * alog((fcec + 1) / (fcec - 1))
            ret = (ret - 2) * yi[1, ii]
            let fc1 = fi[2, ii]
            if fc1 > 0 {
                var ffec = fc / fc1
                var rft: Float
                if zmax - delz - z >= 0 {
                    ffec = max(ffec, xf1)
                    rft = ffec * alog((ffec + 1) / (ffec - 1))
                    rft = 0.5 * yi[2, ii] * (rft - 2)
                } else {
                    // Near the day–night line: lift F1 into F2 and fade its retardation.
                    ffec = max(ffec, xf1)
                    let zn = zmax - delz
                    let hn = 165 + 0.6428 * zn
                    let yn = hn * (yi[2, ii] / hi[2, ii])
                    var rfn = ffec * alog((ffec + 1) / (ffec - 1))
                    rfn = 0.5 * yn * (rfn - 2)
                    let sz = (z - zn) / delz
                    rft = rfn * (1 - sz)
                    let hm = hpf2[ii] - ret
                    let ym = hm / rat[ii]
                    let dhn = (hm - ym) - (hn - yn)
                    if dhn > 0 {
                        let dh = dhn * (1 - sz)
                        hi[2, ii] = (hm - ym) - dh + yi[2, ii]
                        if fc1 - fc > 0 {
                            let y1max = yn + (ym - yn) * (fc1 / fi[3, ii] - 0.834) / 0.166
                            yi[2, ii] = yn + (y1max - yn) * sz
                            hi[2, ii] = (hm - ym) - dh + yi[2, ii]
                        }
                    }
                }
                ret = ret + rft
            }
            hi[3, ii] = hpf2[ii] - ret
            yi[3, ii] = hi[3, ii] / rat[ii]
        }
    }

    /// ESIND: sporadic E critical frequencies (lower decile, median, upper decile).
    func esind() {
        for ii in 1...km {
            hs[ii] = 110
            versy(1, 3, ii)
            fs[1, ii] = gamma[3] * psc[4]
            fs[2, ii] = gamma[2] * psc[4]
            fs[3, ii] = gamma[1] * psc[4]
        }
    }

    /// IONSET: move each layer's parameters into the control-point slots and
    /// make the layers consistent with one another.
    func ionset() {
        if km <= 1 {
            kfx = 1
        } else if km <= 3 {
            kfx = 2
            fi[3, 1] = fi[3, 2]; yi[3, 1] = yi[3, 2]; hi[3, 1] = hi[3, 2]
            fi[1, 2] = fi[1, 3]; yi[1, 2] = yi[1, 3]; hi[1, 2] = hi[1, 3]
            fi[2, 2] = fi[2, 3]; yi[2, 2] = yi[2, 3]; hi[2, 2] = hi[2, 3]
        } else {
            kfx = 3
            fi[3, 1] = fi[3, 2]; yi[3, 1] = yi[3, 2]; hi[3, 1] = hi[3, 2]
            for i in 1...3 {
                fi[i, 2] = fi[i, 3]; yi[i, 2] = yi[i, 3]; hi[i, 2] = hi[i, 3]
            }
            fi[1, 3] = fi[1, 5]; yi[1, 3] = yi[1, 5]; hi[1, 3] = hi[1, 5]
            fi[2, 3] = fi[2, 5]; yi[2, 3] = yi[2, 5]; hi[2, 3] = hi[2, 5]
            fi[3, 3] = fi[3, 4]; yi[3, 3] = yi[3, 4]; hi[3, 3] = hi[3, 4]
        }
        for k in 1...kfx {
            if fi[2, k] > 0 {
                if fi[2, k] - fi[1, k] - 0.2 <= 0 || fi[3, k] - fi[2, k] - 0.2 <= 0 {
                    fi[2, k] = 0
                    yi[2, k] = 0
                    hi[2, k] = 0
                } else {
                    hi[2, k] = min(hi[2, k], hi[3, k])
                }
            }
            let hdif = hi[3, k] - hi[1, k] - 2
            yi[3, k] = min(yi[3, k], hdif)
        }
    }
}
