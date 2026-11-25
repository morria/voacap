import Foundation

// Shared pieces of the system-performance model: SELMOD, ANOIS1, FOBBY,
// ALOSFV, SETLNG, SYSSY, SIGDIS, XLIN, PRBMUF, FNORML, GENFAM, GAIN, PWRDB
// and ALLMODES.

/// An antenna as VOACAP sees it: gain in dBi for each whole MHz from 1 to 30
/// and each whole degree of elevation from 0 to 90, along the path bearing.
public struct VoacapAntenna: Sendable {
    /// gain[(elevation) * 30 + (MHz - 1)]
    public var gain: [Float]
    /// Efficiency in dB for each whole MHz from 1 to 30.
    public var efficiency: [Float]
    public var minFrequency: Float = 2
    public var maxFrequency: Float = 30

    public static func isotropic(gain: Float) -> VoacapAntenna {
        VoacapAntenna(gain: Array(repeating: gain, count: 30 * 91), efficiency: Array(repeating: 0, count: 30))
    }

    /// Build the table from a pattern function of (frequency in MHz, elevation in degrees).
    public static func pattern(_ dBi: (Float, Float) -> Float) -> VoacapAntenna {
        var gain = [Float](repeating: 0, count: 30 * 91)
        for elevation in 0...90 {
            for mhz in 1...30 {
                gain[elevation * 30 + mhz - 1] = dBi(Float(mhz), Float(elevation))
            }
        }
        return VoacapAntenna(gain: gain, efficiency: Array(repeating: 0, count: 30))
    }

    @inline(__always) func array(_ i: Int, _ j: Int) -> Float { gain[(j - 1) * 30 + i - 1] }
}

extension VoacapEngine {
    /// XLIN: linear interpolation in column `k` of two tables.
    func xlin(_ x: Float, _ xn: A2, _ yn: A2, _ ntab: Int, _ k: Int) -> Float {
        if xn[1, k] - x > 0 { return yn[1, k] }
        for j in 1...(ntab - 1) {
            let a = xn[j, k], bb = xn[j + 1, k]
            if a == bb {
                if x == a { return yn[j, k] }
                continue
            } else if a < bb {
                if a - x <= 0 && !(x - bb < 0) { continue }
            } else {
                if a - x < 0 { continue }
                if x - bb <= 0 { continue }
            }
            return yn[j, k] + (x - a) * (yn[j + 1, k] - yn[j, k]) / (bb - a)
        }
        return yn[ntab, k]
    }

    /// The polynomial approximation to the normal distribution's upper tail
    /// that PRBMUF and FNORML share.
    @inline(__always)
    private func normalTail(_ ypx: Float) -> Float {
        var yp = ypx
        if ypx < 0 { yp = -ypx }
        yp = min(5.0, yp)
        var qx = 1 + yp * (0.196854 + yp * (0.115194 + yp * (0.000344 + yp * 0.019527)))
        qx = qx * qx * qx * qx
        qx = 1 / qx
        qx = 0.5 * qx
        return qx
    }

    /// FNORML: cumulative normal distribution.
    func fnorml(_ ypx: Float) -> Float {
        let qx = normalTail(ypx)
        return ypx < 0 ? qx : 1 - qx
    }

    /// PRBMUF: probability that the layer's MUF exceeds the operating frequency.
    func prbmuf(_ fmhz: Float, _ fgo: Float, _ fset: Float, _ il: Int) -> Float {
        var z = fmhz - fgo
        if !(fset > 0) { return z <= 0 ? 1 : 0 }
        var sig = z <= 0 ? fgo * sigl[il] / fset : fgo * sigu[il] / fset
        sig = max(sig, 0.001)
        z = z / sig
        let qx = normalTail(z)
        return z < 0 ? 1 - qx : qx
    }

    /// SELMOD: choose the sample area that controls the ray paths.
    func selmod() {
        if kfx == 1 {
            jmode = 1
        } else if kfx == 2 {
            jmode = fi[1, 1] > fi[1, 2] ? 2 : 1
        } else {
            let delfi = fi[3, 1] - fi[3, 3]
            if abs(delfi) < 0.01 {
                jmode = fi[1, 1] > fi[1, 3] ? 3 : 1
            } else {
                jmode = delfi > 0 ? 3 : 1
            }
        }
    }

    /// ANOIS1: 1 MHz atmospheric noise at the receiver for the two nearest time blocks.
    func anois1() {
        cc = gmtr
        kj = 6
        if cc - 22 < 0 { kj = fint(cc / 4 + 1) }
        tm = Float(4 * kj - 2)
        if cc - tm < 0 {
            jk = kj - 1
        } else if cc - tm == 0 {
            jk = kj
        } else {
            jk = kj + 1
        }
        if jk <= 0 {
            jk = 6
        } else if jk - 6 > 0 {
            jk = 1
        }
        let ceg = rlong < 0 ? 360 + rlong * r2d : rlongd
        let xla = rlat * r2d
        atnu = noisy(kj, xla, ceg)
        atny = noisy(jk, xla, ceg)
    }

    /// FOBBY: oblique frequencies (kHz) for each elevation angle of the reflectrix.
    func fobby(_ k: Int) {
        for iaIndex in 1...nang {
            let del = ang[iaIndex] * d2r
            let rcosd = rz * fcos(del)
            for ih in 1...30 {
                let fv = fvert[ih, k]
                let sphe = rcosd / (rz + htrue[ih, k])
                var sqcos = 1 - sphe * sphe
                sqcos = max(0.000001, sqcos)
                let cphe = fsqrt(sqcos)
                let f = fv / cphe
                freq = f
                ifob[iaIndex, ih, k] = fint(1000 * f)
            }
        }
    }

    /// ALOSFV: deviative loss factor for the high-angle part of each layer.
    func alosfv(_ k: Int) {
        var cf = alosfvCF
        var hm1 = hpmuf[1] - htmuf[1]
        var hz = hi[1, k] - yi[1, k]
        for ih in 1...10 {
            if htrue[ih, k] - htmuf[1] <= 0 {
                afac[ih, k] = 0
            } else {
                var zexp = -2 * (htrue[ih, k] - hz) / yi[1, k]
                zexp = max(zexp, -10)
                cf = 0.2 * fexp(zexp)
                afac[ih, k] = cf * (hprim[ih, k] - htrue[ih, k] - hm1)
            }
        }
        alosfvCF = cf
        if fi[2, k] <= 0 {
            // F2 layer with no F1 ledge; continuity at the E to F2 cusp.
            let a2 = cf
            hz = htrue[11, k]
            hm1 = hprim[13, k] - htrue[13, k]
            for ih in 11...12 {
                var zexp = -2 * (htrue[ih, k] - hz) / yi[3, k]
                zexp = max(zexp, -10)
                cf = a2 * fexp(zexp)
                cf = max(cf, 0.05)
                afac[ih, k] = cf * (hprim[ih, k] - htrue[ih, k] - hm1)
            }
            afac[13, k] = 0
            hm1 = hpmuf[3] - htmuf[3]
            for ih in 14...30 {
                if htrue[ih, k] - htmuf[3] <= 0 {
                    afac[ih, k] = 0
                } else {
                    var zexp = -2 * (htrue[ih, k] - hz) / yi[3, k]
                    zexp = max(zexp, -10)
                    cf = 0.1 * fexp(zexp)
                    cf = max(cf, 0.5)
                    afac[ih, k] = cf * (hprim[ih, k] - htrue[ih, k] - hm1)
                }
            }
        } else {
            // F1 layer; continuity at the E to F1 cusp.
            let a2 = cf
            hz = htrue[11, k]
            hm1 = hprim[13, k] - htrue[13, k]
            for ih in 11...12 {
                var zexp = -2 * (htrue[ih, k] - hz) / yi[2, k]
                zexp = max(zexp, -10)
                cf = a2 * fexp(zexp)
                cf = max(cf, 0.05)
                afac[ih, k] = cf * (hprim[ih, k] - htrue[ih, k] - hm1)
            }
            afac[13, k] = 0
            var htm = htmuf[2]
            var htp = hpmuf[2]
            hm1 = htp - htm
            for ih in 14...20 {
                if htrue[ih, k] - htm <= 0 {
                    afac[ih, k] = 0
                } else {
                    var zexp = -2 * (htrue[ih, k] - hz) / yi[2, k]
                    zexp = max(zexp, -10)
                    cf = 0.1 * fexp(zexp)
                    cf = max(cf, 0.05)
                    afac[ih, k] = cf * (hprim[ih, k] - htrue[ih, k] - hm1)
                }
            }
            // F2 layer above the F1 ledge.
            let a3 = cf
            hz = htrue[21, k]
            hm1 = hprim[23, k] - htrue[23, k]
            for ih in 21...22 {
                var zexp = -2 * (htrue[ih, k] - hz) / yi[3, k]
                zexp = max(zexp, -10)
                cf = a3 * fexp(zexp)
                cf = max(0.05, cf)
                afac[ih, k] = cf * (hprim[ih, k] - htrue[ih, k] - hm1)
            }
            afac[23, k] = 0
            htm = htmuf[3]
            htp = hpmuf[3]
            hm1 = htp - htm
            for ih in 24...30 {
                if htrue[ih, k] - htm <= 0 {
                    afac[ih, k] = 0
                } else {
                    var zexp = -2 * (htrue[ih, k] - hz) / yi[3, k]
                    zexp = max(zexp, -10)
                    cf = 0.1 * fexp(zexp)
                    cf = max(cf, 0.05)
                    afac[ih, k] = cf * (hprim[ih, k] - htrue[ih, k] - hm1)
                }
            }
        }
        for ih in 1...30 { afac[ih, k] = max(0, afac[ih, k]) }
    }

    /// SETLNG: for short paths, copy the sample areas so every slot is filled.
    func setlng() {
        if kfx >= 3 { return }
        func copyArea(_ source: Int, _ target: Int) {
            clat[target] = clat[source]
            clong[target] = clong[source]
            glat[target] = glat[source]
            rd[target] = rd[source]
            gyz[target] = gyz[source]
            rat[target] = rat[source]
            gmdip[target] = gmdip[source]
            clck[target] = clck[source]
            abiy[target] = abiy[source]
            artic[target] = artic[source]
            sigpat[target] = sigpat[source]
            epspat[target] = epspat[source]
            hs[target] = hs[source]
            for il in 1...3 {
                fs[il, target] = fs[il, source]
                fi[il, target] = fi[il, source]
                yi[il, target] = yi[il, source]
                hi[il, target] = hi[il, source]
            }
        }
        func copyIonogram(_ source: Int, _ target: Int) {
            for il in 1...30 {
                hprim[il, target] = hprim[il, source]
                htrue[il, target] = htrue[il, source]
                fvert[il, target] = fvert[il, source]
                afac[il, target] = afac[il, source]
                for iaIndex in 1...40 { ifob[iaIndex, il, target] = ifob[iaIndex, il, source] }
            }
        }
        if kfx <= 1 {
            for target in 2...5 { copyArea(1, target) }
            for target in 2...3 { copyIonogram(1, target) }
        } else {
            for target in 4...5 { copyArea(3, target) }
            copyIonogram(2, 3)
        }
    }

    struct SystemLoss {
        var fm: Float, su: Float, sl: Float, fmp: Float, sup: Float, slp: Float
    }

    /// SYSSY: excess system loss, its deciles and prediction errors from the tables.
    func syssy(g: Float, pt: Float, nn: Int) -> SystemLoss {
        let j = fint(pt + 0.5)
        let gg = abs(g) * 5.729577e1
        var kjIndex = fint(((gg - 40) / 5) + 1.5)
        if kjIndex <= 0 {
            kjIndex = 1
        } else if kjIndex - 9 > 0 {
            kjIndex = 9
        }
        let kk = g >= 0 ? 0 : 8
        let nd = g >= 0 ? 0 : 3
        var lj = fint(Float(j) / 3 + 0.67)
        if lj - 1 < 0 { lj = 8 }
        let lol = lj + kk
        var ld = fint(Float(j) / 6 + 0.67)
        if ld - 1 < 0 { ld = 4 }
        let fm = sys[kjIndex, lol, nn]
        let duValue = sys[kjIndex, lol, nn + 1]
        let dlValue = sys[kjIndex, lol, nn - 1]
        return SystemLoss(fm: fm, su: duValue / 1.28, sl: dlValue / 1.28,
                          fmp: perr[kjIndex, ld, nd + 1], sup: perr[kjIndex, ld, nd + 2],
                          slp: perr[kjIndex, ld, nd + 3])
    }

    /// SIGDIS: absorption parameters and the signal distribution for this path.
    func sigdis() {
        var glav: Float = 0
        let kfxx = kfx
        let xkf = Float(kfxx)
        var xadj: Float = 0, xsu: Float = 0, xsl: Float = 0, xfmp: Float = 0, xsup: Float = 0, xslp: Float = 0
        var ac: Float = 0
        feav = 0
        for k in 1...kfxx {
            abiy[k] = -0.04 + fexp(-2.937 + 0.8445 * fi[1, k])
            abiy[k] = max(0.1, abiy[k])
            let idp = gcdkm > 2500 ? 5 : 2
            let loss = syssy(g: glat[k], pt: clck[k], nn: idp)
            adj = loss.fm
            xadj = xadj + loss.fm
            xsu = xsu + loss.su
            xsl = xsl + loss.sl
            xfmp = xfmp + loss.fmp
            xsup = xsup + loss.sup
            xslp = xslp + loss.slp
            feav = feav + fi[1, k]
            ac = ac + abiy[k]
            artic[k] = adj
            glav = glav + abs(glat[k])
        }
        glav = glav / xkf
        aglat = glav
        acav = ac / xkf
        feav = feav / xkf
        adj = xadj / xkf
        su = xsu / xkf
        sl = xsl / xkf
        ads = xfmp / xkf
        sus = xsup / xkf
        sls = xslp / xkf
        let k = jmode
        htloss = 88
        xve = xlin(90, htrue, fvert, 30, k)
        xve = xve / fi[1, k]
        xnuz = 63.07
        hnu = 4.39
        if feav - 2.0 > 0 {
            afe = 1.359
            bfe = 8.617
        } else if feav - 0.5 <= 0 {
            afe = 0
            bfe = 0
        } else {
            afe = 1.359 * (feav - 0.5) / 1.5
            bfe = 8.617 * (feav - 0.5) / 1.5
        }
        // The tables are for the frequency FTAB.
        glav = glav * r2d
        var ftab: Float
        if glav - 40 <= 0 {
            ftab = yfot[3]
        } else if glav - 50 <= 0 {
            ftab = yfot[3]
            ftab = ftab - (glav - 40) * (ftab - 10) / 10
        } else {
            ftab = 10
        }
        var eslsm: Float = 0
        if ymuf[4] > 0 {
            var pes = prbmuf(ftab, ymuf[4], ymuf[4], 4)
            pes = max(0.1, pes)
            pes = min(0.9, pes)
            eslsm = -10 * alog10(1 - pes)
        }
        var pf2 = prbmuf(ftab, ymuf[3], ymuf[3], 3)
        pf2 = max(0.1, pf2)
        let f2lsm = -10 * alog10(pf2)
        asm = adj - eslsm - f2lsm
        asm = max(asm, 0)
        // Upper decile.
        var pes: Float = 0
        if yfot[4] > 0 {
            pes = prbmuf(ftab, yfot[4], yfot[4], 4)
            pes = max(0.1, pes)
            pes = min(0.9, pes)
        }
        pf2 = prbmuf(ftab, yhpf[3], yhpf[3], 3)
        pf2 = max(0.1, pf2)
        dsu = 1.28 * sl - (10 * alog10(1 - pes) + eslsm) - (10 * alog10(pf2) + f2lsm)
        dsu = max(dsu, 0.5)
        // Lower decile.
        pes = 0
        if yhpf[4] > 0 {
            pes = prbmuf(ftab, yhpf[4], yhpf[4], 4)
            pes = max(0.1, pes)
            pes = min(0.9, pes)
        }
        pf2 = prbmuf(ftab, yfot[3], yfot[3], 3)
        pf2 = max(0.1, pf2)
        dsl = 1.28 * su - (-10 * alog10(1 - pes) - eslsm) - (-10 * alog10(pf2) - f2lsm)
        dsl = max(1, dsl)
    }

    struct NoiseFamily {
        var fa: Float, dua: Float, dla: Float, dms: Float, dus: Float, dls: Float
    }

    /// GENFAM: frequency dependence of atmospheric noise, with its deciles
    /// and prediction errors.
    func genfam(y2: Float, iblk: Int, freq: Float, z: Float) -> NoiseFamily {
        var ibk = iblk
        if y2 < 0 { ibk = ibk + 6 }
        var u1: Float = -0.75
        var x = alog10(freq)
        let u = (8 * powf(2, x) - 11) / 4
        var pz: Float = 0, px: Float = 0, cz: Float = 0
        for pass in 1...2 {
            pz = u1 * fam[1, ibk] + fam[2, ibk]
            px = u1 * fam[8, ibk] + fam[9, ibk]
            for i in 3...7 {
                pz = u1 * pz + fam[i, ibk]
                px = u1 * px + fam[i + 7, ibk]
            }
            if pass == 1 {
                cz = z * pz + px
                cz = z + z - cz
                u1 = u
            }
        }
        let fa = cz * pz + px
        if freq > 20 { x = alog10(20) }
        var v = A1(5)
        for i in 1...5 {
            if i == 5 && freq > 10 { x = 1 }
            var y = dud[1, ibk, i]
            for j in 2...5 { y = y * x + dud[j, ibk, i] }
            v[i] = y
        }
        return NoiseFamily(fa: fa, dua: v[1], dla: v[2], dms: v[5], dus: v[3], dls: v[4])
    }

    /// GAIN for an antenna (`itr` 1 transmitter, 2 receiver): gain in dBi at
    /// elevation `delta` radians, and the antenna efficiency.
    func gain(_ itr: Int, _ delta: Float, _ fmc: Float) -> (rain: Float, eff: Float) {
        let antenna = itr == 1 ? txAntenna : rxAntenna
        guard antenna.minFrequency <= fmc && fmc <= antenna.maxFrequency else { return (0, 0) }
        let i = fint(fmc)
        var ip1 = i + 1
        ip1 = min(30, ip1)
        let xfmc = fmc - Float(i)
        let deltd = delta * r2d
        var j = fint(deltd + 1)
        j = min(j, 91)
        var jp1 = j + 1
        jp1 = min(91, jp1)
        let xdelta = deltd - Float(j - 1)
        let xx = antenna.array(i, j)
        let yx = antenna.array(ip1, j)
        let xy = antenna.array(i, jp1)
        let yy = antenna.array(ip1, jp1)
        let rx = xx + xfmc * (yx - xx)
        let ry = xy + xfmc * (yy - xy)
        let rain = rx + xdelta * (ry - rx)
        let effValue = antenna.efficiency[i - 1] + xfmc * (antenna.efficiency[ip1 - 1] - antenna.efficiency[i - 1])
        return (rain, effValue)
    }

    /// GAIN with a negative first argument: ground reflection loss at sample area `ig`.
    func groundLoss(_ ig: Int, _ delta: Float, _ fmc: Float) -> Float {
        let sigma = sigpat[ig]
        let er = epspat[ig]
        if sigma <= 0 { return 0 }
        let x = 18000 * sigma / fmc
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
        let cv = fsqrt(rho * rho + u * u * s - 2 * rho * u * r * fcos(alpha + 2 * asxv))
            / (rho + u * r + 2 * rho12 * v * q * fcos(alpha * 0.5 + asxv))
        let ch = fsqrt(rho * rho + s - 2 * rho * r * fcos(alpha))
            / (rho + r + 2 * rho12 * q * fcos(alpha * 0.5))
        var rain = 4.3429 * alog(0.5 * (ch * ch + cv * cv))
        rain = abs(rain)
        if delta <= 0.00000001 { rain = 6 }
        return rain
    }

    /// PWRDB: transmitter power in dBW at a frequency.
    func pwrdb(_ freq: Float) -> Float {
        var f = freq
        if f < 2 { f = 2 }
        if f > 30 { f = 30 }
        if txAntenna.minFrequency <= f && f <= txAntenna.maxFrequency { return pwrdba }
        return 30
    }

    /// ALLMODES: reset (flag 0) or append the current hop's modes to the
    /// list for this frequency.
    func allmodes(_ iflg: Int) {
        var ist = 1, lst = 6
        if iflg == 0 {
            nmmod = 0
            for i in 1...20 {
                ztloss[i] = 99999
                ztllow[i] = 999
                ztlhgh[i] = 999
                zhp[i] = -1
            }
            return
        } else if iflg == 1 {
            ist = 1
            hp[2] = hp[1]
            ztloss[2] = 32000
            lst = 1
            iznmode[nmmod + 2] = nmode[2]
            zb[nmmod + 2] = b[2]
        } else if iflg == 2 {
            ist = 4
            lst = 5
        }
        for i in ist...lst where hp[i] > 0 {
            nmmod = nmmod + 1
            ztloss[nmmod] = tloss[i]
            ztllow[nmmod] = tllow[i]
            ztlhgh[nmmod] = tlhgh[i]
            zhp[nmmod] = hp[i]
            zcrel[nmmod] = crel[i]
            zrely[nmmod] = rely[i]
            zhn[nmmod] = hn[i]
            iznmode[nmmod] = nmode[i]
            zsn[nmmod] = sn[i]
            zfldst[nmmod] = fldst[i]
            zsigpow[nmmod] = sigpow[i]
            zb[nmmod] = b[i]
            ztimed[nmmod] = timed[i]
            zabps[nmmod] = abps[i]
            zprob[nmmod] = prob[i]
            zrgain[nmmod] = rgain[i]
            ztgain[nmmod] = tgain[i]
            zfslos[nmmod] = fslos[i]
            zspro[nmmod] = spro[i]
            zeff[nmmod] = eff[i]
            zgrlos[nmmod] = grlos[i]
            zadv[nmmod] = adv[i]
            zobf[nmmod] = obf[i]
        }
    }
}
