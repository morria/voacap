import Foundation

// Sporadic E modes and the combination of all modes into one result per
// frequency: ESMOD, ESREG and RELBIL.

extension VoacapEngine {
    /// Reset one of the spare mode slots.
    private func clearMode(_ i: Int) {
        obf[i] = 1000
        adv[i] = 1000
        fslos[i] = 1000
        tloss[i] = 1000
        abps[i] = 1000
        eff[i] = 0
        grlos[i] = 1000
        rgain[i] = 0
        tgain[i] = 0
        hn[i] = -1
        prob[i] = 0.001
        rely[i] = 0.001
        spro[i] = 0.001
        fldst[i] = -1000
        sigpow[i] = -1000
        sn[i] = -1000
        hp[i] = -1
        b[i] = -1
        nmode[i] = 5
        tllow[i] = 10
        timed[i] = -1
        tlhgh[i] = 10
    }

    /// ESMOD: up to two sporadic E modes (slots 4 and 5).
    func esmod() {
        var k = 1
        for isIndex in 1...km where fs[2, k] - fs[2, isIndex] > 0 { k = isIndex }
        var fsdead = Float(ifob[1, 3, jmode])
        fsdead = fsdead / 1000
        fsdead = min(fsdead, 3)
        for iht in 4...5 {
            clearMode(iht)
            crel[iht] = -1000
        }
        if freq <= fsdead { return }
        if fs[2, k] <= 0 { return }
        let sdmax = 2 * rz * (pio2 - fasin(1 / (1 + hs[k] / rz)))
        let ihsrt = fint(gcdkm / sdmax + 1)
        let ihstp = 2
        if ihsrt > 2 { return }
        var ih = 3
        let ac = 677.2 * acav / (powf(freq + gyz[k], 1.98) + 10.2)
        for ihop in ihsrt...ihstp {
            ih = ih + 1
            let gp = Float(ihop)
            ghop = gcd / gp
            let thet = 0.5 * ghop
            let tans = fsin(thet) / (1 - fcos(thet) + hs[k] / rz)
            let psi = fatan(tans)
            let secs = 1 / fcos(psi)
            let sfvmod = freq / secs
            var esd = fs[2, k] * secs
            let del = pio2 - thet - psi
            let cdel = fcos(del)
            let adel = del * r2d
            if adel < amind { continue }
            var path = 2 * fsin(0.5 * ghop) * (rz + hs[k]) / cdel
            let hop = Float(ihop)
            path = hop * path
            let sflos = 32.45 + 20 * alog10(path * freq)
            let sinp = rz * cdel / (rz + 100)
            let cosp = fsqrt(1 - sinp * sinp)
            var sabps = ac / cosp
            var adx: Float = 0
            // Below the E layer critical frequency, remove the E layer bending above the Es height.
            if fi[1, k] > sfvmod { adx = afe + bfe * alog(sfvmod / fi[1, k]) }
            sabps = sabps + adx
            var pros = prbmuf(freq, esd, ymuf[4], 4)
            pros = min(pros, 0.90)
            let ref = 8.9136 * powf(pros, -0.7)
            esd = fs[1, k] * secs
            var ps = prbmuf(freq, esd, yfot[4], 4)
            ps = min(ps, 0.9)
            let refl = 8.9136 * powf(ps, -0.7)
            esd = fs[3, k] * secs
            ps = prbmuf(freq, esd, yhpf[4], 4)
            ps = min(ps, 0.9)
            let refu = 8.9136 * powf(ps, -0.7)
            tllow[ih] = dsl + hop * (refl - ref)
            tlhgh[ih] = dsu + hop * (ref - refu)
            tllow[ih] = min(tllow[ih], 25)
            tlhgh[ih] = min(tlhgh[ih], 25)
            var y: Float = 0
            for ig in 1...km { y = y + groundLoss(ig, del, freq) }
            let sgrlos = y / Float(km)
            let stgain = gain(1, del, freq).rain
            let receive = gain(2, del, freq)
            let srgain = receive.rain
            eff[ih] = receive.eff
            let xtlos = sflos + hop * (sabps + ref + adx) + (hop - 1) * sgrlos - srgain - stgain + asm
            fldst[ih] = 107.2 + pwrdb(freq) + 20 * alog10(freq) - xtlos - srgain
            sigpow[ih] = pwrdb(freq) - xtlos
            pros = prbmuf(freq, ymuf[4], ymuf[4], 4)
            obf[ih] = 8.9136 * powf(pros, -0.7)
            adv[ih] = 0
            fslos[ih] = sflos
            tloss[ih] = xtlos
            abps[ih] = ac / cosp + adx
            grlos[ih] = sgrlos
            rgain[ih] = srgain
            tgain[ih] = stgain
            hn[ih] = hop
            sn[ih] = sigpow[ih] - rcnse - eff[ih]
            hp[ih] = hs[k]
            b[ih] = adel
            nmode[ih] = 4
            prob[ih] = powf(pros, Float(ihop))
            timed[ih] = path / vofl
        }
    }

    /// ESREG: the mixed Es–F mode. VOACAP presets the slots and returns
    /// without computing one, and so does this.
    func esreg() {
        for im in 6...7 {
            clearMode(im)
            crel[im] = 1000
        }
    }

    /// RELBIL: reliability of each mode, the most reliable mode, the
    /// combined signal from all modes and the reliability of the sum.
    func relbil(_ ifx: Int) {
        let tme: [Float] = [0.0, 0.1257, 0.2533, 0.3853, 0.5244, 0.6745, 0.8416, 1.0364, 1.2815, 1.6449]
        let xeps: Float = 0.05
        let inum = nmmod
        if inum <= 0 { return }
        let du2 = du * du
        let dl2 = dl * dl
        for im in 1...inum {
            if zhp[im] <= 70 {
                zcrel[im] = 0.001
                zrely[im] = 0.001
            } else {
                dslf = ztllow[im]
                dsuf = ztlhgh[im]
                d10r = fsqrt(du2 + dslf * dslf)
                d50r = zsn[im]
                d90r = fsqrt(dl2 + dsuf * dsuf)
                var z = rsn - d50r
                z = z <= 0 ? z / (d10r / 1.28) : z / (d90r / 1.28)
                zrely[im] = 1 - fnorml(z)
                zcrel[im] = 1000
            }
        }
        // Most reliable mode.
        var irmax = 1
        if inum >= 2 {
            for im in 2...inum where zrely[im] > zrely[irmax] { irmax = im }
        }
        var ir = irmax
        var xrel = zrely[ir]
        var xhn = zhn[ir]
        var xsn = zsn[ir]
        for im in 1...inum {
            if im == irmax { continue }
            // Reliability first, but when close prefer fewer hops, then the better median SNR.
            if xrel < 0.00000001 { xrel = 0.00000001 }
            let reltest = abs(zrely[im] - xrel)
            if reltest >= xeps { continue }
            var better = false
            if abs(xhn - zhn[im]) <= xeps {
                if xsn < zsn[im] { better = true }
            } else if xhn > zhn[im] {
                better = true
            }
            if better {
                ir = im
                xhn = zhn[im]
                xsn = zsn[im]
                xrel = zrely[im]
            }
        }
        nrel = ir
        let isLayer = iznmode[ir]
        if inum == 1 {
            reliab[ifx] = zrely[ir]
            dblosl[ifx] = ztllow[ir]
            dblosu[ifx] = ztlhgh[ir]
            dbu[ifx] = zfldst[ir]
            sndb[ifx] = zsn[ir]
            dbw[ifx] = zsigpow[ir]
        } else {
            // Add the signals with random phase: sum the powers in watts.
            var xdslw: Float = 0, xsigs: Float = 0, xdsup: Float = 0, xfld: Float = 0
            var dxsigs: Float = -1000, dxfld: Float = -1000, dxdslw: Float = -1000, dxdsup: Float = -1000
            for iv in 1...inum {
                dxsigs = max(dxsigs, zsigpow[iv])
                dxfld = max(dxfld, zfldst[iv])
                dxdslw = max(dxdslw, zsigpow[iv] - ztllow[iv])
                dxdsup = max(dxdsup, zsigpow[iv] + ztlhgh[iv])
            }
            for im in 1...inum {
                var zexp = 0.1 * (zsigpow[im] - ztllow[im] - dxdslw)
                if zexp > -10 { xdslw = xdslw + powf(10, zexp) }
                zexp = 0.1 * (zsigpow[im] - dxsigs)
                if zexp > -10 { xsigs = xsigs + powf(10, zexp) }
                zexp = 0.1 * (zsigpow[im] + ztlhgh[im] - dxdsup)
                if zexp > -10 { xdsup = xdsup + powf(10, zexp) }
                zexp = 0.1 * (zfldst[im] - dxfld)
                if zexp > -10 { xfld = xfld + powf(10, zexp) }
            }
            var sigmed: Float = -500
            if xsigs > 0 { sigmed = dxsigs + 10 * alog10(xsigs) }
            dblosl[ifx] = xdslw > 0 ? abs(sigmed - 10 * alog10(xdslw) - dxdslw) : 0
            if dblosl[ifx] < 0.2 { dblosl[ifx] = 0.2 }
            if dblosl[ifx] > 30 { dblosl[ifx] = 30 }
            dblosu[ifx] = xdsup > 0 ? abs(dxdsup + 10 * alog10(xdsup) - sigmed) : 0
            if dblosu[ifx] < 0.2 { dblosu[ifx] = 0.2 }
            if dblosu[ifx] > 30 { dblosu[ifx] = 30 }
            dbw[ifx] = sigmed
            let delsig = sigmed - zsigpow[ir]
            dbu[ifx] = xfld > 0 ? dxfld + 10 * alog10(xfld) : -500
            sndb[ifx] = zsn[ir] + delsig
            // Reliability of the sum of the modes.
            d10r = fsqrt(du2 + dblosl[ifx] * dblosl[ifx])
            d50r = sndb[ifx]
            d90r = fsqrt(dl2 + dblosu[ifx] * dblosu[ifx])
            if d10r < 0.2 { d10r = 0.2 }
            if d90r > 30 { d90r = 30 }
            var z = rsn - d50r
            z = z <= 0 ? z / (d10r / 1.28) : z / (d90r / 1.28)
            reliab[ifx] = 1 - fnorml(z)
        }
        gaint[ifx] = ztgain[ir]
        gainr[ifx] = zrgain[ir]
        snrlw[ifx] = d10r
        snrup[ifx] = d90r
        angle[ifx] = zb[ir]
        vhigh[ifx] = zhp[ir]
        delay[ifx] = ztimed[ir]
        // The signal powers were summed, so the transmission loss follows from the sum.
        dblos[ifx] = pwrdb(freq) - dbw[ifx]
        cprob[ifx] = zprob[ir]
        mode[ifx] = isLayer
        nhp[ifx] = fint(zhn[ir])
        xnynois[ifx] = rcnse
        rneff[ifx] = zeff[ir]
        // Power gain needed for the required reliability.
        var itm = abs(lufp - 50) / 5 + 1
        itm = min(itm, 10)
        let tmx = tme[itm - 1] / tme[8]
        if lufp < 50 {
            snpr[ifx] = -(d50r + tmx * d90r) + rsn
        } else {
            snpr[ifx] = -(d50r - tmx * d10r) + rsn
        }
        snxx[ifx] = rsn - snpr[ifx]
    }
}
