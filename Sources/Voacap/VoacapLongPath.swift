import Foundation

// The long-path model: SETTXR, GMLOSS, SELTXR, LNGPAT, CONVH, GETTOP, TABS
// and BABS.

extension VoacapEngine {
    /// SETTXR: losses and antenna gain for every ray at the transmitter and receiver ends.
    func settxr() {
        for jj in 1...2 {
            let k = itxrcp[jj]
            let dend = min(4000, gcdkm)
            for i in 1...45 {
                gml[i, k] = -999
                fhp[i, k] = 999
            }
            var found = false
            for i in 1...45 {
                if hpflx[i, k] - 70 < 0 { break }
                let xhpm = dend / gdflx[i, k]
                if xhpm - 0.9 < 0 { continue }
                if delfx[i, k] - amind < 0 { continue }
                found = true
                break
            }
            if !found {
                // The layer was penetrated: include the over-the-MUF mode.
                let im = modmuf
                delfx[1, k] = delmuf[im]
                hpflx[1, k] = hpmuf[im]
                htflx[1, k] = htmuf[im]
                fvflx[1, k] = fvmuf[im]
                afflx[1, k] = afmuf[im]
                let xhop = Float(nhopmf[im])
                gdflx[1, k] = gcdkm / xhop
                imode[1, k] = modmuf
            }
            let bc = powf(freq + gyz[k], 1.98)
            let ac = 677.2 * abiy[k]
            for i in 1...45 {
                tgainx[i, k] = -10
                andvx[i, k] = 1000
                advx[i, k] = 1000
                grlosx[i, k] = 0
                if hpflx[i, k] < 70 { continue }
                if delfx[i, k] < amind { continue }
                let del = delfx[i, k] * d2r
                let cdel = fcos(del)
                if !(fi[1, k] < fvflx[i, k]) {
                    // D–E layer modes.
                    var xnsq: Float
                    if htloss <= htflx[i, k] {
                        xnsq = 10.2
                    } else {
                        let hnux = 61 + 3 * (htflx[i, k] - 70) / 18
                        xnsq = xnuz * fexp(-2 * (hnux - 60) / hnu)
                    }
                    let heff = min(100, htflx[i, k])
                    var sinp = rz * cdel / (rz + heff)
                    var secp = 1 / fsqrt(1 - sinp * sinp)
                    andvx[i, k] = secp * ac / (bc + xnsq)
                    let xv = max(fvflx[i, k] / fi[1, k], xve)
                    let adx = afe + bfe * alog(xv)
                    sinp = rz * cdel / (rz + hpflx[i, k])
                    secp = 1 / fsqrt(1 - sinp * sinp)
                    advx[i, k] = secp * afflx[i, k] * (powf(fvflx[i, k] + gyz[k], 1.98) + xnsq) / (bc + xnsq) + adx
                    aofx[i, k] = 0
                } else {
                    // F layer modes.
                    let xnsq: Float = 10.2
                    var sinp = rz * cdel / (rz + 100)
                    var secp = 1 / fsqrt(1 - sinp * sinp)
                    andvx[i, k] = secp * ac / (bc + xnsq)
                    sinp = rz * cdel / (rz + hpflx[i, k])
                    secp = 1 / fsqrt(1 - sinp * sinp)
                    advx[i, k] = secp * afflx[i, k] * (powf(fvflx[i, k] + gyz[k], 1.98) + xnsq) / (bc + xnsq)
                    aofx[i, k] = 0
                    if fs[2, k] > 0 {
                        var fsdead = Float(ifob[1, 3, k])
                        fsdead = fsdead / 1000
                        let fmhz = max(freq, fsdead)
                        let sins = rz * cdel / (rz + hs[k])
                        let secs = 1 / fsqrt(1 - sins * sins)
                        let esd = fs[2, k] * secs
                        var pValue = prbmuf(fmhz, esd, ymuf[4], 4)
                        pValue = max(0.1, pValue)
                        pValue = min(0.9, pValue)
                        aofx[i, k] = -10 * alog10(1 - pValue)
                    }
                }
                var y: Float = 0
                for ig in 1...km { y = y + groundLoss(ig, del, freq) }
                grlosx[i, k] = y / Float(km)
                let antenna = gain(jj, del, freq)
                tgainx[i, k] = antenna.rain
                if jj == 2 { efflp[i] = antenna.eff }
                let sphet = rz * cdel / (rz + htflx[i, k])
                var cphet = 1 - sphet * sphet
                cphet = max(0.00000001, cphet)
                cphet = fsqrt(cphet)
                let isLayer = modmuf
                let pros = prbmuf(freq, ymuf[isLayer], ymuf[isLayer], isLayer)
                // Over-the-MUF loss, median.
                let xls = -10 * alog10(pros) / cphet
                andvx[i, k] = andvx[i, k] + xls
                zprob[jj] = pros
            }
        }
    }

    /// GMLOSS: gain minus loss for each ray at both ends.
    func gmloss() {
        for im in 1...7 {
            obf[im] = 1000
            adv[im] = 1000
            fslos[im] = 1000
            tloss[im] = 1000
            abps[im] = 1000
            eff[im] = 0
            grlos[im] = 1000
            rgain[im] = 0
            tgain[im] = 0
            hn[im] = -1
            prob[im] = 0.001
            crel[im] = 0.001
            rely[im] = 0.001
            spro[im] = 0.001
            fldst[im] = -1000
            sigpow[im] = -1000
            sn[im] = -1000
            timed[im] = -1
            hp[im] = -1
            b[im] = -1
            nmode[im] = 5
        }
        settxr()
        for jj in 1...2 {
            let k = itxrcp[jj]
            if iaftxr[k] < 1 { continue }
            for i in 1...iaftxr[k] where hpflx[i, k] - 70 >= 0 && delfx[i, k] - amind >= 0 {
                txrgml[i, jj] = tgainx[i, k] - (andvx[i, k] + advx[i, k])
            }
        }
    }

    /// SELTXR: choose the best take-off ray at each end.
    func seltxr() {
        for jj in 1...2 {
            let k = itxrcp[jj]
            let dend = min(gcdkm, 4000)
            ltxrgm[jj] = 1
            // Find the first good mode.
            var xhpm = dend / gdflx[ltxrgm[jj], k]
            let ihopFirst = fint(xhpm)
            var abandoned = false
            if xhpm - 0.9 < 0 {
                while true {
                    ltxrgm[jj] = ltxrgm[jj] + 1
                    if ltxrgm[jj] > 45 || hpflx[ltxrgm[jj], k] <= 70 || ltxrgm[jj] - iaftxr[k] > 0 {
                        ltxrgm[jj] = ltxrgm[jj] - 1
                        abandoned = true
                        break
                    }
                    if delfx[ltxrgm[jj], k] - amind >= 0 { break }
                }
            }
            if abandoned { continue }
            let hopx = Float(ihopFirst)
            var yhpm = xhpm - hopx
            if yhpm - 0.5 > 0 { yhpm = 1 - yhpm }
            var gmax = tgainx[ltxrgm[jj], k] - xhpm * (andvx[ltxrgm[jj], k] + advx[ltxrgm[jj], k])
            var delmax = abs(delfx[ltxrgm[jj], k] - delopt)
            let ls = ltxrgm[jj]
            if ls > iaftxr[k] { continue }
            for i in ls...iaftxr[k] {
                if hpflx[i, k] - 70 <= 0 { continue }
                let xhop = dend / gdflx[i, k]
                let ihop = fint(xhop)
                if xhop - 0.9 < 0 { continue }
                var yhop = xhop - Float(ihop)
                if yhop - 0.5 > 0 { yhop = 1 - yhop }
                let gnow = tgainx[i, k] - xhop * (andvx[i, k] + advx[i, k])
                gml[i, k] = gnow
                fhp[i, k] = yhop
                let delnow = abs(delfx[i, k] - delopt)
                var select = false
                if abs(gnow - gmax) - gmin <= 0 {
                    // Close in gain minus loss: prefer a whole number of hops, then the best take-off angle.
                    if abs(yhop - yhpm) - ymin <= 0 {
                        select = delnow - delmax < 0
                    } else {
                        select = yhop - yhpm < 0
                    }
                } else {
                    select = gnow - gmax > 0
                }
                if select {
                    ltxrgm[jj] = i
                    xhpm = xhop
                    gmax = gnow
                    delmax = delnow
                    yhpm = yhop
                }
            }
        }
    }

    /// CONVH: geometrical convergence factor and total group path.
    private func convh(gd: Float, phe: Float, del: Float, hp: Float, ray: Float) -> (ch: Float, ptot: Float) {
        let denom = max(abs(fsin(gd)), 0.000001)
        let psi = pio2 - del - phe
        var gm = gd - 2 * psi
        gm = max(gm, 0.001)
        let ptot = 2 * ray + (rz + hp) * gm
        let smallc = (ptot / rz) * fcos(del) / denom
        var ch = 10 * alog10(smallc)
        ch = min(ch, 15)
        return (ch, ptot)
    }

    /// GETTOP: the distance that supports modes travelling above the E layer.
    private func gettop(dikm: Float) -> Float {
        var fpe = A1(3)
        let del = delopt * d2r
        let rcosd = rz * fcos(del)
        for jf in 1...3 {
            let ht = hi[1, jf]
            let fv = fi[1, jf]
            var sphe = rcosd / (rz + ht)
            sphe = sphe * sphe
            sphe = min(sphe, 0.9999)
            fpe[jf] = fv / fsqrt(1 - sphe)
        }
        let fp1 = max(fpe[1], fpe[3])
        let fp2 = fpe[2]
        let fp3 = min(fpe[1], fpe[3])
        if fp1 - fp2 >= 0 { return 0 }
        if freq - fp2 - 0.001 >= 0 { return 0 }
        if freq - fp1 + 0.001 <= 0 { return 0 }
        let gmkm = 0.5 * dikm * (fp2 - freq) * (1 / (fp2 - fp1) + 1 / (fp2 - fp3))
        if gmkm - 1000 <= 0 { return 0 }
        if dikm - gmkm - 1000 <= 0 { return 0 }
        return gmkm
    }

    /// TABS (`absorbing` false) and BABS (`absorbing` true): loss per
    /// kilometre for the parts of the path that miss or pass through the
    /// absorbing region.
    private func pathAbsorption(del: Float, hp: Float, distance: Float, absorbing: Bool) -> Float {
        if distance <= 0 { return 1000 }
        if absorbing ? hp - 70 < 0 : hp - 70 <= 0 { return 1000 }
        var af: Float = absorbing ? acav : 0.1
        let abValue = 10.2 + powf(freq + gyz[3], 1.98)
        af = 677.2 * af / abValue
        let rcosd = rz * fcos(del)
        var sphe = rcosd / (rz + hp)
        let phe = fatan(sphe / fsqrt(1 - sphe * sphe))
        let psi = pio2 - phe - del
        let xmhop = 0.5 * distance / psi
        sphe = rcosd / (rz + 100)
        let scphe = 1 / fsqrt(1 - sphe * sphe)
        let port = xmhop * scphe / (distance * rz)
        return af * port
    }

    /// LNGPAT: put the long-path mode together from the two ends.
    func lngpat() {
        if ltxrgm[1] <= 0 || ltxrgm[2] <= 0 { return }
        let l1 = ltxrgm[1], l2 = ltxrgm[2], k2 = itxrcp[2]
        let del = 0.5 * d2r * (delfx[l1, 1] + delfx[l2, k2])
        let hpx = 0.5 * (hpflx[l1, 1] + hpflx[l2, k2])
        let rcosd = rz * fcos(del)
        let sphe = rcosd / (rz + hpx)
        let phe = fatan(sphe / fsqrt(1 - sphe * sphe))
        let dt = 0.5 * gdflx[l1, 1] / rz
        let dr = 0.5 * gdflx[l2, k2] / rz
        let ray = rz * fcos(phe + del) / sphe
        let convergence = convh(gd: gcd, phe: phe, del: del, hp: hpx, ray: ray)
        let ptot = convergence.ptot
        let free = 36.58 + 20 * alog10(0.6214 * ptot * freq) - convergence.ch
        // Distance travelled in the ionosphere.
        let di = gcd - dt - dr
        let dikm = di * rz
        var dm: Float = 0, am: Float = 0, af: Float = 0, df: Float = 0
        if !(dikm <= -1) {
            let gmkm = gettop(dikm: dikm)
            dm = gmkm / rz
            df = di - dm
            df = max(df, 0)
            am = pathAbsorption(del: del, hp: hpx, distance: dm, absorbing: false)
            af = pathAbsorption(del: del, hp: hpx, distance: df, absorbing: true)
        }
        let xlm = am * dm * rz
        let xlf = af * df * rz
        let xlt = (andvx[l1, 1] + advx[l1, 1]) / 2
        let xlr = (andvx[l2, k2] + advx[l2, k2]) / 2
        let gloss = 0.5 * (grlosx[l1, 1] + grlosx[l2, k2])
        var fhop = df / (dr + dt)
        fhop = max(fhop, 1)
        nmmod = 1
        let tlosl = free + xlt + xlm + xlf + xlr + (fhop - 1) * gloss + asm - tgainx[l2, k2] - tgainx[l1, 1]
        zgrlos[1] = gloss
        zhn[1] = 0.5 * gcd / (dr + dt)
        zhn[1] = max(1.0, zhn[1])
        zhp[1] = hpx
        zrgain[1] = tgainx[l2, k2]
        ztgain[1] = tgainx[l1, 1]
        zeff[1] = efflp[l2]
        ztimed[1] = ptot / vofl
        ztloss[1] = tlosl
        zb[1] = delfx[l1, 1]
        zb[2] = delfx[l2, k2]
        iznmode[1] = imode[l1, 1]
        iznmode[2] = imode[l2, k2]
        zfslos[1] = free
        zadv[1] = (advx[l1, 1] + advx[l2, k2]) * 0.5
        zabps[1] = (andvx[l1, 1] + andvx[l2, k2]) * 0.5
        zobf[1] = (aofx[l1, 1] + aofx[l2, k2]) * 0.5
        zsigpow[1] = pwrdb(freq) - tlosl
        zfldst[1] = 107.2 + pwrdb(freq) + 20 * alog10(freq) - tlosl - zrgain[1]
        zsn[1] = zsigpow[1] - rcnse - zeff[1]
        // Signal deciles from the MUF distribution at each end.
        func deciles(_ isLayer: Int) -> (Float, Float) {
            var pValue = prbmuf(freq, ymuf[isLayer], ymuf[isLayer], isLayer)
            let zlmuf = -10 * alog10(pValue)
            pValue = prbmuf(freq, yfot[isLayer], yfot[isLayer], isLayer)
            let zlfot = -10 * alog10(pValue)
            pValue = prbmuf(freq, yhpf[isLayer], yhpf[isLayer], isLayer)
            let zlhpf = -10 * alog10(pValue)
            let cpr = fvmuf[isLayer] / ymuf[isLayer]
            return ((zlfot - zlmuf) / cpr, (zlmuf - zlhpf) / cpr)
        }
        let first = deciles(iznmode[1])
        dslf = first.0
        dsuf = first.1
        if iznmode[1] != iznmode[2] {
            let second = deciles(iznmode[2])
            dsuf = 0.5 * (second.1 + first.1)
            dslf = 0.5 * (second.0 + first.0)
        }
        ztllow[1] = dsl + zhn[1] * dslf
        ztlhgh[1] = dsu + zhn[1] * dsuf
        ztllow[1] = min(ztllow[1], 25)
        ztlhgh[1] = min(ztlhgh[1], 25)
        zprob[1] = min(zprob[1], zprob[2])
    }
}
