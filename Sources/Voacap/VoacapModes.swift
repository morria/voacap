import Foundation

// Modes for one hop count on the short-path model: FDIST, INMUF and REGMOD.

extension VoacapEngine {
    /// FDIST: find up to six ray sets in the reflectrix for the current hop distance.
    func fdist(_ k: Int) {
        enum Label { case l140, l145, l165, l176, l180, l305, done }
        let dhopkm = ghop * rz
        for im in 1...6 {
            hpmod[im, k] = -1
            itmod[im, k] = 5
        }
        if dhopkm >= dmaxkm[k] { return }
        var ih = 0
        var il = 0
        var label = Label.l140
        loop: while true {
            switch label {
            case .l140:
                il = il + 1
                label = .l145
            case .l145:
                ih = ih + 1
                if ih > 44 { label = .done; continue loop }
                if hpflx[ih + 1, k] <= 0 { label = .done; continue loop }
                let difference = gdflx[ih, k] - gdflx[ih + 1, k]
                if difference < 0 {
                    let e = gdflx[ih, k] - dhopkm
                    if e > 0 { continue loop }
                    if e == 0 { label = .l165; continue loop }
                    let f = dhopkm - gdflx[ih + 1, k]
                    if f > 0 { continue loop }
                    label = f == 0 ? .l176 : .l180
                } else if difference == 0 {
                    if abs(dhopkm - gdflx[ih, k]) <= 0.001 { continue loop }
                    label = .l165
                } else {
                    let e = gdflx[ih, k] - dhopkm
                    if e < 0 { continue loop }
                    if e == 0 { label = .l165; continue loop }
                    let f = dhopkm - gdflx[ih + 1, k]
                    if f < 0 { continue loop }
                    label = f == 0 ? .l176 : .l180
                }
            case .l176:
                ih = ih + 1
                label = .l165
            case .l165:
                delmod[il, k] = delfx[ih, k]
                hpmod[il, k] = hpflx[ih, k]
                htmod[il, k] = htflx[ih, k]
                itmod[il, k] = imode[ih, k]
                afmod[il, k] = afflx[ih, k]
                fvmod[il, k] = fvflx[ih, k]
                label = .l305
            case .l180:
                let dth2 = gdflx[ih + 1, k] - gdflx[ih, k]
                if abs(dth2) - 0.001 <= 0 { label = .l165; continue loop }
                let dth = dhopkm - gdflx[ih, k]
                // Linear interpolation in the table.
                let thet = 0.5 * dhopkm / rz
                let hp1 = hpflx[ih, k]
                let hp2 = hpflx[ih + 1, k]
                let ht1 = htflx[ih, k]
                let ht2 = htflx[ih + 1, k]
                let hpValue = hp1 + (hp2 - hp1) * dth / dth2
                let ht = ht1 + (ht2 - ht1) * dth / dth2
                afmod[il, k] = afflx[ih, k] + (afflx[ih + 1, k] - afflx[ih, k]) * dth / dth2
                hpmod[il, k] = hpValue
                let st = fsin(thet)
                let tanp = st / (1 - fcos(thet) + hpValue / rz)
                let phe = fatan(tanp)
                // Force correct geometry: the radiation angle, then Snell's law for FV.
                delmod[il, k] = (pio2 - phe - thet) * r2d
                var sphi = rz * fcos(delmod[il, k] * d2r) / (rz + ht)
                sphi = 1 - sphi * sphi
                sphi = max(sphi, 0.000001)
                fvmod[il, k] = freq * fsqrt(sphi)
                htmod[il, k] = ht
                itmod[il, k] = imode[ih, k]
                label = .l305
            case .l305:
                if delmod[il, k] >= amind {
                    label = il - 6 < 0 ? .l140 : .done
                } else {
                    hpmod[il, k] = -abs(hpmod[il, k])
                    label = .l145
                }
            case .done:
                break loop
            }
        }
    }

    /// INMUF: add over-the-MUF or zero-distance modes so every layer has a
    /// candidate, then price all the modes with REGMOD.
    func inmuf(_ ihop: Int) {
        let eps: Float = 0.4001
        let k = jmode
        var ireset = false
        var inm = I1(3, -1)
        var osigl = A1(3), osigu = A1(3), oyfot = A1(3), oyhpf = A1(3), oymuf = A1(3)
        var inum = 0
        for im in 1...6 where hpmod[im, k] - 70 >= 0 {
            inum = inum + 1
            inm[itmod[im, k]] = 1
            delmod[inum, k] = delmod[im, k]
            hpmod[inum, k] = hpmod[im, k]
            htmod[inum, k] = htmod[im, k]
            fvmod[inum, k] = fvmod[im, k]
            afmod[inum, k] = afmod[im, k]
            itmod[inum, k] = itmod[im, k]
        }

        var zeroDistance = delmuf[modmuf] >= 89.9
        var finished = false
        if !zeroDistance {
            for ily in 1...3 {
                if ily == 2 && fi[2, k] <= 0 { continue }
                if inum >= 5 { break }
                if inm[ily] < 0 && (freq + eps) >= ymuf[ily] && ihop == nhopmf[ily] {
                    // Over the MUF: insert this layer's MUF mode.
                    inum = inum + 1
                    delmod[inum, k] = delmuf[ily]
                    hpmod[inum, k] = hpmuf[ily]
                    htmod[inum, k] = htmuf[ily]
                    fvmod[inum, k] = fvmuf[ily]
                    afmod[inum, k] = afmuf[ily]
                    itmod[inum, k] = ily
                    inm[ily] = 1
                }
            }
            if ihop == nhopmf[modmuf] {
                if inum > 0 {
                    finished = true
                } else {
                    zeroDistance = true
                }
            }
        }
        if !zeroDistance && !finished {
            // Find a new MUF for this number of hops.
            ireset = true
            for il in 1...3 {
                osigl[il] = sigl[il]
                osigu[il] = sigu[il]
                oyfot[il] = yfot[il]
                oyhpf[il] = yhpf[il]
                oymuf[il] = ymuf[il]
            }
            for jh in 1...3 {
                if inum >= 5 { break }
                if inm[jh] > 0 { continue }
                if jh == 2 && fi[2, k] <= 0 { continue }
                if ihop <= nhopmf[jh] { continue }
                let fv = fvmuf[jh]
                let ht = htmuf[jh]
                let hpver = xlin(fv, fvert, hprim, 30, k)
                let xhp = (hpver - ht) / rz
                var hpnow = hpver
                let hop = Float(ihop)
                let psi = 0.5 * gcd / hop
                var tdel = (fcos(psi) - rz / (rz + hpnow)) / fsin(psi)
                var cdel = 1 / fsqrt(1 + tdel * tdel)
                var sphe = rz * cdel / (rz + ht)
                var secp = 1 / fsqrt(1 - sphe * sphe)
                var esd = fv * secp
                // Correction to Martyn's theorem, as in CURMUF.
                let fob1 = esd
                let xmut = sphe * sphe
                let xfsq = fob1 * fob1 / (fi[jh, k] * fi[jh, k])
                let sph = xfsq * xmut * xhp * (ht + 2 * (rz + ht) * xhp)
                hpnow = hpver + sph
                tdel = (fcos(psi) - rz / (rz + hpnow)) / fsin(psi)
                cdel = 1 / fsqrt(1 + tdel * tdel)
                sphe = rz * cdel / (rz + ht)
                secp = 1 / fsqrt(1 - sphe * sphe)
                esd = fv * secp
                if freq - esd + eps < 0 { continue }
                ymuf[jh] = esd
                sigl[jh] = esd * osigl[jh] / oymuf[jh]
                yfot[jh] = ymuf[jh] - 1.28 * sigl[jh]
                sigu[jh] = esd * osigu[jh] / oymuf[jh]
                yhpf[jh] = ymuf[jh] + 1.28 * sigu[jh]
                inum = inum + 1
                delmod[inum, k] = fatan(tdel) * r2d
                hpmod[inum, k] = hpnow
                htmod[inum, k] = htmuf[jh]
                fvmod[inum, k] = fvmuf[jh]
                afmod[inum, k] = afmuf[jh]
                itmod[inum, k] = jh
                ghop = 2 * psi
            }
        }
        if zeroDistance && inum < 6 {
            // Very short distance: reproduce the vertical ionogram.
            inum = inum + 1
            let fv = freq - 0.001
            var found = 0
            var ihx = 2
            for ih in 2...30 {
                ihx = ih
                let difference = fv - fvert[ih, k]
                if difference < 0 { found = 1; break }
                if difference == 0 { found = 2; break }
            }
            if found == 1 {
                let slope = fvert[ihx, k] - fvert[ihx - 1, k]
                if slope > 0.01 {
                    let slp = (fv - fvert[ihx - 1, k]) / slope
                    delmod[inum, k] = 90
                    hpmod[inum, k] = hprim[ihx - 1, k] + (hprim[ihx, k] - hprim[ihx - 1, k]) * slp
                    htmod[inum, k] = htrue[ihx - 1, k] + (htrue[ihx, k] - htrue[ihx - 1, k]) * slp
                    fvmod[inum, k] = fv
                    afmod[inum, k] = afac[ihx - 1, k] + (afac[ihx, k] - afac[ihx - 1, k]) * slp
                }
            } else if found == 2 {
                delmod[inum, k] = 90
                hpmod[inum, k] = hprim[ihx, k]
                htmod[inum, k] = htrue[ihx, k]
                fvmod[inum, k] = fvert[ihx, k]
                afmod[inum, k] = afac[ihx, k]
            }
            if found != 0 {
                for ily in 1...3 {
                    if ily == 2 && fi[2, k] <= 0 { continue }
                    if htmod[inum, k] <= hi[ily, k] {
                        itmod[inum, k] = ily
                        inm[ily] = 1
                        break
                    }
                }
            }
        }
        regmod()
        if ireset {
            for il in 1...3 {
                sigl[il] = osigl[il]
                sigu[il] = osigu[il]
                yfot[il] = oyfot[il]
                yhpf[il] = oyhpf[il]
                ymuf[il] = oymuf[il]
            }
        }
    }

    /// REGMOD: transmission loss, signal power and signal-to-noise ratio
    /// for each regular-layer mode at the current hop distance.
    func regmod() {
        let k = jmode
        let l = [1, 3, 5][k - 1]
        let ac = 677.2 * acav
        let bc = powf(freq + gyz[l], 1.98)
        let ihop = fint(gcd / ghop + 0.01)
        let hop = Float(ihop)
        hn[7] = -1
        hp[7] = -1
        for im in 1...6 {
            hn[im] = -1
            hp[im] = -1
            if hpmod[im, k] <= 0 { continue }
            // Free space loss relative to an isotropic radiator.
            var del = d2r * delmod[im, k]
            del = min(del, 89.99 * d2r)
            var cdel = fcos(del)
            var psi = ghop * 0.5
            var phe = pio2 - psi - del
            var path = 2 * (hpmod[im, k] + rz * (1 - fcos(psi))) / fcos(phe)
            path = abs(path * hop)
            timed[im] = path / vofl
            fslos[im] = 32.45 + 20 * alog10(path * freq)
            var obfu: Float = 0, obfl: Float = 0
            if fi[1, k] - fvmod[im, k] >= 0 {
                // D–E mode.
                var xnsq: Float
                if htloss - htmod[im, k] <= 0 {
                    xnsq = 10.2
                } else {
                    let hnux = 61 + 3 * (htmod[im, k] - 70) / 18
                    xnsq = xnuz * fexp(-2 * (hnux - 60) / hnu)
                }
                let heff = min(100, htmod[im, k])
                let sinp = rz * cdel / (rz + heff)
                var secp = 1 / fsqrt(1 - sinp * sinp)
                abps[im] = secp * ac / (bc + xnsq)
                let xv = max(fvmod[im, k] / fi[1, k], xve)
                let adx = afe + bfe * alog(xv)
                secp = 1 / fsin(del + psi)
                adv[im] = secp * afmod[im, k] * (powf(fvmod[im, k] + gyz[l], 1.98) + xnsq) / (bc + xnsq) + adx
                obf[im] = 0
            } else {
                // F layer mode.
                let xnsq: Float = 10.2
                let sinp = rz * cdel / (rz + 100)
                var secp = 1 / fsqrt(1 - sinp * sinp)
                abps[im] = secp * ac / (bc + xnsq)
                secp = 1 / fsin(del + psi)
                adv[im] = secp * afmod[im, k] * (powf(fvmod[im, k] + gyz[l], 1.98) + xnsq) / (bc + xnsq)
                // Sporadic E obscuration.
                obf[im] = 0
                if fs[2, k] > 0 {
                    var fsdead = Float(ifob[1, 3, k])
                    fsdead = fsdead / 1000
                    fsdead = min(fsdead, 3)
                    let fmhz = max(freq, fsdead)
                    let sins = rz * cdel / (rz + hs[k])
                    let secs = 1 / fsqrt(1 - sins * sins)
                    var esd = fs[2, k] * secs
                    var pros = prbmuf(fmhz, esd, ymuf[4], 4)
                    pros = min(pros, 0.90)
                    obf[im] = -10 * alog10(1 - pros)
                    esd = fs[1, k] * secs
                    pros = prbmuf(fmhz, esd, yfot[4], 4)
                    pros = min(pros, 0.90)
                    obfu = -10 * alog10(1 - pros)
                    pros = prbmuf(fmhz, esd, yhpf[4], 4)
                    pros = min(pros, 0.9)
                    obfl = -10 * alog10(1 - pros)
                }
            }
            var y: Float = 0
            for ig in 1...km { y = y + groundLoss(ig, del, freq) }
            grlos[im] = y / Float(km)
            tgain[im] = gain(1, del, freq).rain
            let receive = gain(2, del, freq)
            rgain[im] = receive.rain
            eff[im] = receive.eff
            // Only two hops count for the obscuration.
            let hops = min(hop, 2)
            var xtlos = fslos[im] + hop * (abps[im] + adv[im]) + grlos[im] * (hop - 1)
                + hops * obf[im] + asm - rgain[im] - tgain[im]
            let ismod = itmod[im, k]
            let sphet = rz * cdel / (rz + htmod[im, k])
            var cphet = 1 - sphet * sphet
            cphet = max(0.000001, cphet)
            cphet = fsqrt(cphet)
            // MUFday, from the MUF for this hop rather than for the path.
            psi = gcd / 2
            psi = psi / hop
            let cpsi = fcos(psi)
            let spsi = fsin(psi)
            let tanp = spsi / (1 - cpsi + hpmuf[ismod] / rz)
            phe = fatan(tanp)
            del = pio2 - phe - psi
            cdel = fcos(del)
            let sphe = rz * cdel / (rz + htmuf[ismod])
            var xmuf = fvmuf[ismod] / fsqrt(1 - sphe * sphe)
            prob[im] = prbmuf(freq, xmuf, ymuf[ismod], ismod)
            if prob[im] < 0.0001 {
                // Extra loss, up to 24 dB, as MUFday becomes very small.
                var xfac = prob[im] > 0 ? alog10(prob[im]) : -7
                if xfac > -4 { xfac = -4 }
                if xfac < -7 { xfac = -7 }
                let ghlos = (xfac + 4) * 8
                xtlos = xtlos - ghlos
            }
            xmuf = fvmuf[ismod] / cphet
            var pValue = prbmuf(freq, xmuf, ymuf[ismod], ismod)
            if pValue <= 0.000001 { pValue = 0.000001 }
            let xls = -10 * alog10(pValue) / cphet
            xtlos = xtlos + xls * hop
            let cpr = fvmuf[ismod] / ymuf[ismod]
            let fvfot = yfot[ismod] * cpr
            xmuf = fvfot / cphet
            var pf = prbmuf(freq, xmuf, yfot[ismod], ismod)
            if pf <= 0.000001 { pf = 0.000001 }
            let xlsl = -10 * alog10(pf) / cphet
            let fvhpf = yhpf[ismod] * cpr
            xmuf = fvhpf / cphet
            pf = prbmuf(freq, xmuf, yhpf[ismod], ismod)
            if pf <= 0.000001 { pf = 0.000001 }
            let xlsu = -10 * alog10(pf) / cphet
            // Deciles of signal level.
            tllow[im] = dsl + hops * (obfl - obf[im]) + hop * (xlsl - xls)
            tlhgh[im] = dsu + hops * (obf[im] - obfu) + hop * (xls - xlsu)
            tllow[im] = min(tllow[im], 25)
            tlhgh[im] = min(tlhgh[im], 25)
            tloss[im] = xtlos
            fldst[im] = 107.2 + pwrdb(freq) + 20 * alog10(freq) - xtlos - rgain[im]
            sigpow[im] = pwrdb(freq) - xtlos
            sn[im] = sigpow[im] - rcnse - eff[im]
            b[im] = delmod[im, k]
            nmode[im] = ismod
            hp[im] = hpmod[im, k]
            hn[im] = Float(ihop)
        }
    }
}
