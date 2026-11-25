import Foundation

// Noise at the receiver and the ray set for one frequency: GENOIS, PENANG
// and FINDF.

extension VoacapEngine {
    /// GENOIS: combine atmospheric, galactic and man-made noise into the
    /// median noise power and its deciles (ITS Report 87-212, with the
    /// Caruana modification).
    func genois() {
        let dfac: Float = 7.87384, bfac: Float = 30.99872, cfac = Double(Float(5.56765))
        let dume = min(freq, 55)
        // Atmospheric noise, interpolated between the two local time blocks.
        let first = genfam(y2: rlat, iblk: kj, freq: dume, z: atnu)
        let second = genfam(y2: rlat, iblk: jk, freq: dume, z: atny)
        let slop = abs(cc - tm) / 4
        var atnos = first.fa + (second.fa - first.fa) * slop
        let dua = first.dua + (second.dua - first.dua) * slop
        let dla = first.dla + (second.dla - first.dla) * slop
        let sma = first.dms + (second.dms - first.dms) * slop
        let sua = first.dus + (second.dus - first.dus) * slop
        let sla = first.dls + (second.dls - first.dls) * slop
        let au1 = Double(fexp((dua / dfac) * (dua / dfac) + (atnos / 4.34294)))
        let vu1 = au1 * au1 * Double(fexp(dua * dua / bfac) - 1)
        let al1 = Double(fexp((dla / dfac) * (dla / dfac) + (atnos / 4.34294)))
        let vl1 = al1 * al1 * Double(fexp(dla * dla / bfac) - 1)

        // Galactic noise.
        gnos = 52 - 23 * alog10(freq)
        let dug: Float = 2, dlg: Float = 2, smg: Float = 0.5, sug: Float = 0.2, slg: Float = 0.2
        var galactic = true
        var au2 = Double(fexp((dug / dfac) * (dug / dfac) + (gnos / 4.34294)))
        var vu2 = au2 * au2 * Double(fexp(dug * dug / bfac) - 1)
        var al2 = Double(fexp((dlg / dfac) * (dlg / dfac) + (gnos / 4.34294)))
        var vl2 = al2 * al2 * Double(fexp(dlg * dlg / bfac) - 1)
        if freq <= fi[3, kfx] {
            // Galactic noise does not penetrate the ionosphere.
            galactic = false
            gnos = 0
            au2 = 0; vu2 = 0; al2 = 0; vl2 = 0
        }

        // Man-made noise.
        let xnint: [Float] = [76.8, 72.5, 67.2, 53.6, 83.2, 65.2]
        let conn: [Float] = [27.7, 27.7, 27.7, 28.6, 37.5, 29.1]
        let man = noise
        var xnois = Float(man)
        var ma = abs(man)
        znoise = xnois
        if man > 0 {
            // A 3 MHz level in dB below one watt, converted to Fa at 1 MHz and then to the frequency.
            xnois = 204.0 - xnois + 13.22
            xnois = xnois - 27.7 * alog10(freq)
        } else {
            if man == 0 { ma = 4 }
            ma = min(4, ma)
            xnois = xnint[ma - 1] - conn[ma - 1] * alog10(freq)
            znoise = 204.0 - xnint[ma - 1] + conn[ma - 1] * alog10(3)
        }
        let dum: Float = 9.7, dlm: Float = 6, sum: Float = 1.5, smm: Float = 5.4, slm: Float = 1.5
        let au3 = Double(fexp((dum / dfac) * (dum / dfac) + (xnois / 4.34294)))
        let vu3 = au3 * au3 * Double(fexp(dum * dum / bfac) - 1)
        let al3 = Double(fexp((dlm / dfac) * (dlm / dfac) + (xnois / 4.34294)))
        let vl3 = al3 * al3 * Double(fexp(dlm * dlm / bfac) - 1)

        let vu = vu1 + vu2 + vu3
        let vl = vl1 + vl2 + vl3
        let au = au1 + au2 + au3
        let al = al1 + al2 + al3

        let reff = gain(2, 0, freq).eff
        for im in 1...6 { eff[im] = reff }

        atnos = atnos - 204
        gnos = gnos - 204
        xnois = xnois - 204
        let rnse = 4.34294 * alog(powf(10, atnos * 0.1) + powf(10, gnos * 0.1) + powf(10, xnois * 0.1))

        var sigtsqu = log(1 + vu / (au * au))
        var sigtsql = log(1 + vl / (al * al))
        if !(dua <= 12 && dla <= 12) {
            // Keep Spaulding's method from breaking down when the deciles are wide.
            var dxx = Double(rnse)
            dxx = dxx + 204
            let dxx2 = 2 * (log(au) - dxx / 4.34294)
            let dxx3 = 2 * (log(al) - dxx / 4.34294)
            if sigtsqu > dxx2 && dxx2 > 0 { sigtsqu = dxx2 }
            if sigtsql > dxx3 && dxx3 > 0 { sigtsql = dxx3 }
        }
        let xrnse = Float(Double(Float(4.34294)) * (log(au) - sigtsqu / 2) - 204)
        du = Float(cfac * sigtsqu.squareRoot())
        dl = Float(cfac * sigtsql.squareRoot())

        let qpa = powf(10, (atnos - xrnse) * 0.1)
        var qpg = powf(10, (gnos - xrnse) * 0.1)
        if !galactic { qpg = 0 }
        let qpm = powf(10, (xnois - xrnse) * 0.1)
        sigm = fsqrt((qpa * sma) * (qpa * sma) + (qpg * smg) * (qpg * smg) + (qpm * smm) * (qpm * smm))

        var pv = qpa * fexp((dua - du) * 0.23026)
        let sygu1 = (pv * sua) * (pv * sua) + ((pv - qpa) * sma) * ((pv - qpa) * sma)
        pv = qpg * fexp((dug - du) * 0.23026)
        let sygu2 = (pv * sug) * (pv * sug) + ((pv - qpg) * smg) * ((pv - qpg) * smg)
        pv = qpm * fexp((dum - du) * 0.23026)
        let sygu3 = (pv * sum) * (pv * sum) + ((pv - qpm) * smm) * ((pv - qpm) * smm)
        sxgu = fsqrt(sygu1 + sygu2 + sygu3)

        pv = qpa * fexp((dla - dl) * 0.23026)
        let sygl1 = (pv * sla) * (pv * sla) + ((pv - qpa) * sma) * ((pv - qpa) * sma)
        pv = qpg * fexp((dlg - dl) * 0.23026)
        let sygl2 = (pv * slg) * (pv * slg) + ((pv - qpg) * smg) * ((pv - qpg) * smg)
        pv = qpm * fexp((dlm - dl) * 0.23026)
        let sygl3 = (pv * slm) * (pv * slm) + ((pv - qpm) * smm) * ((pv - qpm) * smm)
        sxgl = fsqrt(sygl1 + sygl2 + sygl3)

        rcnse = xrnse
        xnoise = xnois
    }

    /// PENANG: the elevation angle at which the operating frequency
    /// penetrates each layer.
    func penang(_ k: Int) {
        let fmhz = freq
        var ht = htrue[10, k]
        var fv = fvert[10, k]
        var frat = fv / freq
        frat = frat * frat
        if frat - 0.9999 >= 0 {
            delpen[1, k] = 89.9
            delpen[2, k] = 90.0
            delpen[3, k] = 90.0
            return
        }
        var cdel = (rz + ht) * fsqrt(1 - frat) / rz
        delpen[1, k] = cdel - 0.999999 <= 0 ? facos(cdel) * r2d : 0
        if fi[2, k] <= 0 {
            delpen[2, k] = delpen[1, k]
        } else {
            ht = htrue[20, k]
            fv = fvert[20, k]
            frat = fv / freq
            frat = frat * frat
            if frat - 0.9999 >= 0 {
                delpen[2, k] = 89.9
                delpen[3, k] = 90
                return
            }
            cdel = (rz + ht) * fsqrt(1 - frat) / rz
            delpen[2, k] = cdel - 0.999999 <= 0 ? facos(cdel) * r2d : 0
        }
        if fvert[30, k] - freq + 0.0001 >= 0 {
            delpen[3, k] = 89.9
            return
        }
        let r28 = fvert[28, k] / fmhz, r29 = fvert[29, k] / fmhz, r30 = fvert[30, k] / fmhz
        let xm28 = (rz + htrue[28, k]) * fsqrt(1.0 - r28 * r28)
        let xm29 = (rz + htrue[29, k]) * fsqrt(1.0 - r29 * r29)
        let xm30 = (rz + htrue[30, k]) * fsqrt(1.0 - r30 * r30)
        if xm30 - xm29 >= 0 {
            cdel = xm30 - xm28 < 0 ? xm28 / rz : xm30 / rz
        } else {
            cdel = xm29 - xm28 < 0 ? xm28 / rz : xm29 / rz
        }
        if cdel - 0.999999 <= 0 {
            delpen[3, k] = facos(cdel) * r2d
        } else {
            delpen[3, k] = 0
            dskpkm[k] = -1
            dmaxkm[k] = -1
        }
    }

    /// FINDF: every ray at the operating frequency through control point
    /// `k`, with the layer cusps inserted, plus the skip and maximum distance.
    func findf(_ k: Int) {
        enum Label { case l225, l255, l265, l275, l305, l325, l340, l345, l350, l375, l400 }
        let jfhz = fint(1000 * freq)
        dmaxkm[k] = 0
        dskpkm[k] = 10000
        for i in 1...45 {
            hpflx[i, k] = 0
            delfx[i, k] = 0
            gdflx[i, k] = 0
        }
        let fc2 = fi[3, k] * fi[3, k]
        penang(k)
        var iaIndex = 0
        var iaf = 1
        var icusp = -1
        var il = 1
        var ih = 1
        var ilow = 1
        var ihigh = 10
        var fv: Float = 0, hpValue: Float = 0
        var label = Label.l275
        loop: while true {
            switch label {
            case .l225:
                ih = 11
                ilow = 11
                if fi[2, k] <= 0 {
                    il = 3
                    icusp = -1
                    ihigh = 30
                } else {
                    il = 2
                    icusp = -1
                    ihigh = 20
                }
                label = .l275
            case .l255:
                il = 3
                icusp = -1
                ilow = ihigh + 1
                ihigh = 30
                ih = 21
                label = .l275
            case .l265:
                label = il == 1 ? .l225 : (il == 2 ? .l255 : .l400)
            case .l275:
                if delpen[il, k] <= 0 { label = .l265; continue loop }
                if delpen[il, k] - 89.99 > 0 { label = .l400; continue loop }
                iaIndex = iaIndex + 1
                if iaIndex - nang > 0 { label = .l400; continue loop }
                label = delpen[il, k] - ang[iaIndex] <= 0 ? .l345 : .l305
            case .l305:
                if ifob[iaIndex, ilow, k] - jfhz >= 0 { label = .l325; continue loop }
                if ih - ihigh >= 0 { label = .l275; continue loop }
                let difference = ifob[iaIndex, ih, k] - jfhz
                if difference == 0 { label = .l325; continue loop }
                if difference > 0 {
                    ih = ih + 1
                    continue loop
                }
                if ifob[iaIndex, ih + 1, k] - jfhz < 0 {
                    ih = ih + 1
                    continue loop
                }
                label = .l340
            case .l325:
                delfx[iaf, k] = ang[iaIndex]
                htflx[iaf, k] = htrue[ih, k]
                afflx[iaf, k] = afac[ih, k]
                fv = fvert[ih, k]
                hpValue = hprim[ih, k]
                imode[iaf, k] = il
                label = .l375
            case .l340:
                var slopd = Float(ifob[iaIndex, ih + 1, k] - ifob[iaIndex, ih, k])
                slopd = max(1, slopd)
                var slope = Float(jfhz - ifob[iaIndex, ih, k])
                slope = slope / slopd
                htflx[iaf, k] = htrue[ih, k] + slope * (htrue[ih + 1, k] - htrue[ih, k])
                fv = fvert[ih, k] + slope * (fvert[ih + 1, k] - fvert[ih, k])
                delfx[iaf, k] = ang[iaIndex]
                hpValue = hprim[ih, k] + slope * (hprim[ih + 1, k] - hprim[ih, k])
                afflx[iaf, k] = afac[ih, k] + slope * (afac[ih + 1, k] - afac[ih, k])
                imode[iaf, k] = il
                label = .l375
            case .l345:
                // Insert the cusp at the penetration angle.
                delfx[iaf, k] = delpen[il, k]
                htflx[iaf, k] = htrue[ihigh, k]
                afflx[iaf, k] = afac[ihigh, k]
                fv = fvert[ihigh, k]
                hpValue = hprim[ihigh, k]
                iaIndex = iaIndex - 1
                icusp = 0
                imode[iaf, k] = il
                label = .l375
            case .l350:
                if il - 3 >= 0 { label = .l400; continue loop }
                if delpen[il, k] - 89.9 >= 0 { label = .l400; continue loop }
                // Insert the cusp for the next layer.
                delfx[iaf, k] = delfx[iaf - 1, k] + 0.001
                htflx[iaf, k] = htrue[ihigh + 1, k]
                afflx[iaf, k] = afac[ihigh + 1, k]
                fv = fvert[ihigh + 1, k]
                hpValue = hprim[ihigh + 1, k]
                icusp = 1
                imode[iaf, k] = fi[2, k] <= 0 ? 3 : il + 1
                label = .l375
            case .l375:
                // Correct Martyn's theorem for a spherical ionosphere.
                let del = delfx[iaf, k] * d2r
                let rcosd = rz * fcos(del)
                let xfsq = freq * freq / fc2
                let ht = htflx[iaf, k]
                let xmut = 1 - fv * fv / (freq * freq)
                let xhp = (hpValue - ht) / rz
                let sph = xfsq * xmut * xhp * (ht + 2 * (rz + ht) * xhp)
                hpValue = hpValue + sph
                var phe = rcosd / (rz + hpValue)
                phe = fasin(phe)
                let gdr = 2 * rz * (pio2 - del - phe)
                gdflx[iaf, k] = gdr
                hpflx[iaf, k] = hpValue
                fvflx[iaf, k] = fv
                if dskpkm[k] - gdr > 0 {
                    dskpkm[k] = gdr
                    delskp[k] = delfx[iaf, k]
                    htskp[k] = ht
                    hpskp[k] = hpValue
                    fvskp[k] = fv
                    iskp[k] = il
                }
                if dmaxkm[k] - gdr <= 0 && delfx[iaf, k] - amind >= 0 {
                    dmaxkm[k] = gdr
                }
                iaftxr[k] = iaf
                iaf = iaf + 1
                if iaf - 45 > 0 { label = .l400; continue loop }
                label = icusp < 0 ? .l275 : (icusp == 0 ? .l350 : .l265)
            case .l400:
                break loop
            }
        }
    }
}
