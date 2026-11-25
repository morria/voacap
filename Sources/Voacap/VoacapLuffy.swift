import Foundation

// LUFFY: system performance for every frequency at one hour, with the
// short-path model, the long-path model, and VOACAP method 30's smoothing
// between them from 7000 to 10000 km.

extension VoacapEngine {
    /// Electron density profile, ionogram, reflectrix and deviative loss for control point `k`.
    private func prepare(_ k: Int) {
        lecden(k)
        genion(k)
        fobby(k)
        alosfv(k)
    }

    /// The short-path model for the current frequency.
    private func shortPath(k: Int, ihmin: Int) {
        findf(k)
        var ihsrt: Int, ihstp: Int
        if dmaxkm[k] <= 0 {
            // Only one over-the-MUF mode.
            ihsrt = ihmin
            ihstp = ihmin
        } else {
            // Up to three hop counts.
            ihsrt = fint(gcdkm / dmaxkm[k] + 1)
            var ihmax = fint(gcdkm / dskpkm[k])
            if ihsrt < ihmin { ihsrt = ihmin }
            if ihmax < ihsrt { ihmax = ihsrt }
            ihstp = min(ihmax, ihsrt + 2)
            if ihsrt > ihmin { ihsrt = max(ihmin, ihstp - 2) }
        }
        if ihsrt <= ihstp {
            for ihop in ihsrt...ihstp {
                ghop = gcd / Float(ihop)
                fdist(k)
                inmuf(ihop)
                allmodes(100)
            }
        }
        esmod()
        esreg()
        allmodes(2)
    }

    /// The long-path model for the current frequency.
    private func longPath() {
        findf(itxrcp[1])
        findf(itxrcp[2])
        gmloss()
        seltxr()
        lngpat()
    }

    /// LUFFY, as VOACAP runs it for method 30.
    func luffy() {
        let mspecSaved = mspec
        defer { mspec = mspecSaved }
        var ipflag = gcdkm >= gcdlng ? 200 : 100
        if gcdkm >= gcdlng { mspec = 0 }
        selmod()
        anois1()

        // Results of the short (1) and long (2) passes, kept for smoothing.
        var yangle = A2(13, 2), ycprob = A2(13, 2), ydblos = A2(13, 2), ydblosl = A2(13, 2)
        var ydblosu = A2(13, 2), ydelay = A2(13, 2), ydbw = A2(13, 2), yvhigh = A2(13, 2)
        var xdbu = A2(13, 2), xreliab = A2(13, 2), xsndb = A2(13, 2), xsnpr = A2(13, 2)
        var xsnrlw = A2(13, 2), xsnrup = A2(13, 2), ytgain = A2(13, 2), yrgain = A2(13, 2)
        var inhp = I2(13, 2), ymode = I2(13, 2)
        var zdu = A1(13), zdl = A1(13), zangler = A1(13)
        var inmode = I1(13), zmoder = I1(13)

        while true {
            ipfg = ipflag
            var k: Int
            if ipfg == 100 {
                jlong = -1
                k = jmode
                prepare(k)
            } else {
                jlong = 1
                k = 1
                itxrcp[1] = 1
                if !(mspec == 121 && jmode == k) { prepare(k) }
                k = kfx
                itxrcp[2] = kfx
                if !(mspec == 121 && jmode == k) {
                    if kfx == 1 {
                        itxrcp[2] = 2
                    } else {
                        prepare(k)
                    }
                }
            }
            setlng()
            sigdis()
            var ihmin = nhopmf[1]
            if nhopmf[2] > 0 && nhopmf[2] < ihmin { ihmin = nhopmf[2] }
            if nhopmf[3] < ihmin { ihmin = nhopmf[3] }
            let model: Character = ipfg == 100 ? "S" : "L"
            for ifx in 1...12 {
                mdl[ifx] = " "
                nhp[ifx] = -1
                sprob[ifx] = 0
                probmp[ifx] = 0
                if frel[ifx] <= 0 { continue }
                mdl[ifx] = model
                freq = frel[ifx]
                allmodes(0)
                genois()
                if ipfg == 100 {
                    shortPath(k: k, ihmin: ihmin)
                } else {
                    longPath()
                }
                genois()
                relbil(ifx)
                if ipfg != 100 {
                    moder[ifx] = iznmode[2]
                    angler[ifx] = zb[2]
                }
                if mspec == 121 {
                    let idx = (ipflag + 1) / 100
                    yangle[ifx, idx] = angle[ifx]
                    ycprob[ifx, idx] = cprob[ifx]
                    ydblos[ifx, idx] = dblos[ifx]
                    ydblosl[ifx, idx] = dblosl[ifx]
                    ydblosu[ifx, idx] = dblosu[ifx]
                    ydelay[ifx, idx] = delay[ifx]
                    ydbw[ifx, idx] = dbw[ifx]
                    inhp[ifx, idx] = nhp[ifx]
                    yvhigh[ifx, idx] = vhigh[ifx]
                    ymode[ifx, idx] = mode[ifx]
                    xdbu[ifx, idx] = dbu[ifx]
                    xreliab[ifx, idx] = reliab[ifx]
                    xsndb[ifx, idx] = sndb[ifx]
                    xsnpr[ifx, idx] = snpr[ifx]
                    xsnrlw[ifx, idx] = snrlw[ifx]
                    xsnrup[ifx, idx] = snrup[ifx]
                    ytgain[ifx, idx] = ztgain[1]
                    yrgain[ifx, idx] = zrgain[1]
                    if ipflag == 200 {
                        zdu[ifx] = du
                        zdl[ifx] = dl
                        inmode[ifx] = iznmode[1]
                        zangler[ifx] = angler[ifx]
                        zmoder[ifx] = moder[ifx]
                    }
                }
            }
            // Both models are needed only from 7000 to 10000 km.
            if mspec != 121 || gcdkm < 7000 || gcdkm > 10000 { return }
            if ipflag == 100 {
                ipflag = 200
                continue
            }
            break
        }

        // Long path / short path smoothing, from the VOA memo of 15 January 1991.
        for i in 1...12 {
            mdl[i] = " "
            freq = frel[i]
            if freq <= 0 { continue }
            let slpld = ydbw[i, 2] - abs(ydblosl[i, 2])
            let sspld = ydbw[i, 1] - abs(ydblosl[i, 1])
            let delpow = slpld - sspld
            var idx: Int
            if delpow < 0 {
                // The short path is stronger.
                idx = 1
                mdl[i] = "S"
                moder[i] = ymode[i, 1]
                angler[i] = yangle[i, 1]
                gaint[i] = ytgain[i, 1]
                gainr[i] = yrgain[i, 1]
                snxx[i] = rsn - xsnpr[i, 1]
            } else {
                idx = 2
                mdl[i] = "L"
                moder[i] = zmoder[i]
                angler[i] = zangler[i]
                du = zdu[i]
                dl = zdl[i]
                zb[1] = yangle[i, 2]
                ztllow[1] = ydblosl[i, 2]
                ztlhgh[1] = ydblosu[i, 2]
                ztimed[1] = ydelay[i, 2]
                zhn[1] = Float(inhp[i, 2])
                zhp[1] = yvhigh[i, 2]
                zprob[1] = ycprob[i, 2]
                iznmode[1] = inmode[i]
                ztgain[1] = ytgain[i, 2]
                zrgain[1] = yrgain[i, 2]
                zeff[1] = rneff[i]
            }
            mode[i] = ymode[i, idx]
            angle[i] = yangle[i, idx]
            cprob[i] = ycprob[i, idx]
            dblos[i] = ydblos[i, idx]
            dblosl[i] = ydblosl[i, idx]
            dblosu[i] = ydblosu[i, idx]
            delay[i] = ydelay[i, idx]
            dbw[i] = ydbw[i, idx]
            nhp[i] = inhp[i, idx]
            rcnse = xnynois[i]
            vhigh[i] = yvhigh[i, idx]
            dbu[i] = xdbu[i, idx]
            reliab[i] = xreliab[i, idx]
            sndb[i] = xsndb[i, idx]
            snpr[i] = xsnpr[i, idx]
            snrlw[i] = xsnrlw[i, idx]
            snrup[i] = xsnrup[i, idx]
            if delpow >= 0 && gcdkm < 10000 {
                mdl[i] = "M"
                let disint = (gcdkm - 7000) / 3000
                let delx = delpow > 375 ? powf(10, 37.5) : powf(10, delpow / 10)
                let smooth = sspld + 10 * alog10(disint * (delx - 1) + 1)
                zsigpow[1] = smooth + ztllow[1]
                ztloss[1] = pwrdb(freq) - zsigpow[1]
                zsn[1] = zsigpow[1] - rcnse - zeff[1]
                zfldst[1] = 107.2 + pwrdb(freq) + 20 * alog10(freq) - ztloss[1] - zrgain[1]
                nmmod = 1
                genois()
                relbil(i)
            }
        }
    }
}

// MARK: - Running a prediction

/// What VOACAP predicts for one frequency at one hour.
public struct VoacapResult: Sendable {
    public var frequency: Float
    /// Fraction of days in the month the required signal-to-noise ratio is met.
    public var reliability: Float
    /// Median signal-to-noise ratio in dB for a 1 Hz bandwidth.
    public var snr: Float
    /// Median signal power at the receiver, dBW.
    public var signal: Float
    /// Median noise power in a 1 Hz bandwidth, dBW.
    public var noise: Float
    /// Fraction of days the frequency is below the MUF of the mode used.
    public var mufDay: Float
    /// Transmission loss, dB.
    public var loss: Float
    /// Take-off angle of the most reliable mode, degrees.
    public var angle: Float
    public var hops: Int
    /// 1 E, 2 F1, 3 F2, 4 Es.
    public var layer: Int
    /// "S" short-path model, "L" long-path model, "M" smoothed.
    public var model: Character
    /// Signal-to-noise ratio exceeded on the required fraction of days.
    public var snrAtRequiredReliability: Float
}

extension VoacapEngine {
    /// Set the circuit. Longitudes are positive east.
    public func setCircuit(txLat: Float, txLon: Float, rxLat: Float, rxLon: Float, month: Int, ssn: Float) {
        tlatd = txLat
        tlongd = txLon
        rlatd = rxLat
        rlongd = rxLon
        redmap(ssn: ssn, month: month)
        geom()
        sang()
    }

    /// Run VOACAP method 30 for one UT hour (1–24) at up to eleven frequencies.
    /// Returns the circuit MUF and one result per frequency.
    public func predict(hour: Int, frequencies: [Float]) -> (muf: Float, results: [VoacapResult]) {
        geotim(hour)
        virtim()
        timvar()
        f2var()
        esind()
        ionset()
        curmuf()
        frel.fill(0)
        let count = min(frequencies.count, 11)
        for i in 0..<count { frel[i + 1] = frequencies[i] }
        frel[12] = allmuf[it]
        luffy()
        let results = (0..<count).map { index -> VoacapResult in
            let i = index + 1
            return VoacapResult(frequency: frequencies[index], reliability: reliab[i], snr: sndb[i],
                                signal: dbw[i], noise: xnynois[i], mufDay: cprob[i], loss: dblos[i],
                                angle: angle[i], hops: nhp[i], layer: mode[i], model: mdl[i],
                                snrAtRequiredReliability: snxx[i])
        }
        return (allmuf[it], results)
    }
}
