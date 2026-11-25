import Foundation

/// A Swift translation of the VOACAP prediction engine (NTIA/ITS, public
/// domain; this follows the voacapl 16.1207W sources routine by routine).
///
/// Names and array indices are kept from the Fortran so the two can be read
/// side by side. Each COMMON block is a group of stored properties. One
/// engine holds all the state for one circuit, so separate engines can run on
/// separate threads.
public final class VoacapEngine {
    // MARK: Constants (COMMON /CON/)
    let d2r: Float = 0.01745329251
    let dcl: Float = 1.28
    let gama: Float = 0.57721566
    let pi: Float = 3.1415926
    let pi2: Float = 6.283185307
    let pio2: Float = 1.570796326
    let r2d: Float = 57.295779513
    let rz: Float = 6370
    let vofl: Float = 299.79246

    // MARK: Coefficients (COMMON /ONE/, /TWO/, /TWO_AB/, /FONE/)
    let coefficients: VoacapCoefficients
    var ikim = I2(10, 6)
    /// ESLCOF, ESMCOF, ESUCOF, F2COF, FM3COF, ERCOF laid end to end, as the
    /// Fortran EQUIVALENCE has them.
    var cofion = A1(2604)
    let ia = I1([1, 276, 703, 978, 1966, 2407])
    let ib = I1([5, 7, 5, 13, 9, 9])
    let icab = I1([1, 56, 117, 172, 248, 297])
    var ab = A1(318)
    var f2d = A3(16, 6, 6)
    var p = A3(29, 16, 8)
    var abp = A2(2, 9)
    var dud = A3(5, 12, 5)
    var fam = A2(14, 12)
    var sys = A3(9, 16, 6)
    var perr = A3(9, 4, 6)
    var anew = A1(3), bnew = A1(3), achi = A1(2), bchi = A1(2)
    var gamma = A1(6)

    // MARK: Path geometry (COMMON /DON/)
    var alatd: Float = 0, amin: Float = 0.05236, amind: Float = 3
    var btr: Float = 0, btrd: Float = 0, dlong: Float = 0, dmp: Float = 0.85
    var ertr: Float = 1, gcd: Float = 0, gcdkm: Float = 0, pmp: Float = 10
    var pwr: Float = 1, tlat: Float = 0, tlatd: Float = 0, tlong: Float = 0, tlongd: Float = 0
    var rsn: Float = 1, sigtr: Float = 1, rlat: Float = 0, rlatd: Float = 0
    var rlong: Float = 0, rlongd: Float = 0, brtd: Float = 0, flux: Float = 1
    var ssn: Float = 0
    var d90r: Float = 0, d50r: Float = 0, d10r: Float = 0
    var d90s: Float = 0, d50s: Float = 0, d10s: Float = 0
    /// Geomagnetic north pole. voacapl falls back to these 1960-epoch values, so they are used here too.
    var gMagneticLat: Float = 78.5, gMagneticLon: Float = -69.0

    // MARK: Sample areas (COMMON /RON/, /GEOG/, /MFAC/, /ES/, /PSCA/)
    var clat = A1(5), clong = A1(5), glat = A1(5), rd = A1(5)
    var fi = A2(3, 5), yi = A2(3, 5), hi = A2(3, 5)
    var hprim = A2(30, 5), htrue = A2(30, 5), fvert = A2(30, 5), afac = A2(30, 5)
    var km = 1, kfx = 1
    var htr = A2(50, 3), fnsq = A2(50, 3)
    var gyz = A1(5), rat = A1(5), gmdip = A1(5), clck = A1(5), abiy = A1(5), artic = A1(5)
    var sigpat = A1(5, 0.001), epspat = A1(5, 4)
    var f2m3 = A1(5), hpf2 = A1(5), zenang = A1(5), zenmax = A1(5)
    var iedp = -1
    var fsecv = A1(3)
    var fs = A2(3, 5), hs = A1(5)
    var psc = A1([1, 1, 1, 0.7])

    // MARK: Time (COMMON /TIME/, /SSP/, /ION/)
    var it = 1, gmt: Float = 0, utime = A1(24), gmtr: Float = 0, xlmt = A1(24), jtx = 1
    var month = 0
    /// Latitude of the subsolar point at the start and end of each month.
    let sun: [Float] = [-23.05, -17.31, -17.30, -7.89, -7.88, 4.21, 4.26, 14.83, 14.84, 21.93, 21.93, 23.45,
                        23.15, 18.23, 18.20, 8.68, 8.55, -2.86, -2.90, -14.16, -14.20, -21.68, -21.66, -23.45]
    var lufp = 90
    /// Man-made noise at 3 MHz, as a positive number of dB below 1 W.
    var noise = 145
    /// 1 selects the long way round the great circle.
    var npsl = 0

    // MARK: Elevation scan (COMMON /RAYS/)
    let ang = A1([0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,
                  14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44,
                  46, 48, 50, 52, 54, 56, 60, 65, 70, 75, 80, 85, 89.99])
    /// Oblique frequency in kHz for each elevation angle, ionogram point and control point.
    var ifob = I3(40, 30, 5)
    var nang = 0

    // MARK: Gaussian integration (COMMON /DATR/)
    let xt = A1([0.0387724175, 0.1160840707, 0.1926975807, 0.2681521850,
                 0.3419940908, 0.4137792043, 0.4830758017, 0.5494671251,
                 0.6125538897, 0.6719566846, 0.7273182552, 0.7783056514,
                 0.8246122308, 0.8659595032, 0.9020988070, 0.9328128083,
                 0.9579168192, 0.9772599500, 0.9907262387, 0.9982377097])
    let wt = A1([0.0775059480, 0.0770398182, 0.0761103619, 0.0747231691,
                 0.0728865824, 0.0706116474, 0.0679120458, 0.0648040135,
                 0.0613062425, 0.0574397691, 0.0532278470, 0.0486958076,
                 0.0438709082, 0.0387821680, 0.0334601953, 0.0279370070,
                 0.0222458492, 0.0164210584, 0.0104982845, 0.0045212771])
    let nt = 20
    let npl = 1
    let xnpl: Float = 1
    let twdiv: Float = 0.5

    // MARK: Circuit MUFs (COMMON /MUFS/)
    var emuf = A1(24), f1muf = A1(24), f2muf = A1(24), esmuf = A1(24), allmuf = A1(24)
    var fot = A1(24), xluf = A1(24), hpf = A1(24), angmuf = A1(24)
    var modmuf = 0
    var sigl = A1(4), sigu = A1(4), delmuf = A1(4), hpmuf = A1(4), htmuf = A1(4)
    var fvmuf = A1(4), afmuf = A1(4), yfot = A1(4), yhpf = A1(4), ymuf = A1(4)
    var nhopmf = I1(4)

    // MARK: Frequencies (COMMON /FRQ/)
    var frel = A1(29)
    var freq: Float = 0
    var jmode = 0
    var itxrcp = I1(2)

    // MARK: Noise (COMMON /ANOIS/)
    var atnu: Float = -204, atny: Float = -204, cc: Float = 1, tm: Float = 1
    var rcnse: Float = -204, du: Float = 9, dl: Float = 4, sigm: Float = 1.5
    var sxgu: Float = 1, sxgl: Float = 1
    var kj = 1, jk = 1

    // MARK: System loss and reliability terms (COMMON /TON/)
    var adj: Float = 1, ads: Float = 1, gnos: Float = 1, got: Float = 1, rel: Float = 0.01
    var sl: Float = 1, sls: Float = 1, spr: Float = 0.01, su: Float = 1, sus: Float = 1
    var xnoise: Float = 1, znoise: Float = -204

    // MARK: Signal distribution (COMMON /SIGD/)
    var dsl: Float = 0, asm: Float = 0, dsu: Float = 0, aglat: Float = 0
    var dslf: Float = 0, asmf: Float = 0, dsuf: Float = 0, acav: Float = 0, feav: Float = 0
    var afe: Float = 0, bfe: Float = 0, hnu: Float = 0, htloss: Float = 0, xnuz: Float = 0, xve: Float = 0

    // MARK: Candidate modes for one hop count (COMMON /ZON/)
    var abps = A1(7), crel = A1(7), eff = A1(7), fldst = A1(7), grlos = A1(7)
    var hn = A1(7), hp = A1(7), prob = A1(7), rely = A1(7), rgain = A1(7), sigpow = A1(7), sn = A1(7)
    var spro = A1(7), tgain = A1(7), timed = A1(7), tloss = A1(7), b = A1(7), fslos = A1(7), adv = A1(7)
    var obf = A1(7), tllow = A1(7), tlhgh = A1(7)
    var nmode = I1(7)

    // MARK: Every mode found for one frequency (COMMON /allMODE/)
    var zabps = A1(20), zcrel = A1(20), zfldst = A1(20), zhn = A1(20), zhp = A1(20)
    var zprob = A1(20), zrely = A1(20), zrgain = A1(20), zsigpow = A1(20), zsn = A1(20)
    var zspro = A1(20), ztgain = A1(20), ztimed = A1(20), ztloss = A1(20), zb = A1(20), zfslos = A1(20)
    var zgrlos = A1(20), zadv = A1(20), zobf = A1(20)
    var iznmode = I1(22)
    var ztllow = A1(20), ztlhgh = A1(20), zeff = A1(20)
    var nrel = 1, nmmod = 1

    // MARK: Results by frequency (COMMON /SON/)
    var angle = A1(13), angler = A1(13), cprob = A1(13), dblos = A1(13)
    var dblosl = A1(13), dblosu = A1(13), dbu = A1(13), delay = A1(13), dbw = A1(13)
    var nhp = I1(13)
    var xnynois = A1(13), probmp = A1(13), reliab = A1(13), sndb = A1(13), snpr = A1(13)
    var snrlw = A1(13), snrup = A1(13), sprob = A1(13), vhigh = A1(13), rneff = A1(13)
    /// "S" short path, "L" long path, "M" smoothed between the two.
    var mdl = [Character](repeating: " ", count: 14)
    /// Layer of the most reliable mode at each end (1 E, 2 F1, 3 F2, 4 Es, 5 none).
    var mode = I1(13), moder = I1(13)
    var gaint = A1(13), gainr = A1(13), snxx = A1(13)

    // MARK: Reflectrix (COMMON /REFLX/)
    var delfx = A2(45, 3), hpflx = A2(45, 3), htflx = A2(45, 3), gdflx = A2(45, 3), fvflx = A2(45, 3)
    var dskpkm = A1(3), delskp = A1(3), hpskp = A1(3), htskp = A1(3), dmaxkm = A1(3), fvskp = A1(3)
    var iskp = I1(3)
    var imode = I2(45, 3)
    var afflx = A2(45, 3), delpen = A2(3, 5), gml = A2(45, 3), fhp = A2(45, 3)

    // MARK: Modes for one hop (COMMON /MODES/)
    var ghop: Float = 0
    var delmod = A2(6, 3), hpmod = A2(6, 3), htmod = A2(6, 3), fvmod = A2(6, 3), afmod = A2(6, 3)
    var itmod = I2(6, 3)

    // MARK: Long path (COMMON /LPATH/, /LOSX/)
    let gcdlng: Float = 10000
    var txrgml = A2(45, 2)
    let delopt: Float = 3, gmin: Float = 3, ymin: Float = 0.1
    var ltxrgm = I1([1, 1])
    var andvx = A2(45, 3), advx = A2(45, 3), aofx = A2(45, 3), arfx = A2(45, 3)
    var grlosx = A2(45, 3), tgainx = A2(45, 3), tlskm = A2(45, 3)
    var efflp = A1(45)
    var iaftxr = I1(3)

    // MARK: Control
    /// 121 selects smoothing between the short and long path models (method 30).
    var mspec = 121
    var jlong = -1
    var ipfg = 100

    // MARK: Antennas and power (COMMON /CANTENNA/, /PANTENNA/)
    var txAntenna = VoacapAntenna.isotropic(gain: 0)
    var rxAntenna = VoacapAntenna.isotropic(gain: 0)
    /// Transmitter power in dBW.
    var pwrdba: Float = 30

    /// Deviative loss coefficient left over from the last E-layer step of ALOSFV.
    var alosfvCF: Float = 0

    // STATE-INSERT

    public init(coefficients: VoacapCoefficients) {
        self.coefficients = coefficients
        fi[3, 1] = 0.2; fi[3, 2] = 0.2; fi[3, 3] = 0.2; fi[3, 4] = 0.2; fi[3, 5] = 0.2
        for k in 1...5 {
            hi[1, k] = 110; hi[3, k] = 300
            yi[1, k] = 20; yi[3, k] = 100
        }
    }

    @inline(__always) func sunLat(_ i: Int, _ month: Int) -> Float { sun[(month - 1) * 2 + i - 1] }
}

/// The monthly ionospheric and noise coefficient tables, read from the
/// bundled `voacap-coeffs.bin` (built from VOACAP's `coeffNN.asc` and
/// `fof2CCIR.asc`).
public final class VoacapCoefficients: @unchecked Sendable {
    struct Month {
        var ikim: [Int32]
        var dud, fam, sys, fakp, fakabp, xfm3cf, fone, fakmap, abmap, f2d, perr: [Float]
        var xesmcf, xpmap, xeslcf, xesucf, xercof, xf2cof: [Float]
    }
    let months: [Month]

    public init?(data: Data) {
        let perMonth = 11596 * 4
        guard data.count == perMonth * 12 else { return nil }
        var months: [Month] = []
        data.withUnsafeBytes { raw in
            for m in 0..<12 {
                var offset = m * perMonth
                func ints(_ n: Int) -> [Int32] {
                    defer { offset += n * 4 }
                    return (0..<n).map { Int32(littleEndian: raw.loadUnaligned(fromByteOffset: offset + $0 * 4, as: Int32.self)) }
                }
                func floats(_ n: Int) -> [Float] {
                    defer { offset += n * 4 }
                    return (0..<n).map { Float(bitPattern: UInt32(littleEndian: raw.loadUnaligned(fromByteOffset: offset + $0 * 4, as: UInt32.self))) }
                }
                let ikim = ints(60)
                let dud = floats(300), fam = floats(168), sys = floats(864), fakp = floats(2784)
                let fakabp = floats(12), xfm3cf = floats(882), fone = floats(10), fakmap = floats(464)
                let abmap = floats(6), f2d = floats(576), perr = floats(216), xesmcf = floats(854)
                let xpmap = floats(928), xeslcf = floats(550), xesucf = floats(550), xercof = floats(396)
                let xf2cof = floats(1976)
                months.append(Month(ikim: ikim, dud: dud, fam: fam, sys: sys, fakp: fakp, fakabp: fakabp,
                                    xfm3cf: xfm3cf, fone: fone, fakmap: fakmap, abmap: abmap, f2d: f2d,
                                    perr: perr, xesmcf: xesmcf, xpmap: xpmap, xeslcf: xeslcf,
                                    xesucf: xesucf, xercof: xercof, xf2cof: xf2cof))
            }
        }
        self.months = months
    }

    /// The CCIR coefficients that ship with the package.
    public static let bundled: VoacapCoefficients? = {
        guard let url = Bundle.module.url(forResource: "voacap-coeffs", withExtension: "bin"),
              let data = try? Data(contentsOf: url) else { return nil }
        return VoacapCoefficients(data: data)
    }()
}
