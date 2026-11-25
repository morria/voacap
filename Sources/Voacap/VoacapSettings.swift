import Foundation

/// The inputs VOACAP takes besides the circuit itself: the SYSTEM, FPROB and
/// ANTENNA cards of a VOACAP input deck.
public struct VoacapSettings: Sendable {
    /// Transmitter power in watts.
    public var powerWatts = 100.0
    /// Man-made noise at 3 MHz in dB below one watt in 1 Hz: 140 industrial,
    /// 145 residential, 150 rural, 164 remote.
    public var manMadeNoise = 145
    /// Lowest take-off angle considered, in degrees.
    public var minimumAngle = 3.0
    /// Required signal-to-noise ratio in dB for a 1 Hz bandwidth.
    public var requiredSNR = 31.0
    /// Required circuit reliability in percent, used for the
    /// "SNR at required reliability" result.
    public var requiredReliability = 90
    /// Whether sporadic E takes part (FPROB 0.7 for Es rather than 0).
    public var sporadicE = false
    /// Go the long way round the great circle.
    public var longWay = false
    public var transmitAntenna = VoacapAntenna.isotropic(gain: 0)
    public var receiveAntenna = VoacapAntenna.isotropic(gain: 0)
    /// Geomagnetic north pole. voacapl runs with 78.5°N 69°W; the Windows
    /// VOACAP distribution ships 79.5°N 69°W.
    public var magneticPoleLatitude = 78.5
    public var magneticPoleLongitude = -69.0

    public init() {}

    /// The defaults for area maps: 100 W into half-wave dipoles half a
    /// wavelength high and broadside to every path, residential noise.
    public static let area: VoacapSettings = {
        var settings = VoacapSettings()
        let dipole = VoacapAntenna.horizontalDipole(offAzimuth: 0, height: -0.5)
        settings.transmitAntenna = dipole
        settings.receiveAntenna = dipole
        return settings
    }()
}

extension VoacapEngine {
    /// Apply settings. Call before `setCircuit`.
    public func apply(_ settings: VoacapSettings) {
        psc[4] = settings.sporadicE ? 0.7 : 0
        noise = settings.manMadeNoise
        amind = Float(settings.minimumAngle)
        if amind <= 0 { amind = 3 }
        amin = amind * d2r
        rsn = Float(settings.requiredSNR)
        lufp = settings.requiredReliability
        npsl = settings.longWay ? 1 : 0
        pwrdba = 30 + 10 * log10f(Float(max(settings.powerWatts, 1e-6) / 1000))
        txAntenna = settings.transmitAntenna
        rxAntenna = settings.receiveAntenna
        gMagneticLat = Float(settings.magneticPoleLatitude)
        gMagneticLon = Float(settings.magneticPoleLongitude)
    }

    /// Great-circle distance of the current circuit in kilometres.
    public var distanceKm: Double { Double(gcdkm) }
    /// Bearing from transmitter to receiver in degrees.
    public var bearing: Double { Double(btrd) }
}

extension VoacapAntenna {
    /// Read a gain table in the format VOACAP writes to `gainNN.dat`.
    public init?(gainFile text: String) {
        let lines = text.components(separatedBy: .newlines).filter { !$0.isEmpty }.map { Array($0) }
        var antenna = VoacapAntenna.isotropic(gain: 0)
        var n = 2
        for mhz in 1...30 {
            guard n < lines.count, lines[n].count >= 8,
                  let efficiency = Float(String(lines[n][2..<8]).trimmingCharacters(in: .whitespaces)) else { return nil }
            antenna.efficiency[mhz - 1] = efficiency
            var values: [Float] = []
            while values.count < 91 {
                guard n < lines.count else { return nil }
                let line = lines[n]
                n += 1
                var i = 9
                while i + 7 <= line.count && values.count < 91 {
                    guard let value = Float(String(line[i..<i + 7]).trimmingCharacters(in: .whitespaces)) else { return nil }
                    values.append(value)
                    i += 7
                }
                if i == 9 { return nil }
            }
            for elevation in 0...90 { antenna.gain[elevation * 30 + mhz - 1] = values[elevation] }
        }
        self = antenna
    }
}
