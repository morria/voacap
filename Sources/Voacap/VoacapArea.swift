import Foundation

/// VOACAP run from one station to every cell of a polar (bearing by distance)
/// grid, for one month, sunspot number and UT hour.
public struct VoacapArea: Sendable {
    public struct Request: Hashable, Sendable {
        public var latitude: Double
        public var longitude: Double
        /// 1–12
        public var month: Int
        public var ssn: Int
        /// UT hour, 0–23.
        public var hour: Int

        public init(latitude: Double, longitude: Double, month: Int, ssn: Int, hour: Int) {
            self.latitude = latitude
            self.longitude = longitude
            self.month = month
            self.ssn = ssn
            self.hour = hour
        }
    }

    /// The transmitter power the stored signal-to-noise ratios are for.
    public static let referenceWatts = 100.0

    /// Median signal-to-noise ratio in a 1 Hz bandwidth, indexed `[frequency][bearing * rings + ring]`.
    public let snr: [[Float]]
    /// Spread from the median down to the value exceeded on 90% of days, dB.
    public let low: [[Float]]
    /// Spread from the median up to the value exceeded on 10% of days, dB.
    public let high: [[Float]]

    /// Fraction of days the required signal-to-noise ratio (dB in 1 Hz) is met,
    /// as VOACAP's RELBIL computes it.
    public static func reliability(snr: Float, low: Float, high: Float, required: Float) -> Float {
        var z = required - snr
        z = z <= 0 ? z / (low / 1.28) : z / (high / 1.28)
        var yp = abs(z)
        yp = min(5.0, yp)
        var qx = 1 + yp * (0.196854 + yp * (0.115194 + yp * (0.000344 + yp * 0.019527)))
        qx = qx * qx * qx * qx
        qx = 0.5 / qx
        let cumulative = z < 0 ? qx : 1 - qx
        return 1 - cumulative
    }

    /// The point `km` along a great circle from a start point on an initial bearing.
    public static func destination(latitude: Double, longitude: Double, bearing: Double, km: Double) -> (latitude: Double, longitude: Double) {
        let d = km / 6371.0
        let lat1 = latitude * .pi / 180, lon1 = longitude * .pi / 180, brg = bearing * .pi / 180
        let lat2 = asin(min(max(sin(lat1) * cos(d) + cos(lat1) * sin(d) * cos(brg), -1), 1))
        let lon2 = lon1 + atan2(sin(brg) * sin(d) * cos(lat1), cos(d) - sin(lat1) * sin(lat2))
        var degrees = lon2 * 180 / .pi
        degrees = (degrees + 540).truncatingRemainder(dividingBy: 360) - 180
        return (lat2 * 180 / .pi, degrees)
    }

    /// One bearing's worth of circuits. Runs on whichever thread calls it,
    /// with its own engine.
    private static func spoke(_ request: Request, bearing: Double, frequencies: [Float], rings: Int, ringKm: Double,
                              settings: VoacapSettings, coefficients: VoacapCoefficients) -> [[Float]] {
        let engine = VoacapEngine(coefficients: coefficients)
        engine.apply(settings)
        var result = [[Float]](repeating: [Float](repeating: 0, count: rings * 3), count: frequencies.count)
        for ring in 0..<rings {
            let point = destination(latitude: request.latitude, longitude: request.longitude,
                                    bearing: bearing, km: (Double(ring) + 0.5) * ringKm)
            engine.setCircuit(txLat: Float(request.latitude), txLon: Float(request.longitude),
                              rxLat: Float(point.latitude), rxLon: Float(point.longitude),
                              month: request.month, ssn: Float(request.ssn))
            let prediction = engine.predict(hour: request.hour == 0 ? 24 : request.hour, frequencies: frequencies)
            for (index, value) in prediction.results.enumerated() {
                let usable = value.hops > 0 && value.snr.isFinite && engine.snrlw[index + 1].isFinite && engine.snrup[index + 1].isFinite
                result[index][ring * 3] = usable ? value.snr : -300
                result[index][ring * 3 + 1] = usable ? max(engine.snrlw[index + 1], 0.2) : 10
                result[index][ring * 3 + 2] = usable ? max(engine.snrup[index + 1], 0.2) : 10
            }
        }
        return result
    }

    /// Run every circuit, one bearing per task.
    public static func compute(_ request: Request, frequencies: [Float], bearings: Int, rings: Int, ringKm: Double,
                        settings: VoacapSettings = .area, coefficients: VoacapCoefficients) async -> VoacapArea? {
        let cells = bearings * rings
        var snr = [[Float]](repeating: [Float](repeating: -300, count: cells), count: frequencies.count)
        var low = [[Float]](repeating: [Float](repeating: 10, count: cells), count: frequencies.count)
        var high = low
        await withTaskGroup(of: (Int, [[Float]]).self) { group in
            for bearing in 0..<bearings {
                group.addTask(priority: .userInitiated) {
                    if Task.isCancelled { return (bearing, []) }
                    return (bearing, spoke(request, bearing: Double(bearing) * 360 / Double(bearings),
                                           frequencies: frequencies, rings: rings, ringKm: ringKm,
                                           settings: settings, coefficients: coefficients))
                }
            }
            for await (bearing, values) in group where !values.isEmpty {
                for index in 0..<frequencies.count {
                    for ring in 0..<rings {
                        let cell = bearing * rings + ring
                        snr[index][cell] = values[index][ring * 3]
                        low[index][cell] = values[index][ring * 3 + 1]
                        high[index][cell] = values[index][ring * 3 + 2]
                    }
                }
            }
        }
        if Task.isCancelled { return nil }
        return VoacapArea(snr: snr, low: low, high: high)
    }
}
