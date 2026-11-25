import XCTest
@testable import Voacap

/// Runs the engine against output from the original Fortran VOACAP, row by row.
final class ReferenceTests: XCTestCase {
    /// rows[hour][name] = the ten printed columns (MUF column first).
    private func parse(_ name: String) throws -> [Int: [String: [String]]] {
        let url = try XCTUnwrap(Bundle.module.url(forResource: name, withExtension: "out", subdirectory: "Reference"))
        let text = try String(contentsOf: url, encoding: .isoLatin1)
        let names = ["FREQ", "TANGLE", "DELAY", "V HITE", "MUFday", "LOSS", "DBU", "S DBW", "N DBW", "SNR LW", "SNR UP",
                     "SNR", "RPWRG", "REL", "SIG LW", "SIG UP"]
        var rows: [Int: [String: [String]]] = [:]
        var hour: Int?
        for line in text.components(separatedBy: .newlines) {
            let trimmed = line.trimmingCharacters(in: .whitespaces)
            guard let key = names.first(where: { trimmed.hasSuffix(" " + $0) }) else { continue }
            let characters = Array(line)
            guard characters.count >= 56 else { continue }
            if key == "FREQ" {
                hour = Int(Float(String(characters[0..<6]).trimmingCharacters(in: .whitespaces)) ?? 0)
            }
            guard let hour else { continue }
            let columns = stride(from: 6, to: 56, by: 5).map { String(characters[$0..<$0 + 5]).replacingOccurrences(of: " ", with: "") }
            if rows[hour]?[key] == nil { rows[hour, default: [:]][key] = columns }
        }
        return rows
    }

    private func antenna(_ name: String) throws -> VoacapAntenna {
        let url = try XCTUnwrap(Bundle.module.url(forResource: name, withExtension: "dat", subdirectory: "Reference"))
        return try XCTUnwrap(VoacapAntenna(gainFile: String(contentsOf: url, encoding: .isoLatin1)))
    }

    /// Compare every row of every hour. Returns (values compared, mismatches).
    private func compare(_ reference: String, settings: VoacapSettings, tx: (Float, Float), rx: (Float, Float),
                         month: Int, ssn: Float, frequencies: [Float], includeMUFColumn: Bool) throws -> (Int, [String]) {
        let rows = try parse(reference)
        XCTAssertEqual(rows.count, 24)
        let engine = VoacapEngine(coefficients: try XCTUnwrap(VoacapCoefficients.bundled))
        engine.apply(settings)
        engine.setCircuit(txLat: tx.0, txLon: tx.1, rxLat: rx.0, rxLon: rx.1, month: month, ssn: ssn)
        var compared = 0
        var mismatches: [String] = []
        for hour in 1...24 {
            _ = engine.predict(hour: hour, frequencies: frequencies)
            let values: [(String, (Int) -> Float, Float)] = [
                ("TANGLE", { engine.angle[$0] }, 0.15), ("DELAY", { engine.delay[$0] }, 0.15),
                ("V HITE", { engine.vhigh[$0] }, 1), ("MUFday", { engine.cprob[$0] }, 0.015),
                ("LOSS", { engine.dblos[$0] }, 1), ("DBU", { engine.dbu[$0] }, 1), ("S DBW", { engine.dbw[$0] }, 1),
                ("N DBW", { engine.xnynois[$0] + engine.rneff[$0] }, 1), ("SNR", { engine.sndb[$0] }, 1),
                ("RPWRG", { engine.snpr[$0] }, 1), ("REL", { engine.reliab[$0] }, 0.015),
                ("SIG LW", { engine.dblosl[$0] }, 0.15), ("SIG UP", { engine.dblosu[$0] }, 0.15),
                ("SNR LW", { engine.snrlw[$0] }, 0.15), ("SNR UP", { engine.snrup[$0] }, 0.15),
            ]
            for (name, value, tolerance) in values {
                let printed = try XCTUnwrap(rows[hour]?[name], "\(reference) hour \(hour) \(name)")
                for column in (includeMUFColumn ? 0 : 1)...frequencies.count {
                    guard let expected = Float(printed[column]) else { continue }
                    let slot = column == 0 ? 12 : column
                    compared += 1
                    if abs(expected - value(slot)) > tolerance {
                        mismatches.append("\(reference) hour \(hour) \(name) column \(column): VOACAP \(expected), port \(value(slot))")
                    }
                }
            }
        }
        return (compared, mismatches)
    }

    /// tests/p2p/test01 from voacapl: Tangier to Belgrade, 2440 km, isotropic antennas.
    func testShortPathAgainstVoacapTest01() throws {
        var settings = VoacapSettings()
        settings.powerWatts = 500_000
        settings.minimumAngle = 0.1
        settings.requiredSNR = 73
        settings.transmitAntenna = .isotropic(gain: 10)
        settings.receiveAntenna = .isotropic(gain: 20)
        settings.magneticPoleLatitude = 79.5
        let (compared, mismatches) = try compare("test01", settings: settings, tx: (35.80, -5.90), rx: (44.90, 20.50), month: 6, ssn: 100,
                                                 frequencies: [6.07, 7.20, 9.70, 11.85, 13.70, 15.35, 17.73, 21.65, 25.89], includeMUFColumn: true)
        XCTAssertGreaterThan(compared, 3500)
        XCTAssertEqual(mismatches, [])
    }

    /// The input cards echoed at the top of a VOACAP output file.
    private struct Deck {
        var tx: (Float, Float) = (0, 0), rx: (Float, Float) = (0, 0)
        var month = 1
        var ssn: Float = 100
        var settings = VoacapSettings()
        var frequencies: [Float] = []
    }

    private func deck(_ url: URL) throws -> Deck {
        var deck = Deck()
        func signed(_ text: String) -> Float {
            let value = Float(text.dropLast()) ?? 0
            return text.hasSuffix("S") || text.hasSuffix("W") ? -value : value
        }
        for line in try String(contentsOf: url, encoding: .isoLatin1).components(separatedBy: .newlines).prefix(40) {
            let card = line.hasPrefix(" ") ? String(line.dropFirst()) : line
            let fields = card.split(separator: " ").map(String.init)
            switch fields.first {
            case "CIRCUIT" where fields.count >= 5:
                deck.tx = (signed(fields[1]), signed(fields[2]))
                deck.rx = (signed(fields[3]), signed(fields[4]))
            case "MONTH" where fields.count >= 3: deck.month = Int(Float(fields[2]) ?? 1)
            case "SUNSPOT" where fields.count >= 2: deck.ssn = Float(fields[1]) ?? 100
            case "SYSTEM" where fields.count >= 6:
                deck.settings.manMadeNoise = Int(Float(fields[2]) ?? 145)
                deck.settings.minimumAngle = Double(fields[3]) ?? 3
                deck.settings.requiredReliability = Int(Float(fields[4]) ?? 90)
                deck.settings.requiredSNR = Double(fields[5]) ?? 0
            case "ANTENNA" where fields.count >= 3 && fields[1] == "1":
                deck.settings.powerWatts = (Double(fields.last ?? "") ?? 1) * 1000
            case "FREQUENCY":
                let characters = Array(card)
                deck.frequencies = stride(from: 10, to: min(characters.count, 65), by: 5).compactMap {
                    Float(String(characters[$0..<min($0 + 5, characters.count)]).trimmingCharacters(in: .whitespaces))
                }.filter { $0 > 0 }
            default: break
            }
        }
        return deck
    }

    /// Any further reference circuits found in Reference/local, each a
    /// `name.out` with the `name.gain01.dat` and `name.gain02.dat` VOACAP
    /// wrote for it. These are not in the repository; the developer's set
    /// covers the long-path model and the 7000–10000 km smoothing.
    func testLocalCircuitsAgainstVoacapl() throws {
        let urls = Bundle.module.urls(forResourcesWithExtension: "out", subdirectory: "Reference/local") ?? []
        try XCTSkipIf(urls.isEmpty, "No local reference circuits")
        for url in urls {
            let name = url.deletingPathExtension().lastPathComponent
            var deck = try deck(url)
            deck.settings.transmitAntenna = try antenna("local/\(name).gain01")
            deck.settings.receiveAntenna = try antenna("local/\(name).gain02")
            let (compared, mismatches) = try compare("local/\(name)", settings: deck.settings, tx: deck.tx, rx: deck.rx, month: deck.month,
                                                     ssn: deck.ssn, frequencies: deck.frequencies, includeMUFColumn: false)
            XCTAssertGreaterThan(compared, 3000, name)
            XCTAssertEqual(mismatches, [], name)
        }
    }

    /// The built-in dipole reproduces a gain table VOACAP wrote for the same antenna.
    func testDipoleAgainstVoacapGainTable() throws {
        guard let url = Bundle.module.urls(forResourcesWithExtension: "dat", subdirectory: "Reference/local")?.first(where: { $0.lastPathComponent.hasSuffix("gain01.dat") }) else {
            throw XCTSkip("No local gain tables")
        }
        let reference = try XCTUnwrap(VoacapAntenna(gainFile: String(contentsOf: url, encoding: .isoLatin1)))
        let dipole = VoacapAntenna.horizontalDipole(offAzimuth: 0, height: -0.5)
        for index in reference.gain.indices where index % 30 >= 1 {
            XCTAssertEqual(reference.gain[index], dipole.gain[index], accuracy: 0.002)
        }
    }

    func testAreaRunsEverywhere() async throws {
        let coefficients = try XCTUnwrap(VoacapCoefficients.bundled)
        for request in [VoacapArea.Request(latitude: 90, longitude: 0, month: 1, ssn: 0, hour: 0),
                        VoacapArea.Request(latitude: -45.2, longitude: 179.9, month: 7, ssn: 250, hour: 13)] {
            let computed = await VoacapArea.compute(request, frequencies: [2, 7.1, 14.1, 28.3], bearings: 24, rings: 50, ringKm: 400, coefficients: coefficients)
            let area = try XCTUnwrap(computed)
            XCTAssertTrue(area.snr.allSatisfy { $0.allSatisfy(\.isFinite) })
        }
    }
}
