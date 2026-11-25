# Voacap

A Swift translation of the VOACAP HF propagation prediction engine, for iOS, macOS and the other Apple platforms. It has no dependencies and needs no network: the ionospheric coefficient tables ship inside the package (557 KB).

It follows the voacapl 16.1207W Fortran sources routine by routine and implements VOACAP method 30: the short-path model, the long-path model, and the smoothing between them from 7,000 to 10,000 km.

## Use

```swift
import Voacap

var settings = VoacapSettings()
settings.powerWatts = 100
settings.requiredSNR = 14            // dB in 1 Hz; 14 is FT8
settings.transmitAntenna = .horizontalDipole(offAzimuth: 0, height: -0.5)
settings.receiveAntenna = settings.transmitAntenna

let engine = VoacapEngine(coefficients: VoacapCoefficients.bundled!)
engine.apply(settings)
engine.setCircuit(txLat: 40.03, txLon: -105.30, rxLat: 51.5, rxLon: -0.13, month: 11, ssn: 83)

let prediction = engine.predict(hour: 15, frequencies: [7.1, 14.1, 21.2, 28.3])
print(prediction.muf)
for result in prediction.results {
    print(result.frequency, result.reliability, result.snr)
}
```

- Longitudes are positive east. Hours are UT, 1 to 24. Up to eleven frequencies per call, 2 to 30 MHz.
- One `VoacapEngine` holds all the state for one circuit. It is not thread safe; use one engine per thread. `VoacapCoefficients` and `VoacapAntenna` can be shared.
- `VoacapArea.compute` runs a whole bearing-by-distance grid around a station in parallel and returns the median SNR and its spread for each cell; `VoacapArea.reliability` turns those into a reliability for any required SNR.
- Antennas are gain tables (1–30 MHz by 0–90° elevation). Built in: isotropic, VOACAP's type 23 horizontal dipole, a table from any function, and a reader for VOACAP's `gainNN.dat` files.

## What is and isn't included

Included: method 30 point-to-point predictions with the CCIR coefficients (MUF, mode, take-off angle, loss, signal, noise, SNR, reliability and their deciles).

Not included: the URSI coefficients, the other output methods, service probability, multipath probability, VOACAP's other antenna models, and the input-deck and report formats.

The geomagnetic pole defaults to 78.5°N 69°W, which is what voacapl runs with. The Windows VOACAP distribution uses 79.5°N; set `magneticPoleLatitude` to match it.

## Accuracy

`swift test` runs the engine against `tests/p2p/test01.out` from voacapl (Tangier to Belgrade, 24 hours, nine frequencies). Every printed value must agree to within one unit of its last printed digit.

The port was also checked the same way against seven circuits from 5,600 to 16,000 km run with voacapl 0.7.7, covering the long-path model and the smoothing range, with no differences. Those reference files are not in the repository. To add your own, put `name.out` and the `gain01.dat` / `gain02.dat` VOACAP wrote for that run (renamed `name.gain01.dat`, `name.gain02.dat`) in `Tests/VoacapTests/Reference/local/`; the tests pick them up.

## Licence

The Swift code is dedicated to the public domain under CC0 1.0; see `LICENSE`.

## Origin

VOACAP was developed by NTIA/ITS, an agency of the U.S. Government, and is not subject to copyright in the United States. voacapl is J. A. Watson's port to Linux, whose changes are released under CC0. See `NOTICE-VOACAP.txt`. `make-coefficients.py` rebuilds `voacap-coeffs.bin` from voacapl's `coeffNN.asc` and `fof2CCIR.asc`.
