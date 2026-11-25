// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "Voacap",
    platforms: [.iOS(.v15), .macOS(.v12), .tvOS(.v15), .watchOS(.v8), .visionOS(.v1)],
    products: [
        .library(name: "Voacap", targets: ["Voacap"]),
    ],
    targets: [
        .target(name: "Voacap", resources: [.copy("Resources/voacap-coeffs.bin")]),
        .testTarget(name: "VoacapTests", dependencies: ["Voacap"], resources: [.copy("Reference")]),
    ]
)
