{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    wild = {
      url = "github:davidlattimore/wild";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = {
    nixpkgs,
    flake-utils,
    rust-overlay,
    wild,
    ...
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {
        inherit system;
        overlays = [rust-overlay.overlays.default (import wild)];
      };
      # Rustのバージョンを1.91.0に固定
      rustToolchain = pkgs.rust-bin.stable."1.91.1".default.override {
        extensions = ["rust-src" "rust-analyzer" "rustfmt" "clippy"];
      };
    in {
      devShells.default = pkgs.mkShell {
        packages = [
          rustToolchain
          pkgs.wild
          pkgs.gnuplot
        ];
        # 環境変数
        RUSTFLAGS = "-C linker=wild";
        # rust-src のパスを通す: rust-analyzer が標準ライブラリのソースを見つけるために必須
        RUST_SRC_PATH = "${rustToolchain}/lib/rustlib/src/rust/library";
      };
    });
}
