let
  nixpkgs = import <nixpkgs> {  };
in
  with nixpkgs;
  stdenv.mkDerivation {
    name = "arap";
    buildInputs = [
          pkg-config
          cmake eigen
          ceres-solver
          blas
          libGL libGL.dev
          glog
        ] ++ (with xorg; [
          libX11
          libXrandr
          libXinerama
          libXcursor
          xinput
          libXi
          libXext
        ]);

    shellHook = ''
        export EIGEN_DIR="${eigen.out}"
        export CERES_DIR="${ceres-solver.out}"
        export GLOG_DIR="${glog.out}"
    '';
  }
