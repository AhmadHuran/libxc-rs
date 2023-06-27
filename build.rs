// Copyright 2023 Ahmad W. Huran
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use std::env;
use std::path::PathBuf;

fn main() {
    let sysdep = system_deps::Config::new().probe().unwrap();

    let inc:Vec<String> = sysdep
        .all_include_paths()
        .into_iter()
        .map(|b| b.clone().into_os_string().into_string().unwrap())
        .map(|b| format!("-I{}", b))
        .collect();

    let bindings = bindgen::Builder::default()
        .header("libxc.h")
        .clang_arg(&inc[0])
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("libxc-bindings.rs"))
        .expect("Couldn't write bindings!");
}

