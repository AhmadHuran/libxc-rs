# libxc-rs
Wraps around bindgen-generated rust bindings for [libxc](https://tddft.org/programs/libxc/).

Currently, only LDA and GGA functionals with exc vrho and vsigma support.

### Requierments
`libxc >= 6.6.2`

Must be visible to `pkg-config`. If not, compile and install libxc then add to your `.bashrc`

```
export PKG_CONFIG_PATH=/path/to/libxc/pkgconfig/
```
### Installation
In the root of your crate execute

```
cargo add --git https://github.com/AhmadHuran/libxc-rs.git
```

### To-do
- [ ] Add covering tests
- [ ] Add examples
- [ ] Possibly separate the raw bindings as a  -sys crate
