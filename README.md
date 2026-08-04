# yui (結)

`yui` is a collection of libraries for homology computations, with particular focus on *knot homology theories*, written in [Rust](https://www.rust-lang.org).

## Libraries

- [`yui`](./yui/) - The core library.
- [`yui-matrix`](./yui-matrix/) - Sparse matrix computations, based on [nalgebra](https://nalgebra.org).
- [`yui-homology`](./yui-homology/) - Generic homology computations.
- [`yui-link`](./yui-link/) - Knots and links.
- [`yui-kh`](./yui-khovanov/) - Khovanov homology.

## Binaries
- [`ykh`](./bin-ykh/) - Khovanov homology computations.

## Knot/link data

`Link::load(name)` and `Braid::load(name)` (and the `ykh` CLI when given a name like `3_1`) read from an external user-data directory:

1. `$YUI_DATA_DIR` if set.
2. Otherwise the platform user-data dir:
    - macOS: `~/Library/Application Support/yui/`
    - Linux: `${XDG_DATA_HOME:-~/.local/share}/yui/`
    - Windows: `%APPDATA%\yui\`

Inside that dir, data is partitioned by kind, e.g. `<data_dir>/links/3_1.json`, `<data_dir>/braid/3_1.json`.

To populate it from the [KnotInfo](https://knotinfo.org/) database:

```bash
python3 scripts/fetch-knotinfo-data.py            # writes to the default data dir
python3 scripts/fetch-knotinfo-data.py --out DIR  # or to a custom directory
```

## License
`yui` is released under the [MIT license](LICENSE).
