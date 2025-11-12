# ykh

A CLI tool for Khovanov homology computations.

## Setup

* Install [Rust](https://www.rust-lang.org/tools/install).
* Install `ykh` by

```sh
$ cargo install ykh
```

## Usage

Computes the Khovanov homology / complex for the given link. 

```
Usage: ykh <COMMAND> <LINK> [OPTIONS]

Arguments:
  <COMMAND> [kh, ckh, khi, ckhi]
  <LINK>    name or PD-code of the input link.

Options:
  -t, --c-type <C_TYPE>    [default: Z] [possible values: Z, Q, F2, F3]
  -c, --c-value <C_VALUE>  [default: 0]
  -m, --mirror             
  -r, --reduced
  -f, --format <FORMAT>    [default: unicode] [possible values: unicode, tex]
  -h, --help               Print help
```

## Examples

Mode `kh` computes the Khovanov homology of the given input link.
  
```sh
$ ykh kh 3_1
```
```
j\i  -3  -2     -1  0 
 -1                  Z 
 -3                  Z 
 -5       Z           
 -7       (Z/2)       
 -9   Z               
```

`<LINK>` can also be specified by the [PD-code](https://knotinfo.math.indiana.edu/descriptions/pd_notation.html) as:

```sh
$ ykh kh "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]"
```

The following computes the (bigraded) Bar-Natan homology over `Q[H]`.

```
$ ykh kh 3_1 -t Q -c H
```
```
j\i  -2         -1  0 
 -1   .          .   Q[H] 
 -3   .          .   Q[H] 
 -5   (Q[H]/H²)  .   .
```

Mode `ckh` computes the Khovanov complex for the given link, simplified using [Bar-Natan's local algorithm](https://arxiv.org/abs/math/0606318). It can be also computed over the (non-PID) rings `Z[H]`, `Z[T]`, `Z[H, T]`, together with its differential by specifying the option `-d`.

```sh
$ ckh 3_1 -t Z -c H,T -d
```
```
j\i  -3       -2       -1  0 
 -1   .        .        .   Z[H, T] 
 -3   .        .        .   Z[H, T] 
 -5   .        Z[H, T]  .   . 
 -7   Z[H, T]  Z[H, T]  .   . 
 -9   Z[H, T]  .        .   . 

d[-3]: Z[H, T]² -> Z[H, T]²

  ┌         ┐
  │  -H  -2 │
  │ -2T   H │
  └         ┘
```

Mode `khi` and `ckhi` computes the [involutive Khovanov homology / complex](https://arxiv.org/abs/2404.08568) for a strongly invertible knot. Defaults to `-t F2`, which is the only supported type. 

```
$ ykh khi 3_1 
```
```
j\i  0   1   2   3    4 
 9    .   .   .   F₂   F₂ 
 7    .   .   F₂  F₂²  F₂ 
 5    .   .   F₂  F₂   . 
 3    F₂  F₂  .   .    . 
 1    F₂  F₂  .   .    .
```

The following shows the (reduced) involutive Bar-Natan homology for the knot `m(9_46)`. The target knot must start with edge index `1` that lies on the axis of the involution.

```
ykh khi \
"[[18,8,1,7],[13,6,14,7],[12,2,13,1],[8,18,9,17],[5,14,6,15],[2,12,3,11],[16,10,17,9],[15,4,16,5],[10,4,11,3]]" \
-t F2 -c H -r
```
```
j\i  0      1      2          3          4          5          6          7 
 12   .      .      .          .          .          .          (F₂[H]/H)  (F₂[H]/H) 
 10   .      .      .          .          .          .          .          . 
 8    .      .      .          .          (F₂[H]/H)  (F₂[H]/H)  .          . 
 6    .      .      .          (F₂[H]/H)  (F₂[H]/H)  .          .          . 
 4    .      .      .          .          .          .          .          . 
 2    .      F₂[H]  (F₂[H]/H)  .          .          .          .          . 
 0    F₂[H]  .      .          .          .          .          .          .
```

## License
`yui` is released under the [MIT license](LICENSE).
