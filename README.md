# yui (結)

## Setup

* Install [Rust](https://www.rust-lang.org/tools/install).
* Download [`yui`](https://github.com/taketo1024/yui.git).
* Run.

```sh
$ cargo run -r -- [command] [args]
```

### Optional features
Build with `--features [features]` to enable optional features. 
Available options:

* `bigint`: enables `BigInt` to avoid overflow. 
* `poly`  : enables polynomial types.
* `qint`  : enables quadratic integer types (Gaussian and Eisenstein integers)

e.g. 

```sh
$ cargo run -r --features bigint,poly -- [command] [args]
```

## Commands
### Khovanov homology (`kh`)

Computes the Khovanov homology for the given link. 

```
Usage: kh <LINK> [OPTIONS]

Arguments:
  <LINK>  

Options:
  -c, --c-value <C_VALUE>  [default: 0]
  -t, --c-type <C_TYPE>    [default: Z]
  -m, --mirror
  -r, --reduced            
  -b, --bigraded (available only when c = 0.)
```

#### Examples
  
```sh
$ cargo run -r -- kh 3_1 -b
```
```
j\i  -3  -2     -1  0 
 -1                  Z 
 -3                  Z 
 -5       Z           
 -7       (Z/2)       
 -9   Z               
```

`<LINK>` can be specified by the [pd-code](https://knotinfo.math.indiana.edu/descriptions/pd_notation.html) as:

```sh
$ cargo run -r -- kh "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]" -b
```

The following computes the (bigraded) Bar-Natan homology over `Q[H]`.

```
$ cargo run --features poly -- kh 3_1 -t Q -c H
```
```
[-3]: 0
[-2]: (Q[H]/H²)
[-1]: 0
[0]: Q[H]²
```

### Khovanov complex (`ckh`)

```
Usage: ckh [OPTIONS] <LINK>

Arguments:
  <LINK>  

Options:
  -c, --c-value <C_VALUE>  [default: 0]
  -t, --c-type <C_TYPE>    [default: Z] [possible values: Z, Q, F2, F3, Gauss, Eisen]
  -m, --mirror             
  -r, --reduced            
  -b, --bigraded           
```

#### Examples

```sh
$ cargo run -r -- ckh 3_1 -t Z -b
```
```
j\i  -3  -2  -1  0 
 -1               Z 
 -3               Z 
 -5       Z        
 -7   Z   Z        
 -9   Z            

C[(-3, -9)]: Z -> 0
[[]]

C[(-3, -7)]: Z -> Z
[[-2]]

...
```

`ckh` can be computed over the ring `Z[H]`.

```sh
$ cargo run -r --features poly -- ckh 3_1 -t Z -c H
```
```
C[-3]: Z[H]² -> Z[H]²
[[-H, -2],
 [0, H]]

C[-2]: Z[H]² -> 0
[[]]

C[-1]: 0 -> Z[H]²
[[]]

C[0]: Z[H]² -> 0
[[]]
```

### ss-invariant (`ss`)

See [this paper](https://arxiv.org/abs/2211.02494) for detail.

```
Usage: yui ss [OPTIONS] --c-value <C_VALUE> <LINK>

Arguments:
  <LINK>  

Options:
  -c, --c-value <C_VALUE>  
  -t, --c-type <C_TYPE>    [default: Z] [possible values: Z, Q, F2, F3, Gauss, Eisen]
  -m, --mirror             
  -r, --reduced            
```

The following computes the [s-invariant](https://knotinfo.math.indiana.edu/descriptions/rasmussen_invariant.html) over `Q`.

```
$ cargo run --features poly -- ss 3_1 -t Q -c H
-2
```

## License
`yui` is released under the [MIT license](LICENSE).
