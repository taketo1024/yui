# yui (結)

## Khovanov homology (`kh`)

```
Usage: yui kh <NAME> [OPTIONS]

Arguments:
  <NAME>  

Options:
  -l, --link <LINK>        
  -c, --c-value <C_VALUE>  [default: 0]
  -t, --c-type <C_TYPE>    [default: i64]
      --reduced            
      --bigraded           
```

### Example 1 
  
```sh
$ cargo run -r -- kh 3_1 --bigraded

j\i  -3  -2     -1  0 
 -1                  Z 
 -3                  Z 
 -5       Z           
 -7       (Z/2)       
 -9   Z               
```

Note: Any link can be given by the pd-code via `-l [[e1, e2, e3, e4], ...]` option. 

### Example 2 (reduced)

```sh
$ cargo run -r -- kh 3_1 --bigraded --reduced

j\i  -3  -2  -1  0 
 0                 
 -2               Z 
 -4                
 -6       Z        
 -8   Z            
```

### Example 3 (Lee homology)

```sh
$ cargo run -r -- kh 3_1 -c 2

H[-3]: 0
H[-2]: (Z/2)²
H[-1]: 0
H[0]: Z²
```

Note: Lee = `c = 2`, Bar-Natan = `c = 1`.

### Example 4 (various coefficients)

```sh
$ cargo run -r -- kh 3_1 -c "(1,1)" -t gauss

H[-3]: 0
H[-2]: (Z[i]/1 + i)²
H[-1]: 0
H[0]: Z[i]²
```

## Licence
TODO
