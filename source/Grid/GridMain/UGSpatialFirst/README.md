# UGSpatialFirst

`UGSpatialFirst` is a Uniform Grid implementation whose physical and scratch
arrays use the layout `(i, j, k, variable, block)`. It implements the same
Flash-X Grid interfaces as `Grid/GridMain/UG`.

Select it with the `Grid=UGSpatialFirst` setup variable and the mandatory
`-index-reorder` option. For example:

```sh
python3 bin/setup.py Sod -auto Grid=UGSpatialFirst -index-reorder \
  -maxblocks=1 -nofbs parallelIO=True --2d
```

The index-reorder option keeps Grid API clients consistent with the
spatial-first four-dimensional block pointers returned by `Grid_tile`.
