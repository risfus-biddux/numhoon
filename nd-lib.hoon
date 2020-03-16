!:
:: in dojo, to load library into subject: /+  nd-lib
:: in dojo, to create an nd-array: =my-nd (make-nd-array:nd-lib flat=(gulf 0 23) dims=~[4 3 2])
:: todo: when using =/ to initialize as the bunt value, just use =|
:: todo: access-value should wrap something like get-value-index
:: so in matrix-mult we can call get-value-index without having to pass everything
:: todo: ohhh i don't think you need to limo if you have already declared its type as a list.
:: actually... at least in defining strides we do need to.. so no.
|%
+$  nd-array  [flat=(list @) dims=(list @) strides=(list @)]     :: define type
:: dot product of plain list right now since we just need for calcs
:: not of some kind of 1-dimensional nd-array
::
++  dot-product
  |=  [a=(list @) b=(list @)]
  ^-  @ 
  (roll (turn (zip a b) mul-list) add)
++  zip                                       :: like haskell
  |=  [a=(list @) b=(list @)]
  ?>  =((lent a) (lent b))
  |-
  ?~  a  ~
  [i=~[-:a -:b] t=$(a +:a, b +:b)]
++  mul-list                                  :: multiply two elements of a 2-list
  |=  l=(list @)
  ^-  @
  ?>  =((lent l) 2)
  (mul -:l +<:l)
::  need to combine row and col
::  and additionally, can probably find a math trick to
::  not need to pull individually. would need to refactor the constructor?
::
++  get-row
  |=  [nd=nd-array index=@]
  ^-  (list @)
  ?>  =((lent dims.nd) 2)
  ?>  (lth index (snag 0 dims.nd))
  =/  list-out=(list @)  ~
  =+  curr-col=0
  |-  ^-  (list @)
  ?:  =(curr-col (snag 1 dims.nd))
    (flop list-out)
  %=  $
    list-out  [(access-value nd ~[index curr-col]) list-out]
    curr-col  +(curr-col)
  ==
++  get-col
  |=  [nd=nd-array index=@]
  ^-  (list @)
  ?>  =((lent dims.nd) 2)
  ?>  (lth index (snag 1 dims.nd))
  =/  list-out=(list @)  ~
  =+  curr-row=0
  |-  ^-  (list @)
  ?:  =(curr-row (snag 0 dims.nd))
    (flop list-out)
  %=  $
    list-out  [(access-value nd ~[curr-row index]) list-out]
    curr-row  +(curr-row)
  ==
++  access-value
  :: access value of an nd-array from indices
  ::
  |=  [nd=nd-array indices=(list @)]
  =<
  ?>  (check-boundedness indices dims.nd)
  ^-  @
  (snag (dot-product indices strides.nd) flat.nd)
  |%
  ++  check-boundedness                         :: is each index less than its dim?
    |=  [indices=(list @) dims=(list @)]
    ?>  =((lent indices) (lent dims))
    |-  ^-  ?
    ?~  indices  %.y
    ?:  (lth -:indices -:dims)
      $(indices +:indices, dims +:dims)
    %.n
  --
++  make-nd-array
  :: constructor, which also pre-calculates strides.
  :: note that these strides are just index-based, not memory-based as in numpy.
  :: in numpy, say you use a dtype with size 4 bytes, then the strides would be
  :: 4x what they are here.
  ::
  |=  [flat=(list @) dims=(list @)]
  ?>  =(~ (find ~[0] dims))
  ?>  =((lent flat) (roll dims mul))
  =/  dims-calc=(list @)  +:dims                :: can't do these...
  =.  dims-calc  (flop dims-calc)               :: ...in the same line
  =/  strides=(list @)  (limo ~[1])
  :: build strides list
  ::
  |-  ^-  nd-array
  ?~  dims-calc
    [flat=flat dims=dims strides=strides]
  %=  $
    strides  [(mul -:dims-calc -:strides) strides]
    :: throw away head of dims-calc. not used in product for later indices
    ::
    dims-calc  +:dims-calc
  ==
++  scalar-mult
  |=  [nd=nd-array scalar=@]
  =.  flat.nd
  (turn flat.nd |=(a/@ (mul a scalar)))
  nd
++  matrix-mult
  :: 2d only
  ::
  |=  [x=nd-array y=nd-array]
  ?>  &(=((lent dims.x) 2) =((lent dims.y) 2))
  ?>  =((snag 1 dims.x) (snag 0 dims.y))
  =/  flat-out=(list @)  ~
  =/  dims-out=(list @)  ~[(snag 0 dims.x) (snag 1 dims.y)]
  =+  current-row=0
  |-  ^-  nd-array
  ?:  =(current-row (snag 0 dims-out))
  ::  actual line we return after we're done with last row
  ::
    (make-nd-array flat=(flop flat-out) dims=dims-out)
  =+  current-col=0
  =/  flat-out-calc=(list @)  ~
  =.  flat-out-calc
    |-  ^-  (list @)
    ?:  =(current-col (snag 1 dims-out))
      flat-out-calc
    =/  x-row=(list @)  (get-row x current-row)
    =/  y-col=(list @)  (get-col y current-col)
    %=  $
      flat-out-calc  [(dot-product x-row y-col) flat-out-calc]
      current-col  +(current-col)
    ==
  %=  $
    flat-out  (weld flat-out-calc flat-out)
    current-row  +(current-row)
  ==
--
