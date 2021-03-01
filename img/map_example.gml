graph [
  directed 1
  node [
    id 0
    label "ENSG00000251380"
    version 37
    coord "chr5:135444214-135447348(-)"
    color "#1a9850"
    status "latest"
  ]
  node [
    id 1
    label "ENSG00000253827"
    version 5
    coord "chr5:134779908-134788089(-)"
    color "#fc8d59"
    status "mapped"
  ]
  node [
    id 2
    label "ENSG00000255833"
    version 37
    coord "chr5:135444226-135452351(-)"
    color "#1a9850"
    status "latest"
  ]
  node [
    id 3
    label "ENSG00000255290"
    version 6
    coord "chr5:134779908-134788089(-)"
    color "#fc8d59"
    status "mapped"
  ]
  node [
    id 4
    label "ENSG00000254763"
    version 6
    coord "chr5:134779908-134788089(-)"
    color "#fc8d59"
    status "mapped"
  ]
  edge [
    source 1
    target 0
    label "(5>>6)(0.38>>1.00)"
  ]
  edge [
    source 1
    target 3
    label "(5>>6)(1.00>>1.00)"
  ]
  edge [
    source 1
    target 4
    label "(5>>6)(1.00>>1.00)"
  ]
  edge [
    source 3
    target 0
    label "(6>>7)(0.38>>1.00)"
  ]
  edge [
    source 3
    target 2
    label "(6>>7)(1.00>>1.00)"
  ]
  edge [
    source 4
    target 0
    label "(6>>7)(0.38>>1.00)"
  ]
  edge [
    source 4
    target 2
    label "(6>>7)(1.00>>1.00)"
  ]
]
