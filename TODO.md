# TODO LIST:

-   add enterotype to the network of metabolites and OTU

-   add food/diet module and link it to the gut OTU and gut microbiota metabolites module;

    > ref: Noronha, A., Modamio, J., Jarosz, Y., Guerard, E., Sompairac, N., Preciat, G., ... & Thiele, I. (2019). The Virtual Metabolic Human database: integrating human and gut microbiome metabolism with nutrition and disease. *Nucleic acids research*, *47*(D1), D614-D624.
    >
    > website: [Virtual Metabolic Human (vmh.life)](https://www.vmh.life/)

-   Establish three modules and their relationship: gut microbiota (OTU); gut microbiota metabolites; food/diet;

### API:

1.  microbe
2.  metabolites
3.  gene
4.  food items

### class:

vmh:

1.  type (character): microbe, metabolites, gene, food items

2.  net (data.table):

    -   microbe
    -   metabolites
    -   gene
    -   food items
    -   reaction

### methods

corr:

1.  mSet, bSet
2.  network

?
