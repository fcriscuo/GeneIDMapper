# GeneIDMapper
A utility to resolve gene name from chromosome coordinates. Based on RxJava

Utilizes an early release of David Moten's rtree API. This API has been significantly 
modified since the GeneIDMapper application was originally developed. 

Utilizes a early version of RxJava. The Maven pom file still references the original Netflix repository.

Creates an R-tree data structure where known genes are represented by zero-height Rectangles based on the chromosome number, 
start position, and length. Start & end positions for genes on negative strands are reversed. All locations are normalized to
a positive strand orientation. Given a variant's location (i.e. chromosome and position), the application will resolve the 
gene name associated with that position. Positions that cannot be resolved to a gene are considered intergenic.

The application does not take advantage of Java 8 capabilities.

TODO: refactor to Java 8
      refactor to the latest RxJava release (https://github.com/ReactiveX/RxJava)
      refactor to the latest rtree release (https://github.com/davidmoten/rtree)
      update gene location maps