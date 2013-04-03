package org.biojava3.core.sequence.location.template;

import java.util.List;

import org.biojava3.core.sequence.Strand;

public interface Location extends Iterable<Location> {

  int getMin();

  int getMax();

  Strand getStrand();

  List<Location> getSubLocations();

  boolean isComplex();

}
