package org.biojava.bio.program.tagvalue;

public interface BoundaryFinder {
  public boolean dropBoundaryValues();
  public boolean isBoundaryStart(Object value);
  public boolean isBoundaryEnd(Object value);
}
