package org.biojava3.core.sequence.location;

import java.util.Arrays;
import java.util.List;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.AbstractLocation;
import org.biojava3.core.sequence.location.template.Location;

/**
 * Very basic implementation of the Location interface.
 *
 * @author ayates
 */
public class SimpleLocation extends AbstractLocation {

  public SimpleLocation(int start, int end, Strand strand) {
    super(start, end, strand);
  }

  public SimpleLocation(int start, int end, Strand strand, boolean circular,
      List<Location> subLocations) {
    super(start, end, strand, circular, subLocations);
  }

  public SimpleLocation(int start, int end, Strand strand, boolean circular,
      Location... subLocations) {
    super(start, end, strand, circular, subLocations);
  }

  public SimpleLocation(int start, int end, Strand strand, boolean circular) {
    super(start, end, strand, circular);
  }

  public SimpleLocation(int start, int end, Strand strand,
      List<Location> subLocations) {
    super(start, end, strand, subLocations);
  }

  public SimpleLocation(int start, int end, Strand strand,
      Location... subLocations) {
    super(start, end, strand, Arrays.asList(subLocations));
  }
}