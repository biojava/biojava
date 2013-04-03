package org.biojava3.core.sequence.location;

import static org.biojava3.core.sequence.Strand.UNDEFINED;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.Location;
import org.junit.Test;

public class LocationTest {

  @Test
  public void testSubLocations() {
    List<SimpleLocation> expected = Arrays.asList(
      new SimpleLocation(1, 10, Strand.UNDEFINED),
      new SimpleLocation(2, 6, Strand.UNDEFINED),
      new SimpleLocation(4, 5, Strand.UNDEFINED),
      new SimpleLocation(6, 7, Strand.UNDEFINED),
      new SimpleLocation(11, 20, UNDEFINED)
    );

    Location location = new SimpleLocation(1, 20, UNDEFINED,
      new SimpleLocation(1, 10, UNDEFINED,
          new SimpleLocation(2,6,UNDEFINED, new SimpleLocation(4,5,UNDEFINED)),
          new SimpleLocation(6,7,UNDEFINED)
      ),
      new SimpleLocation(11, 20, UNDEFINED)
    );

    List<Location> actual = location.getAllSubLocations();

    assertEquals("Checking sublocations iterate as expected", toStr(expected), toStr(actual));
  }

  @Test(expected=IllegalArgumentException.class)
  public void badLocations() {
    new SimpleLocation(10, 1, Strand.UNDEFINED);
  }

  private <L extends Location> String toStr(List<L> locations) {
    StringBuilder sb = new StringBuilder();
    for(L l: locations) {
      sb.append(l).append("|");
    }
    return sb.toString();
  }
}
