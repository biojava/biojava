package org.biojava.dasobert.das2.io;

import java.io.InputStream;
import org.biojava.dasobert.dasregistry.DasCoordinateSystem;

public interface DASRegistryCoordinatesReader {

    /** read a DAS2 coordinates response and return a list of coordinate systems.
	 *
	 */
	public DasCoordinateSystem[] readRegistryCoordinates(InputStream stream);

}
