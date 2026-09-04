/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These
 * should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the
 * biojava-l mailing list, visit the home page at:
 *
 * http://www.biojava.org/
 */
package demo;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.io.density.DensityMapCache;
import org.biojava.nbio.structure.io.density.DensityMapKind;
import org.biojava.nbio.structure.io.density.DensityMapResult;

/**
 * Shows 1CBS with its electron density drawn around the bound retinoic acid.
 * <p>
 * Both maps are displayed: the 2mFo-DFc map in blue at 1 sigma, which should hug
 * the ligand closely, and the mFo-DFc difference map as a red and green pair at
 * 3 sigma, which for a well refined structure should show very little.
 * <p>
 * Once the window is up, the map can be manipulated from the Rasmol command box
 * at the bottom, for instance
 * <pre>
 * isosurface ID "bj_density_2fofc" delete
 * </pre>
 * to remove just the blue surface. Pressing Reset Display keeps the maps, since
 * they are folded into the saved state when they are drawn.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DemoShowElectronDensity {

	/**
	 * @param args an optional PDB ID to display instead of the default
	 * @throws Exception if the structure or the map could not be fetched
	 */
	public static void main(String[] args) throws Exception {
		String id = args.length > 0 ? args[0] : "1cbs";

		Structure structure = StructureIO.getStructure(id);
		StructureAlignmentJmol viewer = new StructureAlignmentJmol();
		viewer.setStructure(structure);
		viewer.evalString("select all; cartoon on; color chain; "
				+ "select ligand; wireframe 0.16; spacefill 0.4; color cpk;");

		DensityMapCache cache = new DensityMapCache();
		System.out.println("Density cache: " + cache.getCachePath());

		// Clipping to the ligand keeps the surface readable and, for a large map,
		// keeps the contouring quick enough not to stall the interface.
		for (DensityMapKind kind : new DensityMapKind[] {DensityMapKind.TWO_FO_FC, DensityMapKind.FO_FC}) {
			cache.findDensityMap(new PdbId(id), kind).ifPresent(map -> {
				System.out.printf("%-8s from %-20s %s (%,d bytes)%n",
						map.getKind(), map.getSource(), map.getFile().getName(), map.getFileSizeBytes());
				viewer.getJmolPanel().loadDensityMap(map, "{ligand}", 5.0);
			});
		}

		DensityMapResult any = cache.findDensityMap(new PdbId(id), DensityMapKind.AUTO).orElse(null);
		if (any == null) {
			System.out.println("No density is available for " + id
					+ " - try an entry with deposited structure factors, such as 1cbs.");
		}
	}
}
