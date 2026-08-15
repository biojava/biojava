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

import java.io.IOException;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.io.density.DensityMapCache;
import org.biojava.nbio.structure.io.density.DensityMapKind;
import org.biojava.nbio.structure.io.density.DensityMapResult;
import org.biojava.nbio.structure.io.density.NoDensityMapException;

/**
 * Fetches electron density and cryo-EM maps, printing which source answered.
 * <p>
 * The three entries chosen exercise the three outcomes the fallback chain has to
 * handle:
 * <ul>
 * <li><b>1cbs</b> &mdash; an X-ray structure with deposited structure factors,
 * served by the first source tried.</li>
 * <li><b>6hu9</b> &mdash; a cryo-EM structure. Every X-ray source has nothing for
 * it, so the chain resolves the associated EMDB entry instead and picks up the
 * author-recommended contour level along the way.</li>
 * <li><b>4hhb</b> &mdash; deposited in 1984 without structure factors, so no
 * source has anything. This is a normal outcome, not an error, and the exception
 * says which sources were tried and why each declined.</li>
 * </ul>
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DemoFetchElectronDensity {

	/**
	 * @param args ignored
	 * @throws IOException if a server could not be reached at all
	 */
	public static void main(String[] args) throws IOException {
		DensityMapCache cache = new DensityMapCache();
		System.out.println("Caching under: " + cache.getCachePath());
		System.out.println();

		show(cache, "1cbs", DensityMapKind.TWO_FO_FC);
		show(cache, "1cbs", DensityMapKind.FO_FC);
		show(cache, "6hu9", DensityMapKind.AUTO);
		show(cache, "4hhb", DensityMapKind.AUTO);
	}

	private static void show(DensityMapCache cache, String id, DensityMapKind kind) throws IOException {
		System.out.printf("%s (%s)%n", id, kind);
		try {
			DensityMapResult result = cache.getDensityMap(new PdbId(id), kind);
			System.out.printf("    source   : %s%n", result.getSource());
			System.out.printf("    format   : %s%s%n", result.getFormat(),
					result.isRenderable() ? "" : "  (needs an FFT before display)");
			System.out.printf("    kind     : %s%n", result.getKind());
			System.out.printf("    file     : %s (%,d bytes)%n", result.getFile(), result.getFileSizeBytes());
			System.out.printf("    cached   : %s%n", result.isFromCache());
			if (result.getEmdbId() != null) {
				System.out.printf("    EMDB     : %s%n", result.getEmdbId());
			}
			if (result.getRecommendedContourLevel() != null) {
				System.out.printf("    contour  : %s (author recommended)%n", result.getRecommendedContourLevel());
			}
			if (result.getContourInSigma() != null) {
				System.out.printf("    in sigma : %.2f%n", result.getContourInSigma());
			}
		} catch (NoDensityMapException e) {
			System.out.printf("    no map available%n");
			e.getAttempts().forEach((source, reason) ->
					System.out.printf("      %-24s %s%n", source, reason));
		}
		System.out.println();
	}
}
