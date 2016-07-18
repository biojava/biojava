package org.biojava.nbio.structure.align.quaternary;

import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.cluster.Subunit;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.cluster.SubunitExtractor;

/**
 * Quaternary Structure Alignment (QS-Align). The algorithm takes as input two
 * protein structures at the quaternary structure level (multiple interacting
 * chains) and calculates the equivalent cross chains and the optimal
 * superposition of the complexes, together with alignment quality scores.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class QsAlign {

	public static QsAlignResult align(Structure s1, Structure s2,
			SubunitClustererParameters cParams, QsAlignParameters aParams) {
		return align(SubunitExtractor.extractSubunits(s1, cParams),
				SubunitExtractor.extractSubunits(s2, cParams), cParams, aParams);
	}

	public static QsAlignResult align(List<Subunit> s1, List<Subunit> s2,
			SubunitClustererParameters cParams, QsAlignParameters aParams) {
		
		return null;
	}

}
