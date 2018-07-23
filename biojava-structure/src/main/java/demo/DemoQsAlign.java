/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.quaternary.QsAlign;
import org.biojava.nbio.structure.align.quaternary.QsAlignParameters;
import org.biojava.nbio.structure.align.quaternary.QsAlignResult;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;

/**
 * Demo on how to use programatically {@link QsAlign} for the alignment of
 * quaternary structures.
 * <p>
 * Small oligomers: proliferating cell nuclear antigens (1PLR, 3HI8, 3IFV),
 * photosynthetic reaction centers (2JIY, 1DXR)
 * <p>
 * Big oligomers: cytochrome bc1 complexes (1bcc, 1kb9, 1qcr), phycocyanin
 * (2VML, 2BV8), bacterial ribosome (1FJG, 4V54).
 * 
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class DemoQsAlign {

	public static void main(String[] args) throws IOException,
			StructureException {

		// Align two trimeric DNA clamps
		Structure s1 = StructureIO.getStructure("1bcc");
		Structure s2 = StructureIO.getStructure("1kb9");

		// Select the parameters for clustering and alignment
		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QsAlignParameters alignParams = new QsAlignParameters();

		QsAlignResult result = QsAlign
				.align(s1, s2, clusterParams, alignParams);

		System.out.println(result);

	}
}
