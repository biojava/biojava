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
package org.biojava.nbio.structure.symmetry.core;

/**
 * Method used for symmetry perception in the {@link QuatSymmetryDetector}.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public enum SymmetryPerceptionMethod {
	NO_ROTATION,
	C2_ROTATION,
	ROTATION,
	ROTO_TRANSLATION;
}
