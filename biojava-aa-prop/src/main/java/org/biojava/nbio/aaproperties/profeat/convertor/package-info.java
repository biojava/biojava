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
/**
 * Set of classes that enable the conversion protein sequences into various attributes.
 * The seven different attributes are<p/>
 * Hydrophobicity (Polar, Neutral, Hydrophobicity)<br/>
 * Normalized van der Waals volume (Range 0 - 2.78, 2.95 - 4.0, 4.03 - 8.08)<br/>
 * Polarity (Value 4.9 - 6.2, 8.0 - 9.2, 10.4 - 13.0)<br/>
 * Polarizability (Value 0 - 1.08, 0.128 - 0.186, 0.219 - 0.409)<br/>
 * Charge (Positive, Neutral, Negative)<br/>
 * Secondary structure (Helix, Strand, Coil)<br/>
 * Solvent accessibility (Buried, Exposed, Intermediate)<br/>
 *
 * @author Koh Chuan Hock
 *
 * @version 2011.08.22
 * @since 3.0.2
 */
package org.biojava.nbio.aaproperties.profeat.convertor;

