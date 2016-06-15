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
/**<p>
 * BioJava provide a module biojava-protein-disorder for prediction disordered regions
 * from a protein sequence. Biojava-protein-disorder module for now contains one method
 * for the prediction of disordered regions. This method is based on the Java
 * implementation of RONN predictor.
 * <p>
 * This code has been originally developed for use with <a href="http://www.compbio.dundee.ac.uk/jabaws">JABAWS</a>.
 * We call this code JRONN. JRONN is based on the C implementation of RONN algorithm and uses the same model data,
 * therefore gives the same predictions. JRONN based on <a href="http://www.strubi.ox.ac.uk/RONN">RONN</a>
 * version 3.1 which is still current in time of writing (August 2011).
 * Main motivation behind JRONN development was providing an implementation of RONN more
 * suitable to use by the automated analysis pipelines and web services.
 * Robert Esnouf has kindly allowed us to explore the RONN code and share the results with the community.
 * <p>
 * Original version of RONN is described in Yang,Z.R., Thomson,R., McMeil,P. and Esnouf,R.M. (2005)
 * RONN: the bio-basis function neural network technique applied to the detection of natively
 * disordered regions in proteins. Bioinformatics 21: 3369-3376
 * <p>
 * Examples of use are provided below. For more information please refer to JronnExample testcases.

 * @author Peter Troshin
 *
 * @version 1.0 January 2010
 * @since 3.0.2
 */
package org.biojava.nbio.ronn;

