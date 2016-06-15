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
package org.biojava.nbio.structure.domain.pdp;

public class PDPParameters {
	public static final int MAXLEN = 3200;
	public static final int MAXDOM = 30;
	public static final int MAX_CUTS = 80;
	public static final int MAXSIZE = 350;
	public static final int MIN_DOMAIN_LENGTH = 35;
	public static final int ENDS = 12 ;
	public static final int ENDSEND = 9;
	public static final float RG1 = 1.0f;
	public static final float RG  = 0.0f;
	public static final float TD1 = 40.f;
	public static final float TD  = 25.f;
	public static final float DBL = .05f;
	public static final int MAXCONT = 900;

	public static final float CUT_OFF_VALUE=.50f; /* decide to cut */
	public static final float CUT_OFF_VALUE1=.29f; /* decide to combine */
	public static final float CUT_OFF_VALUE2=.44f; /* decide to double cut */
	public static final float CUT_OFF_VALUE1S=.19f; /* decide to combine small domains */
	public static final float CUT_OFF_VALUE1M=.21f; /* decide to combine medium domains */


}
