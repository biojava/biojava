/**
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
package org.biojava.bio.program.sax;

import java.util.HashMap;
import java.util.StringTokenizer;

import org.xml.sax.SAXException;

/**
 * A Helper class for parsing summary lines of Blast-like
 * output. For example:
 * <p>
 * GB_PR9:HUMPSULTRA  Begin: 1 End: 213
 * !L25275 Human estrogen sulfotransferase mRNA, ...          381  e-103
 * <p>
 * The first token from the left side of the first line becomes
 * the start of the description.
 * For GCG, 2 tokens are taken from the right hand side of
 * the second line.  In the above case, the 2 right-most tokens
 * are 381 and e-103.
 * <p>
 * Primary author - <ul>
 * <li>Mathieu Wiepert (Mayo Foundation)
 * </ul>>
 * Copyright &copy; 2001 Mayo Foundation
 *
 *
 * This code released to the biojava project, April 2001
 * under the LGPL license.
 *
 * @author Mayo Foundation
 * @author Greg Cox
 * @version 0.1
 *
 */
final class GCGBlastSummaryLineHelper implements SummaryLineHelperIF {
    private StringBuffer oHitDescription;
    private String previousHitID;
    private String previousScore;
    private String previousEValue;
    public GCGBlastSummaryLineHelper() {
    }

    /**
     * GCG Summary lines come in pairs and look like this
     *
     * GB_PR4:AK027092  Begin: 685 End: 988
     * !AK027092 Homo sapiens cDNA: FLJ23439 fis, clone...    504  e-140
     *
     * Because of the two lines per hit, there are two states
     * for the helper to deal with.
     *
     * They are parsed according to the following rules:
     *
     * The first line becomes the beginning of the description.
     *
     * The second line is handled like this:
     * From the left, tokenizing on white space, extracting the first
     * token (above, this would be "!AK027092") and places it as a
     * String in the object Buffer.
     *
     * From the right, and tokenizing on white space, looks for
     * a specified number of tokens, which it places as
     * Strings in the object map
     *
     * @param poLine     -
     * @param poMap  A HashMap of name-value pairs to be
     * be interpreted by the calling class. The first two
     * items in the map will be the HitId and the HitDescription.
     * Subsequent will be attribute name-values pairs such as
     * Score, E-value.
     */
    public void parse(String poLine, HashMap poMap,
              BlastLikeVersionSupport poVersion) throws SAXException {

        int iGrab = 2; //number of tokens to take from the right
        int iCount;
        boolean firstLine; //state variable
        if (poLine.startsWith("!")) {
            firstLine = false;
        } else {
            firstLine = true;
            oHitDescription= new StringBuffer();
        }
        StringTokenizer oSt = new StringTokenizer(poLine);

        //GCG-blast all flavors - two tokens:
        //first is score
        //next is Evalue
        //These tokens are on the seocnd line.
        //NOTE:  For GCG, if the next match is for the same
        //Accession number then, the score and evalue are not
        //repeated

        if (!firstLine) {
            //populate Map...
            iCount = oSt.countTokens() - iGrab - 1;
            //first token is the hit id
            String hitID = oSt.nextToken();
            poMap.put("hitId",hitID);
            //oHitDescription.setLength(0);

            for (int i = 0; i < iCount; i++) {
                oHitDescription.append(oSt.nextToken());
                oHitDescription.append(" ");
            }
            poMap.put("hitDescription",oHitDescription.substring(0));

            //If the Previous hitID is the same as this, then
            //GCG does not repeat the evalue.  So use the previous e-value
            //now collect score and e-value
            if (hitID.equals(previousHitID)) {
            System.out.println("curr: " + hitID + " prev: " + previousHitID);
	            poMap.put("score",previousScore);
    	        poMap.put("expectValue",previousEValue);
    	    } else {
    	    	previousHitID = hitID;
    	    	previousScore = oSt.nextToken();
    	    	previousEValue = oSt.nextToken();
            	poMap.put("score",previousScore);
            	poMap.put("expectValue",previousEValue);
            }

        } else {
            //initalize description
            oHitDescription.setLength(0);
            oHitDescription.append(oSt.nextToken());
            oHitDescription.append(" ");
        }
    }
}

