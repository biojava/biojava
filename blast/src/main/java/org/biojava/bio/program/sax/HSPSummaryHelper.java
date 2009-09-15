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
package org.biojava.bio.program.sax;

import java.util.HashMap;
import java.util.StringTokenizer;

import org.xml.sax.SAXException;


/**
 * A Helper class for parsing HSP summary sections...
 * e.g. things like
 *
 * Score =  326 bits (826), Expect = 2e-89
 *
 * Should be passed in a comma separated list.
 *
 * Primary author -
 *                 Simon Brocklehurst (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Colin Hardman      (CAT)
 *                 Stuart Johnston    (CAT)
 *                 Mathieu Wiepert    (Mayo Foundation)
 *
 * Copyright 2000 Cambridge Antibody Technology Group plc.
 * 
 *
 * This code released to the biojava project, May 2000
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @version 0.1
 *
 */
final class HSPSummaryHelper {

    private HSPSummaryHelper() {
    }

    /**
     * Takes a comma separated String of Blast-like HSP summary
     * information e.g.
     *
     * Score =  210 bits (454), Expect = 6e-53, Frame = +2 / +2
     *
     * and sets the contents of the HashMap of key value pairs viz. e.g.
     *
     * key: bitScore       -     value: 210
     * key: expectValue    -     value: 6e-53
     *
     * @param poLine     A String representation of the HSP Summary.
     * @param HashMap    To put key/value info in.
     * @param oVersion a <code>BlastLikeVersionSupport</code> value
     * @exception SAXException if an error occurs
     */
    public static void parse(String poLine, HashMap poMap,
              BlastLikeVersionSupport poVersion)
    throws SAXException {

    int iProgram = poVersion.getProgram();

    if ( (iProgram == BlastLikeVersionSupport.NCBI_BLASTN) ||
         (iProgram == BlastLikeVersionSupport.NCBI_BLASTX) ||
         (iProgram == BlastLikeVersionSupport.NCBI_BLASTP) ||
         (iProgram == BlastLikeVersionSupport.NCBI_TBLASTN) ||
         (iProgram == BlastLikeVersionSupport.NCBI_TBLASTX) ) {

        parseNCBIBlast(poLine,poMap,poVersion);
        return;
    }

    if ( (iProgram == BlastLikeVersionSupport.WU_BLASTN) ||
         (iProgram == BlastLikeVersionSupport.WU_BLASTX) ||
         (iProgram == BlastLikeVersionSupport.WU_BLASTP) ||
         (iProgram == BlastLikeVersionSupport.WU_TBLASTN) ||
         (iProgram == BlastLikeVersionSupport.WU_TBLASTX) ) {

        parseWUBlast(poLine,poMap,poVersion);
        return;
    }

    if ( (iProgram == BlastLikeVersionSupport.GCG_BLASTN)) {
        //Similar enough to use NCBI parser for now
        parseNCBIBlast(poLine,poMap,poVersion);
        return;
    }
    //If get here, then program is not supported.
    throw (new SAXException(
        "Failed attempting to parse an HSP Summary because program ".
        concat(poVersion.getProgramString()).
            concat(" is not supported.")));
    }

    static void parseNCBIBlast(String poLine, HashMap poMap,
                   BlastLikeVersionSupport poVersion) {

    String oToken;
    String oToken2;
    String oKey;
    String oValue;
    char[]  aoTmpArray;
    StringTokenizer oSt;
    StringTokenizer oSt2;
    StringTokenizer oTmp;
    StringBuffer    oTmpBuffer = new StringBuffer();
    poMap.clear();

    //System.out.println(">>>>" + poLine);

    //Tokenize on commas, and make lower case...
    oSt = new StringTokenizer(poLine.toLowerCase(),",");

    while (oSt.hasMoreTokens()) {

        oToken = oSt.nextToken().trim();

        oSt2 = new StringTokenizer(oToken);

        while (oSt2.hasMoreTokens()) {

        oToken2 = oSt2.nextToken().trim();

        //now grab info on a case-by-case basis
        //and put into HashMap...

        //NCBI-BLAST, WU-BLAST
        if (oToken2.equals("score")) {

            oKey = "score";
            //assume "Token = value ..."
            oSt2.nextToken(); //skip =
            oValue = oSt2.nextToken(); //grab score

            poMap.put(oKey,oValue);
            break;
        }

        //NCBI Blast, WU-BLAST
        if (oToken2.startsWith("expect")) {
            //could be "expect" or "expect(2) etc."
            oKey = "expectValue";
            //assume " Token = value"
            oSt2.nextToken(); //skip =
            oValue = oSt2.nextToken();

            poMap.put(oKey,oValue);
            break;
        }

        //NCBI Blast, WU-BLAST

        if (oToken2.equals("identities")) {

            //assume " identities = 129/168 (76%)"
            oSt2.nextToken(); //skip =

            //this token is 129/157
            oTmp = new StringTokenizer(oSt2.nextToken(),"/");
            oKey = "numberOfIdentities";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);

            oKey = "alignmentSize";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);

            //here next token is (76%)

            oTmp = new StringTokenizer(oSt2.nextToken(),"(%)");
            oKey = "percentageIdentity";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);

            break;
        }

        //NCBI Blast, WU-BLAST
        if (oToken2.equals("positives")) {

            //assume " positives = 129/168 (76%)"
            oSt2.nextToken(); //skip =

            //this token is 129/157
            oTmp = new StringTokenizer(oSt2.nextToken(),"/");
            oKey = "numberOfPositives";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);


            //here next token is like (76%)

            oTmp = new StringTokenizer(oSt2.nextToken(),"(%)");
            oKey = "percentagePositives";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);


            break;
        }


        //NCBI Blast, WU-BLAST
        if (oToken2.equals("strand")) {

            //assume " strand = plus / minus"
            oSt2.nextToken(); //skip =

            //this token is "plus"

            oKey = "queryStrand";
            oValue = oSt2.nextToken();
            poMap.put(oKey,oValue);

            oSt2.nextToken(); //skip "/"

            oKey = "hitStrand";
            oValue = oSt2.nextToken();
            poMap.put(oKey,oValue);

            break;
        }



        if (oToken2.equals("frame")) {

            //assume " Frame = +3 " for blastx and tblastn
            //assume " Frame = +3 / -1" for tblastx

            oSt2.nextToken(); //skip =

            if (poVersion.getProgram() ==
            BlastLikeVersionSupport.NCBI_BLASTX) {
            oKey = "queryFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);
            break;
            }


            if (poVersion.getProgram() ==
            BlastLikeVersionSupport.NCBI_TBLASTN) {
            oKey = "hitFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);
            break;
            }

            if (poVersion.getProgram() ==
            BlastLikeVersionSupport.NCBI_TBLASTX) {
            oKey = "queryFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);

            //skip "/"

            oSt2.nextToken();

            oKey = "hitFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);

            break;
            }


        }




        } //end loop over

        //System.out.println(oToken);
    } //end loop over "score = 119 bits" type tokens

    }


    static void parseWUBlast(String poLine, HashMap poMap,
                   BlastLikeVersionSupport poVersion) {

    String oToken;
    String oToken2;
    String oKey;
    String oValue;
    char[]  aoTmpArray;
    StringTokenizer oSt;
    StringTokenizer oSt2;
    StringTokenizer oTmp;
    StringBuffer    oTmpBuffer = new StringBuffer();
    poMap.clear();

    //System.out.println(">>>>" + poLine);

    //Tokenize on commas, and make lower case...
    oSt = new StringTokenizer(poLine.toLowerCase(),",");

    while (oSt.hasMoreTokens()) {

        oToken = oSt.nextToken().trim();

        oSt2 = new StringTokenizer(oToken);

        while (oSt2.hasMoreTokens()) {

        oToken2 = oSt2.nextToken().trim();

        //now grab info on a case-by-case basis
        //and put into HashMap...

        //NCBI-BLAST, WU-BLAST
        if (oToken2.equals("score")) {

            oKey = "score";
            //assume "Token = value ..."
            oSt2.nextToken(); //skip =
            oValue = oSt2.nextToken(); //grab score

            poMap.put(oKey,oValue);
            break;
        }

        //NCBI Blast, WU-BLAST
        if (oToken2.startsWith("expect")) {
            //could be "expect" or "expect(2) etc."
            oKey = "expectValue";
            //assume " Token = value"
            oSt2.nextToken(); //skip =
            oValue = oSt2.nextToken();

            poMap.put(oKey,oValue);
            break;
        }

        //NCBI Blast, WU-BLAST

        if (oToken2.equals("identities")) {

            //assume " identities = 129/168 (76%)"
            oSt2.nextToken(); //skip =

            //this token is 129/157
            oTmp = new StringTokenizer(oSt2.nextToken(),"/");
            oKey = "numberOfIdentities";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);

            oKey = "alignmentSize";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);

            //here next token is (76%)

            oTmp = new StringTokenizer(oSt2.nextToken(),"(%)");
            oKey = "percentageIdentity";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);

            break;
        }

        //NCBI Blast, WU-BLAST
        if (oToken2.equals("positives")) {

            //assume " positives = 129/168 (76%)"
            oSt2.nextToken(); //skip =

            //this token is 129/157
            oTmp = new StringTokenizer(oSt2.nextToken(),"/");
            oKey = "numberOfPositives";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);


            //here next token is like (76%)

            oTmp = new StringTokenizer(oSt2.nextToken(),"(%)");
            oKey = "percentagePositives";
            oValue = oTmp.nextToken();
            poMap.put(oKey,oValue);


            break;
        }


        //NCBI Blast, WU-BLAST
        if (oToken2.equals("strand")) {

            //assume " strand = plus / minus"
            oSt2.nextToken(); //skip =

            //this token is "plus"

            oKey = "queryStrand";
            oValue = oSt2.nextToken();
            poMap.put(oKey,oValue);

            oSt2.nextToken(); //skip "/"

            oKey = "hitStrand";
            oValue = oSt2.nextToken();
            poMap.put(oKey,oValue);

            break;
        }



        if (oToken2.equals("frame")) {

            //assume " Frame = +3 " for blastx and tblastn
            //assume " Frame = +3 / -1" for tblastx

            oSt2.nextToken(); //skip =

            if (poVersion.getProgram() ==
            BlastLikeVersionSupport.WU_BLASTX) {
            oKey = "queryFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);
            break;
            }


            if (poVersion.getProgram() ==
            BlastLikeVersionSupport.WU_TBLASTN) {
            oKey = "hitFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);
            break;
            }

            if (poVersion.getProgram() ==
            BlastLikeVersionSupport.WU_TBLASTX) {
            oKey = "queryFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);

            //skip "/"

            oSt2.nextToken();

            oKey = "hitFrame";
            aoTmpArray = oSt2.nextToken().toCharArray();
            oTmpBuffer.setLength(0);
            if (aoTmpArray[0] == '+') {
                oTmpBuffer.append("plus");
            } else {
                oTmpBuffer.append("minus");
            }
            oTmpBuffer.append(aoTmpArray[1]);
            oValue = oTmpBuffer.substring(0);
            poMap.put(oKey,oValue);

            break;
            }


        }




        } //end loop over

        //System.out.println(oToken);
    } //end loop over "score = 119 bits" type tokens

    }



}

