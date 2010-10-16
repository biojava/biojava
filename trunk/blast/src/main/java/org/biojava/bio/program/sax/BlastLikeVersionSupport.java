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

import java.util.StringTokenizer;

/**
 * A Helper class for checking Blast-like program version support.
 * 
 * When the parsing mode is STRICT,
 * parsing will not be attempted on output from
 * unsupported versions of programs.
 * When the parsing mode is LAZY, parsing will be
 * attempted if the program is recognized, but
 * the particular version of the program is
 * not supported. In this case, incorrect results
 * may be obtained.
 *
 * Primary author -
 *                 Simon Brocklehurst (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Colin Hardman      (CAT)
 *                 Stuart Johnston    (CAT)
 *                 Mathieu Wiepert    (Mayo Foundation)
 *                 Keith James        (Sanger Institute)
 *                 Travis Banks
 *
 * Copyright 2000 Cambridge Antibody Technology Group plc.
 *
 * This code released to the biojava project, May 2000
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Travis Banks
 * @version 0.1
 *
 */
final class BlastLikeVersionSupport {

    public static final int    UNKNOWN           = 0;

    public static final int    NCBI_BLASTN       = 1;
    public static final int    NCBI_BLASTP       = 2;
    public static final int    NCBI_BLASTX       = 3;
    public static final int    NCBI_TBLASTN      = 4;
    public static final int    NCBI_TBLASTX      = 5;

    public static final int    WU_BLASTN         = 11;
    public static final int    WU_BLASTP         = 12;
    public static final int    WU_BLASTX         = 13;
    public static final int    WU_TBLASTN        = 14;
    public static final int    WU_TBLASTX        = 15;

    public static final int    HMMER             = 21;

    public static final int    GCG               = 31; 
    public static final int    GCG_BLASTN        = 32; 

    // NCBI Blast
    public static final int    V2_0_11           = 100;
    public static final int    V2_2_2            = 101;
    public static final int    V2_2_3            = 102;
    public static final int	   V2_2_11			 = 103;
    public static final int    V2_2_15           = 104;
    public static final int    V2_2_18           = 105;
    // GCG Blast
    public static final int    V2_0_10           = 150;
    public static final int    V2_1_2            = 151;

    // WU Blast
    public static final int    V2_0A19MP_WASHU   = 200;

    // HMMER
    public static final int    V2_0              = 300;
    public static final int    V2_1_1            = 301;

    public static final int STRICT = 0;
    public static final int LAZY   = 1;

    //set default parsing mode
    private static   int    iMode            = BlastLikeVersionSupport.STRICT;
    private static   int    iProgram         = BlastLikeVersionSupport.UNKNOWN;
    private static   int    iVersion         = BlastLikeVersionSupport.UNKNOWN;
    private static   String oProgramString   = "unknown";
    private static   String oVersionString   = "unknown";
    private static   String oProgramStub;

    public BlastLikeVersionSupport() {
    }
    /**
     * If parsing mode is strict, then check program and version
     * must be supported.
     * If parsing mode is lazy, then unsupported versions of
     * supported programs will return true.
     *
     * Given a program name, and a version checks if the software
     * is supported. Returns true if it is, false if not.
     *
     * @param poProgram  A String representation of the program name
     * @param poVersion  A String representation of the version
     * @return boolean   -
     */
     public boolean isSupported() {
         
         //Check version support for NCBI Blast
         if ( (iProgram == NCBI_BLASTN) ||
              (iProgram == NCBI_BLASTX) ||
              (iProgram == NCBI_BLASTP) ||
              (iProgram == NCBI_TBLASTN) ||
              (iProgram == NCBI_TBLASTX) ) {

             if (iVersion == V2_0_11 || iVersion == V2_2_2 || iVersion == V2_2_3 
                     || iVersion == V2_2_11 || iVersion == V2_2_15 || iVersion == V2_2_18) {
                 return true;
             }

             //if get here, program version is unsupported
             //return false if mode is strict, true if LAZY

             if (this.getMode() == BlastLikeVersionSupport.STRICT) {
                 return false;
             }
             if (this.getMode() == BlastLikeVersionSupport.LAZY) {
                 return true;
             }
         } //end if NCBI Blast

         //Check version support for WU_BLAST
         if ( (iProgram == WU_BLASTN) ||
              (iProgram == WU_BLASTX) ||
              (iProgram == WU_BLASTP) ||
              (iProgram == WU_TBLASTN) ||
              (iProgram == WU_TBLASTX) ) {

             if (iVersion == V2_0A19MP_WASHU) {
                 return true;
             }
         }

         //Check version support for HMMER
         if (iProgram == HMMER) {
             if (iVersion == V2_0 || iVersion == V2_1_1 ) {
                 return true;
             }
         }

         //Check version support for GCG
         if (iProgram == GCG_BLASTN) {
             if ((iVersion == V2_0_10)||
                 (iVersion == V2_1_2)){
                 return true;
             }

             //if get here, program version is unsupported
             //return false if mode is strict, true if LAZY

             if (this.getMode() == BlastLikeVersionSupport.STRICT) {
                 return false; 
             }
             if (this.getMode() == BlastLikeVersionSupport.LAZY) {
                 return true; 
             }
         } //end if GCG Blast

         //if get here, the program version is unsupported
         //return false if mode is strict, true if LAZY

         if (this.getMode() == BlastLikeVersionSupport.STRICT) {
             return false; 
         }
         if (this.getMode() == BlastLikeVersionSupport.LAZY) {
             return true; 
         }

         //if get here, program is unsupported because
         //program type is not recognized
         return false;
     }

    /**
     * Describe 'isStartOfDataSet' method here.
     * @param poLine line to process     
     * @return boolean -
     */
    public boolean isStartOfDataSet(String poLine) {

        if ( (poLine.startsWith("BLAST")) ||
             (poLine.startsWith("HMMER")) || 
             (poLine.startsWith("TBLAST")) ||
             (poLine.startsWith("!!")) ) {  //GCG check  
            if (poLine.startsWith("!!")) {
                //Flag to indicate program type, not needed for
                //other programs
                iProgram = GCG;
            }

            return true;
        }           

        //if get here, not the start of a new dataset
        return false;
    }

    public int getProgram() {
        return iProgram;
    }

    public int getVersion() {
        return iVersion;
    }

    public String getProgramString() {
        return oProgramString;
    }

    public String getVersionString() {
        return oVersionString;
    }

    /**
     * Assign program and version from by parsing a line
     * from the raw output.
     * @param poLine line to process
     * @return true if format recognised (could be wrong versioN), false
     * if the format not recognised at all.
     */
    public boolean assignProgramAndVersion(String poLine) {

        //Take first two tokens only (this is OK for current programs)
        boolean tFormatFound = false;

        StringTokenizer oSt = new StringTokenizer(poLine);

        //deal with Blast, WU-blast, GCG and HMMER

        //first potentially identify program e.g. blastn
        oProgramStub   = oSt.nextToken().toLowerCase();
        oVersionString = oSt.nextToken().toLowerCase();

        //if it's a blast, then choose ncbi-blast, wu-blast or GCG
        if (oProgramStub.indexOf("blast") != -1) {
            if ((oVersionString.indexOf("washu") == -1)) {
                //here if GCG or NCBI-BLAST
                if (iProgram==GCG) { 
                    oProgramString = "gcg-".concat(oProgramStub);
                } else {
                    oProgramString = "ncbi-".concat(oProgramStub);
                }
            } else {
                //here if WU-BLAST
                oProgramString = "wu-".concat(oProgramStub);
            }
        }

        //if it's hmmer
        if (oProgramStub.indexOf("hmmer") != -1) {
            oProgramString = oProgramStub;
        }

        //NCBI blast
        if (oProgramString.equals("ncbi-blastn")) {
            iProgram = NCBI_BLASTN;
            tFormatFound = true;
        }
        if (oProgramString.equals("ncbi-blastx")) {
            iProgram = NCBI_BLASTX;
            tFormatFound = true;
        }
        if (oProgramString.equals("ncbi-blastp")) {
            iProgram = NCBI_BLASTP;
            tFormatFound = true;
        }
        if (oProgramString.equals("ncbi-tblastn")) {
            iProgram = NCBI_TBLASTN;
            tFormatFound = true;
        }
        if (oProgramString.equals("ncbi-tblastx")) {
            iProgram = NCBI_TBLASTX;
            tFormatFound = true;
        }
        if (oVersionString.equals("2.0.11")) {
            iVersion = V2_0_11;
        }
        if (oVersionString.equals("2.2.2")) {
            iVersion = V2_2_2;
        }
        if (oVersionString.equals("2.2.3")) {
            iVersion = V2_2_3;
        }
        if (oVersionString.equals("2.2.11")) {
            iVersion = V2_2_11;
        }
        if (oVersionString.equals("2.2.15")) {
            iVersion = V2_2_15;
        }
        if (oVersionString.equals("2.2.18")) {
            iVersion = V2_2_18;
        }

        //wu-blast
        if (oProgramString.equals("wu-blastn")) {
            iProgram = WU_BLASTN;
            tFormatFound = true;
        }
        if (oProgramString.equals("wu-blastx")) {
            iProgram = WU_BLASTX;
            tFormatFound = true;
        }
        if (oProgramString.equals("wu-blastp")) {
            iProgram = WU_BLASTP;
            tFormatFound = true;
        }
        if (oProgramString.equals("wu-tblastn")) {
            iProgram = WU_TBLASTN;
            tFormatFound = true;
        }
        if (oProgramString.equals("wu-tblastx")) {
            iProgram = WU_TBLASTX;
            tFormatFound = true;
        }
        if (oVersionString.equals("2.0a19mp-washu")) {
            iVersion = V2_0A19MP_WASHU;
            tFormatFound = true;
        }

        //hmmer
        if (oProgramString.equals("hmmer")) {
            iProgram = HMMER;
            tFormatFound = true;
        }
        if (oVersionString.equals("2.0")) {
            iVersion = V2_0;
            tFormatFound = true;
        }
        if (oVersionString.equals("2.1.1")) {
            iVersion = V2_0;
            tFormatFound = true;
        }

        //GCG
        if (oProgramString.equals("gcg-blastn")) {
            iProgram = GCG_BLASTN;  //First pass assigned this to GCG
            //Second pass is more specific.
            tFormatFound = true;
        }

        //GCG Version
        if (oVersionString.equals("2.0.10")) {
            iVersion = V2_0_10;
        }
        //GCG Version
        if (oVersionString.equals("2.1.2")) {
            iVersion = V2_1_2;
        }

        if (!tFormatFound) {
            return false;
        }
        return true;
    }

    /**
     * Set the parsing mode to STRICT or LAZY.
     *
     * @param piMode     Should be one of ParsingMode.STRICT
     * or ParsingMode.LAZY.
     */
    public void setMode(int piMode) {
        iMode = piMode;
    }

    /**
     * Get parsing mode. Typically will be compared to
     * ParsingMode.STRICT or ParsingMode.LAZY to
     * see if parsing should continue or not.
     *
     * @return int   The current parsing mode.
     */
    public int getMode() {
        return iMode;
    }
}

