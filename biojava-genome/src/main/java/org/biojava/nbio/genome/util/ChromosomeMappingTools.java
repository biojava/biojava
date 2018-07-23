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
package org.biojava.nbio.genome.util;

import com.google.common.collect.Range;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.biojava.nbio.genome.parsers.genename.ChromPos;
import org.biojava.nbio.genome.parsers.genename.GeneChromosomePosition;
import org.biojava.nbio.genome.parsers.twobit.TwoBitFacade;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

/**
 *  A class that can map chromosomal positions to mRNA (coding sequence) positions.
 *
 *  @author Andreas Prlic
 */

public class ChromosomeMappingTools {

    private static final Logger logger = LoggerFactory.getLogger(ChromosomeMappingTools.class);

    private static final String newline = System.getProperty("line.separator");

    public static final String CHROMOSOME = "CHROMOSOME";
    public static final String CDS = "CDS";

    private static int base = 1;
    public static void setCoordinateSystem(int baseInt) {
        base = baseInt;
    }

    /** 
     * Pretty print the details of a GeneChromosomePosition to a String
     *
     * @param chromosomePosition
     * @return
     */
    public static String formatExonStructure(GeneChromosomePosition chromosomePosition ){
        if ( chromosomePosition.getOrientation() == '+')
            return formatExonStructureForward(chromosomePosition);
        return formatExonStructureReverse(chromosomePosition);
    }

    private static String formatExonStructureForward(GeneChromosomePosition chromPos) {

        StringWriter s = new StringWriter();

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds   = chromPos.getExonEnds();

        int cdsStart = chromPos.getCdsStart();
        int cdsEnd   = chromPos.getCdsEnd();

        boolean inCoding = false;
        int codingLength = 0;

        for (int i = 0; i < exonStarts.size(); i++) {

            int start = exonStarts.get(i);
            int end = exonEnds.get(i);

            if (start <= cdsStart +1 && end >= cdsStart+1) {

                inCoding = true;
                codingLength += (end - cdsStart);
                s.append("     UTR         : ").append(format(start)).append(" - ").append(format(cdsStart));
                s.append(newline);
                s.append(" ->  Exon        : ").append(format(cdsStart + 1)).append(" - ").append(format(end)).append(" | ").append(Integer.toString(end - cdsStart)).append(" | ").append(Integer.toString(codingLength)).append(" | ").append(Integer.toString(codingLength % 3));
                s.append(newline);

            } else if (start <= cdsEnd && end >= cdsEnd) {
                //logger.debug(" <-- CDS end at: " + cdsEnd );
                inCoding = false;
                codingLength += (cdsEnd - start);

                s.append(" <-  Exon        : ").append(format(start + 1)).append(" - ").append(format(cdsEnd)).append(" | ").append(Integer.toString(cdsEnd - start)).append(" | ").append(Integer.toString(codingLength)).append(" | ").append(Integer.toString(codingLength % 3));
                s.append(newline);
                s.append("     UTR         : " + (cdsEnd +1) + " - " + format(end));
                s.append(newline);

            } else if (inCoding) {
                // full exon is coding
                codingLength += (end - start);

                s.append("     Exon        : ").append(format(start + 1)).append(" - ").append(format(end)).append(" | ").append(Integer.toString(end - start)).append(" | ").append(Integer.toString(codingLength)).append(" | ").append(Integer.toString(codingLength % 3));
                s.append(newline);
            }
        }
        s.append("Coding Length: ");
        s.append((codingLength-3)+"");
        s.append(newline);
        return s.toString();
    }

    private static String formatExonStructureReverse(GeneChromosomePosition chromPos) {
        StringWriter s = new StringWriter();

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds   = chromPos.getExonEnds();


        int cdsStart = chromPos.getCdsStart();
        int cdsEnd   = chromPos.getCdsEnd();

        // logger.debug("CDS START:" +format(cdsStart) + " - " + format(cdsEnd));

        boolean inCoding = false;
        int codingLength = 0;

        if (cdsEnd < cdsStart) {
            int tmp = cdsEnd;
            cdsEnd = cdsStart;
            cdsStart = tmp;
        }

        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {

            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }

            if (start <= cdsEnd && end >= cdsEnd) {
                inCoding = true;

                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }
                codingLength += (cdsEnd - tmpstart);

                s.append("     UTR         :" + format(cdsEnd + 1) + " | " + format(end));
                s.append(newline);
                if (tmpstart == start)
                    s.append(" ->  ");
                else
                    s.append(" <-> ");
                s.append("Exon        :").append(format(tmpstart + 1)).append(" - ").append(format(cdsEnd)).append(" | ").append(Integer.toString(cdsEnd - tmpstart)).append(" | ").append(Integer.toString(codingLength)).append(" | ").append(Integer.toString(codingLength % 3));
                s.append(newline);
                // single exon with UTR on both ends
                if (tmpstart != start)
                    s.append("     UTR         :" + format(cdsStart ) + " - " + format(start + 1));
                s.append(newline);

            } else if (start <= cdsStart && end >= cdsStart) {
                inCoding = false;
                codingLength += (end - cdsStart);

                s.append(" <-  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + (end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));
                s.append(newline);
                s.append("     UTR         : " + format(start+1) + " - " + format(cdsStart ));
                s.append(newline);


            } else if (inCoding) {
                // full exon is coding
                codingLength += (end - start);

                s.append("     Exon        : " + format(start+1) + " - " + format(end) + " | " + (end - start) + " | " + codingLength + " | " + (codingLength % 3));
                s.append(newline);
            } else {
                // e.g. see UBQLN3
                s.append(" no translation! UTR: ").append(format(start)).append(" - ").append(format(end));
                s.append(newline);
            }
        }

        s.append("CDS length: ").append(Integer.toString(codingLength - 3));
        s.append(newline);

        return s.toString();
    }

    /**
     * Get the length of the CDS in nucleotides.
     *
     * @param chromPos
     * @return length of the CDS in nucleotides.
     */
    public static int getCDSLength(GeneChromosomePosition chromPos) {

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds = chromPos.getExonEnds();

        int cdsStart = chromPos.getCdsStart();
        int cdsEnd = chromPos.getCdsEnd();

        int codingLength;
        if (chromPos.getOrientation().equals('+'))
            codingLength = ChromosomeMappingTools.getCDSLengthForward(exonStarts, exonEnds, cdsStart, cdsEnd);
        else
            codingLength = ChromosomeMappingTools.getCDSLengthReverse(exonStarts, exonEnds, cdsStart, cdsEnd);
        return codingLength;
    }

    /**
     * Maps the position of a CDS nucleotide back to the genome
     *
     * @param cdsNucleotidePosition
     * @return a ChromPos object
     */
    public static ChromPos getChromosomePosForCDScoordinate(int cdsNucleotidePosition, GeneChromosomePosition chromPos) {

        logger.debug(" ? Checking chromosome position for CDS position " + cdsNucleotidePosition);

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds = chromPos.getExonEnds();

        logger.debug(" Exons:" + exonStarts.size());

        int cdsStart = chromPos.getCdsStart();
        int cdsEnd = chromPos.getCdsEnd();


        ChromPos chromosomePos = null;

        if (chromPos.getOrientation().equals('+'))

            chromosomePos = ChromosomeMappingTools.getChromPosForward(cdsNucleotidePosition, exonStarts, exonEnds, cdsStart, cdsEnd);
        else
            chromosomePos = ChromosomeMappingTools.getChromPosReverse(cdsNucleotidePosition, exonStarts, exonEnds, cdsStart, cdsEnd);

        logger.debug("=> CDS pos " + cdsNucleotidePosition + " for " + chromPos.getGeneName() + " is on chromosome at  " + chromosomePos);
        return chromosomePos;

    }

    /** 
     * Returns a nicely formatted representation of the position
     *
     * @param chromosomePosition
     * @return
     */
    private static String format(int chromosomePosition){
        return String.format("%,d", chromosomePosition);
    }

    /**
     * Get the CDS position mapped on the chromosome position
     *
     * @param exonStarts
     * @param exonEnds
     * @param cdsStart
     * @param cdsEnd
     * @return
     */
    public static ChromPos getChromPosReverse(int cdsPos, List<Integer> exonStarts, List<Integer> exonEnds, int cdsStart, int cdsEnd) {

        boolean inCoding = false;
        int codingLength = 0;

        if (cdsEnd < cdsStart) {
            int tmp = cdsEnd;
            cdsEnd = cdsStart;
            cdsStart = tmp;
        }

        int lengthExons = 0;

        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {

            logger.debug("Exon #" + (i+1) + "/" + exonStarts.size());
            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            lengthExons += end - start;

            logger.debug("     is " + cdsPos + " part of Reverse exon? " + format(start+1) + " - " + format(end) + " | " + (end - start+1));
            logger.debug("     CDS start: " + format(cdsStart+1) + "-" + format(cdsEnd) + " coding length counter:" + codingLength);

            if (start+1 <= cdsEnd && end >= cdsEnd) {

                // FIRST EXON
                inCoding = true;

                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }

                // here one of the few places where we don't say start+1
                int check = codingLength + cdsEnd - tmpstart ;

                logger.debug("First Exon    | " + (check) + " | " + format(start+1) + " " + format(end) + " | " + (cdsEnd - tmpstart) + " | " + cdsPos );


                if ( ( check > cdsPos)  )  {
                    int tmp = cdsPos - codingLength ;
                    logger.debug(" -> found position in UTR exon:  " + format(cdsPos) + " " + format(tmpstart+1) + " tmp:" + format(tmp) + " cs:" + format(cdsStart+1) + " ce:" + format(cdsEnd) + " cl:" + codingLength);
                    return new ChromPos((cdsEnd - tmp), -1) ;
                }

                // don't add 1 here
                codingLength += (cdsEnd - tmpstart );

                boolean debug = logger.isDebugEnabled();

                if ( debug ) {

                    StringBuffer b = new StringBuffer();

                    b.append("     UTR         :" + format(cdsEnd + 1) + " - " + format(end) + newline);
                    if (tmpstart == start)
                        b.append(" ->  ");
                    else
                        b.append(" <-> ");
                    b.append("Exon        :" + format(tmpstart + 1) + " - " + (cdsEnd) + " | " + format(cdsEnd - tmpstart + 1) + " - " + codingLength + " | " + (codingLength % 3) + newline);

                    // single exon with UTR on both ends
                    if (tmpstart != start)
                        b.append("     UTR         :" + format(cdsStart) + " - " + format(start + 1) + newline);

                    logger.debug(b.toString());
                }
            } else if (start <= cdsStart && end >= cdsStart) {

                // LAST EXON
                inCoding = false;

                if (codingLength + end - cdsStart >= cdsPos) {

                    // how many remaining coding nucleotides?
                    int tmp = codingLength + end - cdsStart - cdsPos ;

                    logger.debug("cdl: " +codingLength + " tmp:" + tmp + " cdsStart: " + format(cdsStart));

                    logger.debug(" -> XXX found position noncoding exon:  cdsPos:" + cdsPos + " s:" + format(start + 1) + " tmp:" + format(tmp) + " cdsStart:" + (cdsStart + 1) + " codingLength:" + codingLength + " cdsEnd:" + format(cdsEnd));

                    return new ChromPos((cdsStart + tmp),-1);
                }

                codingLength += (end - cdsStart);

                logger.debug(" <-  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + format(end - cdsStart+1) + " | " + codingLength + " | " + (codingLength % 3));
                logger.debug("     UTR         : " + format(start+1) + " - " + format(cdsStart ));

            } else if (inCoding) {

                if (codingLength + end - start -1  >= cdsPos) {

                    int tmp = cdsPos - codingLength ;

                    if ( tmp > (end - start ) ) {
                        tmp = (end - start );
                        logger.debug("changing tmp to " + tmp);
                    }
                    logger.debug("     " + cdsPos + " " + codingLength + " | " + (cdsPos - codingLength) + " | " + (end -start) + " | " + tmp);
                    logger.debug("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start) + " | " + codingLength + " | " + (codingLength % 3));
                    logger.debug(" ->  RRR found position coding exon:  " + cdsPos + " " + format(start+1) + " " + format(end) + " " + tmp + " " + format(cdsStart+1) + " " + codingLength);

                    return new ChromPos((end - tmp),cdsPos %3);
                }
                // full exon is coding
                codingLength += (end - start) ;

                logger.debug("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start+1) + " | " + codingLength + " | " + (codingLength % 3));
            } else {
                // e.g. see UBQLN3
                logger.debug(" no translation!");
            }
            logger.debug("     coding length: " + codingLength + "(phase:" + (codingLength % 3) + ") CDS POS trying to map:" + cdsPos);
        }

        logger.debug("length exons: " + lengthExons);
        // could not map, or map over the full length??
        return new ChromPos(-1,-1);

    }

    /**
     * Get the CDS position mapped onto the chromosome position
     *
     * @param exonStarts
     * @param exonEnds
     * @param cdsStart
     * @param cdsEnd
     * @return
     */
    public static ChromPos getChromPosForward(int cdsPos, List<Integer> exonStarts, List<Integer> exonEnds, int cdsStart, int cdsEnd) {
        boolean inCoding = false;
        int codingLength = 0;

        @SuppressWarnings("unused")
		int lengthExons = 0;
        // map forward
        for (int i = 0; i < exonStarts.size(); i++) {

            // start can include UTR
            int start = exonStarts.get(i);
            int end = exonEnds.get(i);

            lengthExons += end - start;

            if (start <= cdsStart +1 && end >= cdsStart+1) {
                // first exon with UTR
                if (codingLength + (end - cdsStart-1) >= cdsPos) {
                    // we are reaching our target position
                    int tmp = cdsPos - codingLength;


                    logger.debug(cdsStart + " | " + codingLength + " | " + tmp);
                    logger.debug(" -> found position in UTR exon:  #"+(i+1)+ " cdsPos:" + cdsPos +
                            " return:"+(cdsStart +1 + tmp) +" start:" + format(start + 1) + " " + format(tmp) + " " + cdsStart + " " + codingLength);

                    // we start 1 after cdsStart...
                    return new ChromPos((cdsStart +1 + tmp),-1);
                }
                inCoding = true;
                codingLength += (end - cdsStart);

                logger.debug("     UTR         : " + format(start+1) + " - " + (cdsStart ));
                logger.debug(" ->  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + format(end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));

            } else if (start+1 <= cdsEnd && end >= cdsEnd) {
                // LAST EXON with UTR
                //logger.debug(" <-- CDS end at: " + cdsEnd );
                inCoding = false;
                if (codingLength + (cdsEnd - start-1) >= cdsPos) {
                    int tmp = cdsPos - codingLength;

                    logger.debug(" <-  Exon        : " + format(start+1) + " - " + format(cdsEnd) + " | " + format(cdsEnd - start) + " | " + codingLength + " | " + (codingLength % 3));
                    logger.debug("     UTR         : " + format(cdsEnd + 1) + " - " + format(end));
                    logger.debug( codingLength + " | " + tmp + " | " + format(start+1));
                    logger.debug(" -> chromPosForward found position in non coding exon:  " + cdsPos + " " + format(start+1) + " " + format(tmp) + " " + format(cdsStart) + " " + codingLength);

                    return new ChromPos((start +1 + tmp),cdsPos%3);
                }
                codingLength += (cdsEnd - start-1);

                logger.debug(" <-  Exon        : " + format(start+1) + " - " + format(cdsEnd) + " | " + format(cdsEnd - start) + " | " + codingLength + " | " + (codingLength % 3));
                logger.debug("     UTR         : " + format(cdsEnd + 1) + " - " + format(end));


            } else if (inCoding) {
                // A standard coding Exon
                // tests for the maximum length of this coding exon
                if (codingLength + (end - start -1)  >= cdsPos) {

                    // we are within the range of this exon
                    int tmp = cdsPos - codingLength ;

                    logger.debug("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start) + " | " + tmp + " | " + codingLength);
                    logger.debug(" -> found chr position in coding exon #" + (i+1) + ":  cdsPos:" + format(cdsPos) + " s:" + format(start) + "-" + format(end) + " tmp:" + format(tmp) + " cdsStart:" + format(cdsStart) + " codingLength:" + codingLength);

                    return new ChromPos((start +1 + tmp),cdsPos%3);
                }
                // full exon is coding
                codingLength += (end - start );

                logger.debug("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start) + " | " + codingLength + " | " + (codingLength % 3));
            }
        }
        return new ChromPos(-1,-1);
    }

    /**
     * Get the length of the coding sequence
     *
     * @param exonStarts
     * @param exonEnds
     * @param cdsStart
     * @param cdsEnd
     * @return
     */
    public static int getCDSLengthReverse(List<Integer> exonStarts, List<Integer> exonEnds, int cdsStart, int cdsEnd) {

        int codingLength = 0;

        if (cdsEnd < cdsStart) {
            int tmp = cdsEnd;
            cdsEnd = cdsStart;
            cdsStart = tmp;
        }
        cdsStart = cdsStart + base;

        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {

            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            start = start + base;

            if ((start < cdsStart && end < cdsStart) || (start > cdsEnd && end > cdsEnd))
                continue;

            if (start < cdsStart)
                start = cdsStart;

            if (end > cdsEnd)
                end = cdsEnd;

            codingLength += (end - start + 1);
        }
        return codingLength - 3;
    }

    /**
     * Get the length of the coding sequence
     *
     * @param exonStarts
     * @param exonEnds
     * @param cdsStart
     * @param cdsEnd
     * @return
     */
    public static int getCDSLengthForward(List<Integer> exonStarts, List<Integer> exonEnds, int cdsStart, int cdsEnd) {

        int codingLength = 0;

        for (int i = 0; i < exonStarts.size(); i++) {

            int start = exonStarts.get(i)+base;
            int end = exonEnds.get(i);

            if ( (start < cdsStart+base && end < cdsStart) || (start > cdsEnd && end > cdsEnd) )
                continue;

            if (start < cdsStart+base)
                start = cdsStart+base;

            if (end > cdsEnd)
                end = cdsEnd;

            codingLength += (end - start + 1);
        }
        return codingLength - 3;
    }

    /** 
     * Extracts the exon boundaries in CDS coordinates. (needs to be divided by 3 to get AA positions)
     *
     * @param chromPos
     * @return
     */
    public static List<Range<Integer>> getCDSExonRanges(GeneChromosomePosition chromPos){
        if ( chromPos.getOrientation() == '+')
            return getCDSExonRangesForward(chromPos,CDS);
        return getCDSExonRangesReverse(chromPos,CDS);
    }

    /** Extracts the boundaries of the coding regions in chromosomal coordinates
     *
     * @param chromPos
     * @return
     */
    public static List<Range<Integer>> getChromosomalRangesForCDS(GeneChromosomePosition chromPos){
        if ( chromPos.getOrientation() == '+')
            return getCDSExonRangesForward(chromPos,CHROMOSOME);
        return getCDSExonRangesReverse(chromPos,CHROMOSOME);
    }

    private static List<Range<Integer>> getCDSExonRangesReverse(GeneChromosomePosition chromPos, String responseType) {

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds   = chromPos.getExonEnds();

        List<Range<Integer>> data = new ArrayList<>();
        int cdsStart = chromPos.getCdsStart();
        int cdsEnd   = chromPos.getCdsEnd();

        boolean inCoding = false;
        int codingLength = 0;

        if (cdsEnd < cdsStart) {
            int tmp = cdsEnd;
            cdsEnd = cdsStart;
            cdsStart = tmp;
        }

        java.lang.StringBuffer s =null;

        boolean debug = logger.isDebugEnabled();

        if ( debug)
            s = new StringBuffer();

		//int lengthExons = 0;

        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {

            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            //lengthExons += end - start;
            //s.append("Reverse exon: " + end + " - " + start + " | " + (end - start));
            //s.append(newline);

            if (start <= cdsEnd && end >= cdsEnd) {
                inCoding = true;

                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }
                codingLength += (cdsEnd - tmpstart);
                if ( debug ) {
                    s.append("     UTR         :").append(format(cdsEnd + 1)).append(" | ").append(format(end));
                    s.append(newline);
                    if (tmpstart == start)
                        s.append(" ->  ");
                    else
                        s.append(" <-> ");
                    s.append("Exon        :").append(format(tmpstart + 1)).append(" - ").append(format(cdsEnd)).append(" | ").append(cdsEnd - tmpstart).append(" | ").append(codingLength).append(" | ").append(codingLength % 3);
                    s.append(newline);
                    // single exon with UTR on both ends
                    if (tmpstart != start)
                        s.append("     UTR         :").append(format(cdsStart)).append(" - ").append(format(start + 1));
                    s.append(newline);
                }


                Range<Integer> r ;
                if ( responseType.equals(CDS))
                    r = Range.closed(0,codingLength);
                else
                    r = Range.closed(tmpstart,cdsEnd);

                data.add(r);

            } else if (start <= cdsStart && end >= cdsStart) {
                inCoding = false;

                Range<Integer> r;
                if ( responseType.equals(CDS))
                    r = Range.closed(codingLength,codingLength+(end-cdsStart));
                else
                    r = Range.closed(cdsStart+1,end);

                data.add(r);

                codingLength += (end - cdsStart);
                if  (debug) {
                    s.append(" <-  Exon        : " + format(cdsStart + 1) + " - " + format(end) + " | " + (end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));
                    s.append(newline);
                    s.append("     UTR         : ").append(format(start + 1)).append(" - ").append(format(cdsStart));
                    s.append(newline);
                }
            } else if (inCoding) {
                // full exon is coding
                Range<Integer> r;
                if ( responseType.equals(CDS))
                    r = Range.closed(codingLength,codingLength+(end-start));
                else
                    r = Range.closed(start,end);
                data.add(r);

                codingLength += (end - start);
                if  (debug) {
                    s.append("     Exon        : " + format(start + 1) + " - " + format(end) + " | " + (end - start) + " | " + codingLength + " | " + (codingLength % 3));
                    s.append(newline);
                }
            } else {
                // e.g. see UBQLN3
                if ( debug ) {
                    s.append(" no translation! UTR: " + format(start) + " - " + format(end));
                    s.append(newline);
                }
            }
        }
        if ( debug ) {
            s.append("CDS length: ").append(Integer.toString(codingLength - 3));
            s.append(newline);
            logger.debug(s.toString());
        }

        return data;
    }

    private static List<Range<Integer>> getCDSExonRangesForward(GeneChromosomePosition chromPos, String responseType) {

        List<Range<Integer>> data = new ArrayList<>();
        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds   = chromPos.getExonEnds();

        int cdsStart = chromPos.getCdsStart();
        int cdsEnd   = chromPos.getCdsEnd();

        boolean inCoding = false;
        int codingLength = 0;

        for (int i = 0; i < exonStarts.size(); i++) {

            int start = exonStarts.get(i);
            int end = exonEnds.get(i);

            if (start <= cdsStart && end >= cdsStart) {

                inCoding = true;
                codingLength += (end - cdsStart);

                Range<Integer> r;
                if ( responseType.equals(CDS))
                    r = Range.closed(0,codingLength);
                else
                    r = Range.closed(cdsStart,end);
                data.add(r);

            } else if (start <= cdsEnd && end >= cdsEnd) {
                //logger.debug(" <-- CDS end at: " + cdsEnd );
                inCoding = false;

                Range<Integer> r;
                if ( responseType.equals(CDS))
                    r = Range.closed(codingLength,codingLength+(cdsEnd-start));
                else
                    r = Range.closed(start,cdsEnd);
                data.add(r);
                codingLength += (cdsEnd - start);

            } else if (inCoding) {
                // full exon is coding
                Range<Integer> r;
                if ( responseType.equals(CDS))
                    r = Range.closed(codingLength,codingLength+(end-start));
                else
                    r = Range.closed(start,end);
                data.add(r);
                codingLength += (end - start);
            }
        }
        return data;
    }

    /**
     * I have a genomic coordinate, where is it on the mRNA
     *
     * @param coordinate
     * @param chromosomePosition
     * @return
     */
    public static int getCDSPosForChromosomeCoordinate(int coordinate, GeneChromosomePosition chromosomePosition) {

        if ( chromosomePosition.getOrientation() == '+')
        	return getCDSPosForward(coordinate,
                    chromosomePosition.getExonStarts(),
                    chromosomePosition.getExonEnds(),
                    chromosomePosition.getCdsStart(),
                    chromosomePosition.getCdsEnd());

        return getCDSPosReverse(coordinate,
                chromosomePosition.getExonStarts(),
                chromosomePosition.getExonEnds(),
                chromosomePosition.getCdsStart(),
                chromosomePosition.getCdsEnd());
    }
    
	/** 
	 * Converts the genetic coordinate to the position of the nucleotide on the mRNA sequence for a gene 
	 * living on the forward DNA strand.
	 * 
	 * @param chromPos The genetic coordinate on a chromosome 
     * @param exonStarts The list holding the genetic coordinates pointing to the start positions of the exons (including UTR regions)  
     * @param exonEnds The list holding the genetic coordinates pointing to the end positions of the exons (including UTR regions)
     * @param cdsStart The start position of a coding region
     * @param cdsEnd The end position of a coding region
     * 
     * @return the position of the nucleotide base on the mRNA sequence corresponding to the input genetic coordinate (base 1)
	 * 
	 * @author Yana Valasatava
	 */
    public static int getCDSPosForward(int chromPos, List<Integer> exonStarts, List<Integer> exonEnds,
            int cdsStart, int cdsEnd) {
    	
    	// the genetic coordinate is not in a coding region
        if ( (chromPos < (cdsStart+base) ) || ( chromPos > (cdsEnd+base) ) ) {
        	logger.debug("The "+format(chromPos)+" position is not in a coding region");
            return -1;
        }
        
        logger.debug("looking for CDS position for " +format(chromPos));
        
        // map the genetic coordinates of coding region on a stretch of a reverse strand
        List<Range<Integer>> cdsRegions = getCDSRegions(exonStarts, exonEnds, cdsStart, cdsEnd);
        
        int codingLength = 0;
        int lengthExon = 0;
        for (Range<Integer> range : cdsRegions) {
        	
		    int start = range.lowerEndpoint();
		    int end = range.upperEndpoint();
		    
		    lengthExon = end - start;

		    if (start+base <= chromPos && end >= chromPos ) {
		    	return codingLength + (chromPos-start);
		    }
	        else { 
	        	codingLength += lengthExon;
	        }
        }
        return -1;
    }
    
	/** 
	 * Converts the genetic coordinate to the position of the nucleotide on the mRNA sequence for a gene 
	 * living on the reverse DNA strand.
	 * 
	 * @param chromPos The genetic coordinate on a chromosome 
     * @param exonStarts The list holding the genetic coordinates pointing to the start positions of the exons (including UTR regions)  
     * @param exonEnds The list holding the genetic coordinates pointing to the end positions of the exons (including UTR regions)
     * @param cdsStart The start position of a coding region
     * @param cdsEnd The end position of a coding region
     * 
     * @return the position of the nucleotide base on the mRNA sequence corresponding to the input genetic coordinate (base 1)
	 * 
	 * @author Yana Valasatava
	 */
    public static int getCDSPosReverse(int chromPos, List<Integer> exonStarts, List<Integer> exonEnds,
            int cdsStart, int cdsEnd) {
    	
    	// the genetic coordinate is not in a coding region
        if ( (chromPos < (cdsStart+base)) || ( chromPos > (cdsEnd+base) ) ) {
        	logger.debug("The "+format(chromPos)+" position is not in a coding region");
            return -1;
        }
        
        logger.debug("looking for CDS position for " +format(chromPos));
                
        // map the genetic coordinate on a stretch of a reverse strand
        List<Range<Integer>> cdsRegions = getCDSRegions(exonStarts, exonEnds, cdsStart, cdsEnd);
        
        int codingLength = 0;
        int lengthExon = 0;
        for ( int i=cdsRegions.size()-1; i>=0; i-- ) {
        	
		    int start = cdsRegions.get(i).lowerEndpoint();
		    int end = cdsRegions.get(i).upperEndpoint();
		    
		    lengthExon = end - start;
		    // +1 offset to be a base 1
		    if (start+base <= chromPos && end >= chromPos ) {
		    	return codingLength + (end-chromPos+1);
		    }
	        else { 
	        	codingLength += lengthExon;
	        }
        }
        return -1;
    }
    
    /** 
     * Extracts the exons boundaries in CDS coordinates corresponding to the forward DNA strand.
     *
     * @param origExonStarts The list holding the genetic coordinates pointing to the start positions of the exons (including UTR regions)
     * @param origExonEnds The list holding the genetic coordinates pointing to the end positions of the exons (including UTR regions)
     * @param cdsStart The start position of a coding region
     * @param cdsEnd The end position of a coding region
     * 
     * @return the list of genetic positions corresponding to the exons boundaries in CDS coordinates
     */
    public static List<Range<Integer>> getCDSRegions(List<Integer> origExonStarts, List<Integer> origExonEnds, int cdsStart, int cdsEnd) {
    	
        // remove exons that are fully landed in UTRs
        List<Integer> exonStarts = new ArrayList<Integer>(origExonStarts);
        List<Integer> exonEnds = new ArrayList<Integer>(origExonEnds);
        
        int j=0;
        for (int i = 0; i < origExonStarts.size(); i++) {
        	if ( ( origExonEnds.get(i) < cdsStart) || ( origExonStarts.get(i) > cdsEnd) ) {
        		exonStarts.remove(j);
        		exonEnds.remove(j);
        	}
        	else {
        		j++;
        	}
        }
        
        // remove untranslated regions from exons
        int nExons = exonStarts.size();
        exonStarts.remove(0);
        exonStarts.add(0, cdsStart);
        exonEnds.remove(nExons-1);
        exonEnds.add(cdsEnd);
    	
        List<Range<Integer>> cdsRegion = new ArrayList<Range<Integer>>();
        for ( int i=0; i<nExons; i++ ) {
        	Range<Integer> r = Range.closed(exonStarts.get(i), exonEnds.get(i));
        	cdsRegion.add(r);
        }
		return cdsRegion;
    }
    
    /** 
     * Extracts the DNA sequence transcribed from the input genetic coordinates.
     *
     * @param twoBitFacade the facade that provide an access to a 2bit file
     * @param gcp The container with chromosomal positions
     * 
     * @return the DNA sequence transcribed from the input genetic coordinates
     */
    public static DNASequence getTranscriptDNASequence(TwoBitFacade twoBitFacade, GeneChromosomePosition gcp) throws Exception {
    	return getTranscriptDNASequence(twoBitFacade,gcp.getChromosome(),gcp.getExonStarts(), gcp.getExonEnds(), gcp.getCdsStart(), gcp.getCdsEnd(), gcp.getOrientation());
    }
    
    /** 
     * Extracts the DNA sequence transcribed from the input genetic coordinates.
     *
     * @param chromosome the name of the chromosome
     * @param exonStarts The list holding the genetic coordinates pointing to the start positions of the exons (including UTR regions)  
     * @param exonEnds The list holding the genetic coordinates pointing to the end positions of the exons (including UTR regions)
     * @param cdsStart The start position of a coding region
     * @param cdsEnd The end position of a coding region
     * @param orientation The orientation of the strand where the gene is living
     * 
     * @return the DNA sequence transcribed from the input genetic coordinates
     */
	public static DNASequence getTranscriptDNASequence(TwoBitFacade twoBitFacade, String chromosome, List<Integer> exonStarts, List<Integer> exonEnds, int cdsStart, int cdsEnd, Character orientation) throws Exception {

		List<Range<Integer>> cdsRegion = getCDSRegions(exonStarts, exonEnds, cdsStart, cdsEnd);

		String dnaSequence = "";
		for (Range<Integer> range : cdsRegion) {
			String exonSequence = twoBitFacade.getSequence(chromosome,range.lowerEndpoint(), range.upperEndpoint());
            dnaSequence += exonSequence;
		}
		if (orientation.equals('-')) {
            dnaSequence = new StringBuilder(dnaSequence).reverse().toString();
			DNASequence dna = new DNASequence(dnaSequence);
			SequenceView<NucleotideCompound> compliment = dna.getComplement();
            dnaSequence = compliment.getSequenceAsString();
		}
		return new DNASequence(dnaSequence.toUpperCase());
	}
}
