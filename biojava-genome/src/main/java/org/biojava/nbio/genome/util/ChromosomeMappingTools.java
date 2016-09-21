package org.biojava.nbio.genome.util;

import com.google.common.collect.Range;

import org.biojava.nbio.genome.parsers.genename.ChromPos;
import org.biojava.nbio.genome.parsers.genename.GeneChromosomePosition;


import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

/**
 *  A class that can map chromosomal positions to mRNA (coding sequence) positions.
 *
 *  @author Andreas Prlic
 */

public class ChromosomeMappingTools {

    public static final boolean debug = false;
    private static final String newline = System.getProperty("line.separator");

    public static String formatExonStructure(GeneChromosomePosition chromosomePosition ){
        if ( chromosomePosition.getOrientation() == '+')
            return formatExonStructureForward(chromosomePosition);

        return formatExonStructureReverse(chromosomePosition);

    }

    public static String printHTMLExonStructure(GeneChromosomePosition chromosomePosition ){
        if ( chromosomePosition.getOrientation() == '+')
            return printHTMLExonStructureForward(chromosomePosition);

        return printHTMLExonStructureReverse(chromosomePosition);

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


            //s.append("forward exon: " + start + " - " + end + " | " + (end - start));
            //s.append(newline);

            if (start <= cdsStart +1 && end >= cdsStart+1) {

                inCoding = true;
                codingLength += (end - cdsStart);
                s.append("     UTR         : " + format(start) + " - " + format(cdsStart ));
                s.append(newline);
                s.append(" ->  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + (end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));
                s.append(newline);

            } else if (start <= cdsEnd && end >= cdsEnd) {
                //System.out.println(" <-- CDS end at: " + cdsEnd );
                inCoding = false;
                codingLength += (cdsEnd - start);

                s.append(" <-  Exon        : " + format(start+1) + " - " + format(cdsEnd) + " | " + (cdsEnd - start) + " | " + codingLength + " | " + (codingLength % 3));
                s.append(newline);
                s.append("     UTR         : " + (cdsEnd +1) + " - " + format(end));
                s.append(newline);


            } else if (inCoding) {
                // full exon is coding
                codingLength += (end - start);

                s.append("     Exon        : " + format(start+1) + " - " + format(end) + " | " + (end - start) + " | " + codingLength + " | " + (codingLength % 3));
                s.append(newline);
            }
            //			if ( inCoding )
            //				System.out.println("exon phase at end:" + (codingLength % 3));
            //
            //			System.out.println("   coding length: " + codingLength);


        }
        s.append("Coding Length: ");
        s.append((codingLength-3)+"");
        s.append(newline);
        return s.toString();
    }


    private static String printHTMLExonStructureForward(GeneChromosomePosition chromPos) {
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


            //s.append("forward exon: " + start + " - " + end + " | " + (end - start));
            //s.append(newline);

            if (start <= cdsStart && end >= cdsStart) {

                inCoding = true;
                codingLength += (end - cdsStart);
                s.append("<tr><th>Region</th><th>start</th><th>end</th><th>region length</th><th>phase at end</th></tr>");
                s.append("<tr><td>UTR</td><td>"  + format(start)      + "</td><td>" +  format(cdsStart) + "</td><td></td><td></td></tr>");
                s.append(newline);
                s.append("<tr><td>Exon</td><td>" +  showGenePosLink(chromPos, (cdsStart + 1)) + "</td><td>" +  showGenePosLink(chromPos, end)     + "</td><td>" + (end - cdsStart) + "</td><td>" + (codingLength % 3)+"</td></tr>");
                s.append(newline);

            } else if (start <= cdsEnd && end >= cdsEnd) {
                //System.out.println(" <-- CDS end at: " + cdsEnd );
                inCoding = false;
                codingLength += (cdsEnd - start);

                s.append("<tr><td>Exon</td><td>" +  showGenePosLink(chromPos, (start + 1)) + "</td><td>" +  showGenePosLink(chromPos,cdsEnd) + "</td><td>" + (cdsEnd - start) +  "</td><td>" + (codingLength % 3)+"</td></tr>");
                s.append(newline);
                s.append("<tr><td>UTR</td><td>" +  format(cdsEnd +1) + "</td><td>" +  format(end)+"</td><td></td><td></td></tr>");
                s.append(newline);


            } else if (inCoding) {
                // full exon is coding
                codingLength += (end - start);

                s.append("<tr><td>Exon</td><td>" +  showGenePosLink(chromPos, (start + 1)) + "</td><td>" +  showGenePosLink(chromPos, end) + "</td><td>" + (end - start) + "</td><td>" +  (codingLength % 3)+"</td></tr>");
                s.append(newline);
            }
            //			if ( inCoding )
            //				System.out.println("exon phase at end:" + (codingLength % 3));
            //
            //			System.out.println("   coding length: " + codingLength);


        }
        //s.append("Length of coding sequence: ");
        //s.append((codingLength-3)+" nucleotides.");
        s.append(newline);
        return s.toString();
    }

    private static String showGenePosLink(GeneChromosomePosition chromPos, Integer pos ) {

        String spos = format(pos);

        StringBuffer buf = new StringBuffer();
        buf.append("<a href=\"/pdb/gene/");
        buf.append(chromPos.getGeneName());
        buf.append("?chromosome=");
        buf.append(chromPos.getChromosome());
        buf.append("&range=");
        buf.append(pos);
        buf.append("\">");
        buf.append(spos);
        buf.append("</a>");

        return buf.toString();
    }

    private static String printHTMLExonStructureReverse(GeneChromosomePosition chromPos) {
        StringWriter s = new StringWriter();

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds   = chromPos.getExonEnds();


        int cdsStart = chromPos.getCdsStart();
        int cdsEnd   = chromPos.getCdsEnd();

        // System.out.println("CDS START:" +format(cdsStart) + " - " + format(cdsEnd));

        boolean inCoding = false;
        int codingLength = 0;

        if (cdsEnd < cdsStart) {
            int tmp = cdsEnd;
            cdsEnd = cdsStart;
            cdsStart = tmp;
        }
//        System.out.println("CDS START:" +format(cdsStart) + " - " + format(cdsEnd));

        int lengthExons = 0;
        s.append("<tr><th>Region</th><th>start</th><th>end</th><th>region length</th><th>phase at end</th></tr>");
        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {

            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            lengthExons += end - start;
            //s.append("Reverse exon: " + end + " - " + start + " | " + (end - start));
            //s.append(newline);

            if (start <= cdsEnd && end >= cdsEnd) {
                inCoding = true;


                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }
                codingLength += (cdsEnd - tmpstart);

                s.append("<tr><td><div data-toggle=\"tooltip\" data-placement=\"top\" title=\"Untranslated Region\">UTR</div></td><td>" + format(cdsEnd + 1) + "</td><td>" + format(end)+"</td><td></td><td></td></tr>");
                s.append(newline);

                s.append("<tr><td>Exon</td><td>" + showGenePosLink(chromPos,(tmpstart+1)) + "</td><td>" + showGenePosLink(chromPos, cdsEnd) + "</td><td>" + (cdsEnd - tmpstart) + "</td><td>"  + (codingLength % 3)+"</td></tr>");
                s.append(newline);
                // single exon with UTR on both ends
                if (tmpstart != start)
                    s.append("<tr><td><div data-toggle=\"tooltip\" data-placement=\"top\" title=\"Untranslated Region\">UTR</div></td><td>" + format(cdsStart) + "</td><td>" + format(start + 1) + "</td><td></td><td></td></tr>");
                s.append(newline);

            } else if (start <= cdsStart && end >= cdsStart) {
                inCoding = false;
                codingLength += (end - cdsStart);

                s.append("<tr><td>Exon</td><td>" + showGenePosLink(chromPos,(cdsStart+1)) + "</td><td>" + showGenePosLink(chromPos, end) + "</td><td>" + (end - cdsStart) +"</td><td>"+ (codingLength % 3)+"</td></tr>");
                s.append(newline);
                s.append("<tr><td><div data-toggle=\"tooltip\" data-placement=\"top\" title=\"Untranslated Region\">UTR</div></td><td>" + format(start+1) + "</td><td>" + format(cdsStart)+"</td><td></td><td></td></tr>");
                s.append(newline);


            } else if (inCoding) {
                // full exon is coding
                codingLength += (end - start);

                s.append("<tr><td>Exon</td><td>" + showGenePosLink(chromPos,(start+1)) + "</td><td>" + showGenePosLink(chromPos, end) + "</td><td>" + (end - start) + "</td><td>"  + (codingLength % 3)+"</td></tr>");
                s.append(newline);
            } else {
                // e.g. see UBQLN3
                s.append("<tr><td data-toggle=\"tooltip\" data-placement=\"right\" title=\"Untranslated Region\">UTR</td><td>" + format(start) + "</td><td>" + format(end)+"</td><td></td><td></td></tr>");
                s.append(newline);
            }
        }

        // s.append("Length coding sequence: " + (codingLength-3) + " nucleotides");
        s.append(newline);

        return s.toString();
    }


    private static String formatExonStructureReverse(GeneChromosomePosition chromPos) {
        StringWriter s = new StringWriter();

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds   = chromPos.getExonEnds();


        int cdsStart = chromPos.getCdsStart();
        int cdsEnd   = chromPos.getCdsEnd();

        // System.out.println("CDS START:" +format(cdsStart) + " - " + format(cdsEnd));

        boolean inCoding = false;
        int codingLength = 0;

        if (cdsEnd < cdsStart) {
            int tmp = cdsEnd;
            cdsEnd = cdsStart;
            cdsStart = tmp;
        }
//        System.out.println("CDS START:" +format(cdsStart) + " - " + format(cdsEnd));

        int lengthExons = 0;

        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {

            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            lengthExons += end - start;
            //s.append("Reverse exon: " + end + " - " + start + " | " + (end - start));
            //s.append(newline);

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
                s.append("Exon        :" + format(tmpstart+1) + " - " + format(cdsEnd) + " | " + (cdsEnd - tmpstart) + " | " + codingLength + " | " + (codingLength % 3));
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
                s.append(" no translation! UTR: " + format(start) + " - " + format(end));
                s.append(newline);
            }
        }

        s.append("CDS length: " + (codingLength-3));
        s.append(newline);

        return s.toString();
    }

    /**
     * Get the length of the CDS in base pairs
     *
     * @param chromPos
     * @return
     */
    public static int getCDSLength(GeneChromosomePosition chromPos) {

        if (debug) {

            System.out.println(chromPos);

            System.out.println("chromosomal information: ");
            //
            System.out.println("Gene:" + chromPos.getGeneName());
            System.out.println("  Transcription (including UTRs): " + chromPos.getTranscriptionStart() + " - " + chromPos.getTranscriptionEnd() + " (length:" + (chromPos.getTranscriptionEnd() - chromPos.getTranscriptionStart()) + ")");
            System.out.println("  Orientation: " + chromPos.getOrientation());
            System.out.println("  CDS: " + (chromPos.getCdsStart()) + " - " + chromPos.getCdsEnd() + " (length: " + (chromPos.getCdsEnd() - chromPos.getCdsStart()) + ")");
        }
        //System.out.println("======");
        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds = chromPos.getExonEnds();

        if (debug)
            System.out.println("Exons:" + exonStarts.size());

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
     * maps the position of a CDS nucleotide back to the genome
     *
     * @param cdsNucleotidePosition
     * @return
     */
    public static ChromPos getChromosomePosForCDScoordinate(int cdsNucleotidePosition, GeneChromosomePosition chromPos) {
        if (debug)
            System.out.println(" ? Checking chromosome position for CDS position " + cdsNucleotidePosition);

        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds = chromPos.getExonEnds();
        if (debug)
            System.out.println(" Exons:" + exonStarts.size());

        int cdsStart = chromPos.getCdsStart();
        int cdsEnd = chromPos.getCdsEnd();


        ChromPos chromosomePos = null;

        if (chromPos.getOrientation().equals('+'))

            chromosomePos = ChromosomeMappingTools.getChromPosForward(cdsNucleotidePosition, exonStarts, exonEnds, cdsStart, cdsEnd);
        else
            chromosomePos = ChromosomeMappingTools.getChromPosReverse(cdsNucleotidePosition, exonStarts, exonEnds, cdsStart, cdsEnd);
        if (debug)
            System.out.println("=> CDS pos " + cdsNucleotidePosition + " for " + chromPos.getGeneName() + " is on chromosome at  " + chromosomePos);
        return chromosomePos;

    }

    private static String format(int chromosomePosition){
        // returns a nicely formatted representation of the position

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
    public static ChromPos getChromPosReverse(int cdsPos, List<Integer> exonStarts,
                                              List<Integer> exonEnds, int cdsStart, int cdsEnd) {

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
            if (debug)
                System.out.println("Exon #" + (i+1) + "/" + exonStarts.size());
            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            lengthExons += end - start;
            if (debug) {
                System.out.println("     is " + cdsPos + " part of Reverse exon? " + format(start+1) + " - " + format(end) + " | " + (end - start+1));
                System.out.println("     CDS start: " + format(cdsStart+1) + "-" + format(cdsEnd) + " coding length counter:" + codingLength);
            }
            if (start+1 <= cdsEnd && end >= cdsEnd) {


                // FIRST EXON
                inCoding = true;


                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }

                // here one of the few places where we don't say start+1
                int check = codingLength + cdsEnd - tmpstart ;
                if ( debug) {
                    System.out.println("First Exon    | " + (check) + " | " + format(start+1) + " " + format(end) + " | " + (cdsEnd - tmpstart) + " | " + cdsPos );
                }

                if ( ( check > cdsPos)  )
                {

                    int tmp = cdsPos - codingLength ;

                    if ( debug) {
                        System.out.println(" -> found position in UTR exon:  " + format(cdsPos) + " " + format(tmpstart+1) + " tmp:" + format(tmp) + " cs:" + format(cdsStart+1) + " ce:" + format(cdsEnd) + " cl:" + codingLength);
                    }
                    return new ChromPos((cdsEnd - tmp), -1) ;
                }


                // don't add 1 here
                codingLength += (cdsEnd - tmpstart );
                if (debug) {
                    System.out.println("     UTR         :" + format(cdsEnd + 1) + " - " + format(end));
                    if (tmpstart == start)
                        System.out.print(" ->  ");
                    else
                        System.out.print(" <-> ");
                    System.out.println("Exon        :" + format(tmpstart+1) + " - " + (cdsEnd) + " | " + format(cdsEnd - tmpstart+1) + " - " + codingLength + " | " + (codingLength % 3));

                    // single exon with UTR on both ends
                    if (tmpstart != start)
                        System.out.println("     UTR         :" + format(cdsStart ) + " - " + format(start+1));
                }
            } else if (start <= cdsStart && end >= cdsStart) {

                // LAST EXON
                inCoding = false;

                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }


                if (codingLength + end - cdsStart >= cdsPos) {

                    // how many remaining coding nucleotides?
                    int tmp = codingLength + end - cdsStart - cdsPos ;

                    if ( debug) {

                        System.out.println("cdl: " +codingLength + " tmp:" + tmp + " cdsStart: " + format(cdsStart));

                        System.out.println(" -> XXX found position noncoding exon:  cdsPos:" + cdsPos + " s:" + format(start + 1) + " tmp:" + format(tmp) + " cdsStart:" + (cdsStart + 1) + " codingLength:" + codingLength + " cdsEnd:" + format(cdsEnd));
                    }
                    return new ChromPos((cdsStart + tmp),-1);
                }

                codingLength += (end - cdsStart);
                if (debug) {
                    System.out.println(" <-  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + format(end - cdsStart+1) + " | " + codingLength + " | " + (codingLength % 3));
                    System.out.println("     UTR         : " + format(start+1) + " - " + format(cdsStart ));
                }
            } else if (inCoding) {


                if (codingLength + end - start -1  >= cdsPos) {

                    int tmp = cdsPos - codingLength ;

                    if ( tmp > (end - start ) ) {

                        tmp = (end - start );
                        if ( debug ) {
                            System.out.println("changing tmp to " + tmp);
                        }
                    }
                    if ( debug) {
                        System.out.println("     " + cdsPos + " " + codingLength + " | " + (cdsPos - codingLength) + " | " + (end -start) + " | " + tmp);
                        System.out.println("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start) + " | " + codingLength + " | " + (codingLength % 3));
                        System.out.println(" ->  RRR found position coding exon:  " + cdsPos + " " + format(start+1) + " " + format(end) + " " + tmp + " " + format(cdsStart+1) + " " + codingLength);
                    }
                    return new ChromPos((end - tmp),cdsPos %3);
                }

                // full exon is coding
                codingLength += (end - start) ;
                if (debug)
                    System.out.println("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start+1) + " | " + codingLength + " | " + (codingLength % 3));
            } else {
                // e.g. see UBQLN3
                if (debug)
                    System.out.println(" no translation!");
            }

            //if ( inCoding )
            //	System.out.println("     exon phase at end:" + ((codingLength) % 3));
            if (debug)
                System.out.println("     coding length: " + codingLength + "(phase:" + (codingLength % 3) + ") CDS POS trying to map:" + cdsPos);


        }
        if (debug)
            System.out.println("length exons: " + lengthExons);
        // could not map, or map over the full length??
        return new ChromPos(-1,-1);

        //return codingLength - 3;
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
    public static ChromPos getChromPosForward(int cdsPos, List<Integer> exonStarts, List<Integer> exonEnds,
                                              int cdsStart, int cdsEnd) {
        boolean inCoding = false;
        int codingLength = 0;

        int lengthExons = 0;
        // map forward
        for (int i = 0; i < exonStarts.size(); i++) {


            // start can include UTR
            int start = exonStarts.get(i);
            int end = exonEnds.get(i);

            lengthExons += end - start;

            //if (debug)
            //System.out.println("forward exon: " + start + " - " + end + " | " + (end - start));

            if (start <= cdsStart +1 && end >= cdsStart+1) {
                // first exon with UTR
                if (codingLength + (end - cdsStart-1) >= cdsPos) {
                    // we are reaching our target position
                    int tmp = cdsPos - codingLength;

                    if ( debug) {
                        System.out.println(cdsStart + " | " + codingLength + " | " + tmp);
                        System.out.println(" -> found position in UTR exon:  #"+(i+1)+ " cdsPos:" + cdsPos +
                                " return:"+(cdsStart +1 + tmp) +" start:" + format(start + 1) + " " + format(tmp) + " " + cdsStart + " " + codingLength);
                    }
                    // we start 1 after cdsStart...
                    return new ChromPos((cdsStart +1 + tmp),-1);
                }
                inCoding = true;
                codingLength += (end - cdsStart);
                if (debug) {
                    System.out.println("     UTR         : " + format(start+1) + " - " + (cdsStart ));
                    System.out.println(" ->  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + format(end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));
                }
            } else if (start+1 <= cdsEnd && end >= cdsEnd) {
                // LAST EXON with UTR
                //System.out.println(" <-- CDS end at: " + cdsEnd );
                inCoding = false;
                if (codingLength + (cdsEnd - start-1) >= cdsPos) {
                    int tmp = cdsPos - codingLength;

                    if ( debug) {
                        System.out.println(" <-  Exon        : " + format(start+1) + " - " + format(cdsEnd) + " | " + format(cdsEnd - start) + " | " + codingLength + " | " + (codingLength % 3));
                        System.out.println("     UTR         : " + format(cdsEnd + 1) + " - " + format(end));
                        System.out.println( codingLength + " | " + tmp + " | " + format(start+1));
                        System.out.println(" -> chromPosForward found position in non coding exon:  " + cdsPos + " " + format(start+1) + " " + format(tmp) + " " + format(cdsStart) + " " + codingLength);
                    }
                    return new ChromPos((start +1 + tmp),cdsPos%3);
                }
                codingLength += (cdsEnd - start-1);
                if (debug) {
                    System.out.println(" <-  Exon        : " + format(start+1) + " - " + format(cdsEnd) + " | " + format(cdsEnd - start) + " | " + codingLength + " | " + (codingLength % 3));
                    System.out.println("     UTR         : " + format(cdsEnd + 1) + " - " + format(end));
                }

            } else if (inCoding) {
                // A standard coding Exon
                // tests for the maximum length of this coding exon
                if (codingLength + (end - start -1)  >= cdsPos) {

                    // we are within the range of this exon
                    int tmp = cdsPos - codingLength ;

                    if ( debug) {
                        System.out.println("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start) + " | " + tmp + " | " + codingLength);
                        System.out.println(" -> found chr position in coding exon #" + (i+1) + ":  cdsPos:" + format(cdsPos) + " s:" + format(start) + "-" + format(end) + " tmp:" + format(tmp) + " cdsStart:" + format(cdsStart) + " codingLength:" + codingLength);
                    }
                    return new ChromPos((start +1 + tmp),cdsPos%3);
                }
                // full exon is coding
                codingLength += (end - start );
                if (debug)
                    System.out.println("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start) + " | " + codingLength + " | " + (codingLength % 3));
            }
            //			if ( inCoding )
            //				System.out.println("exon phase at end:" + (codingLength % 3));
            //
            //			System.out.println("   coding length: " + codingLength);


        }

        //System.out.println("length exons: " + lengthExons);
        //return codingLength - 3;

        // could not map!

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
    public static int getCDSLengthReverse(List<Integer> exonStarts,
                                          List<Integer> exonEnds, int cdsStart, int cdsEnd) {

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

            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            lengthExons += end - start;

            if (start <= cdsEnd && end >= cdsEnd) {
                inCoding = true;


                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }
                codingLength += (cdsEnd - tmpstart);
                if (debug) {
                    System.out.println("     UTR         :" + (cdsEnd + 1) + " - " + (end));
                    if (tmpstart == start)
                        System.out.print(" ->  ");
                    else
                        System.out.print(" <-> ");
                    System.out.println("Exon        :" + tmpstart + " - " + cdsEnd + " | " + (cdsEnd - tmpstart) + " | " + codingLength + " | " + (codingLength % 3));
                    // single exon with UTR on both ends
                    if (tmpstart != start)
                        System.out.println("     UTR         :" + (cdsStart - 1) + " - " + start);
                }
            } else if (start <= cdsStart && end >= cdsStart) {
                inCoding = false;
                codingLength += (end - cdsStart);
                if (debug) {
                    System.out.println(" <-  Exon        : " + cdsStart + " - " + end + " | " + (end - cdsStart) + " | " + (codingLength-3)  + " | " + (codingLength % 3));
                    System.out.println("     UTR         : " + start + " - " + (cdsStart - 1));
                }

            } else if (inCoding) {
                // full exon is coding
                codingLength += (end - start);
                if (debug)
                    System.out.println("     Exon        : " + start + " - " + end + " | " + (end - start) + " | " + codingLength + " | " + (codingLength % 3));
            } else {
                // e.g. see UBQLN3
                if (debug)
                    System.out.println(" no translation!");
            }

            //if ( inCoding )
            //	System.out.println("     exon phase at end:" + ((codingLength) % 3));

            //System.out.println("   coding length: " + codingLength + "(phase:"+(codingLength %3)+")");


        }
        if (debug)
            System.out.println("length exons: " + lengthExons + " codin length: " + (codingLength - 3));
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
    public static int getCDSLengthForward(List<Integer> exonStarts, List<Integer> exonEnds,
                                          int cdsStart, int cdsEnd) {
        boolean inCoding = false;
        int codingLength = 0;

        int lengthExons = 0;
        // map forward
        for (int i = 0; i < exonStarts.size(); i++) {

            int start = exonStarts.get(i);
            int end = exonEnds.get(i);
            lengthExons += end - start;

            if (debug)
                System.out.println("forward exon: " + (start+1) + " - " + end + " | " + (end - start));

            if (start+1 <= cdsStart +1 && end >= cdsStart+1) {

                inCoding = true;
                codingLength += (end - cdsStart);
                if (debug) {
                    System.out.println("     UTR         : " + start + " - " + (cdsStart ));
                    System.out.println(" ->  Exon        : " + (cdsStart+1) + " - " + end + " | " + (end - cdsStart+1) + " | " + codingLength + " | " + (codingLength % 3));
                }
            } else if (start+1 <= cdsEnd && end >= cdsEnd) {
                //System.out.println(" <-- CDS end at: " + cdsEnd );
                inCoding = false;
                codingLength += (cdsEnd - start);
                if (debug) {
                    System.out.println(" <-  Exon        : " + (start +1)+ " - " + cdsEnd + " | " + (cdsEnd - start+1) + " | " + codingLength + " | " + (codingLength % 3));
                    System.out.println("     UTR         : " + cdsEnd + 1 + " - " + end);
                }

            } else if (inCoding) {
                // full exon is coding
                codingLength += (end - start);
                if (debug)
                    System.out.println("     Exon        :" + (start+1) + " - " + end + " | " + (end - start+1) + " | " + codingLength + " | " + (codingLength % 3));
            }
            //			if ( inCoding )
            //				System.out.println("exon phase at end:" + (codingLength % 3));
            //
            //			System.out.println("   coding length: " + codingLength);


        }
        if (debug) {
            System.out.println("length exons: " + lengthExons);
            System.out.println("CDS length:" + (codingLength-3));
        }
        return codingLength-3 ;
    }



    /** Extracts the exon boundaries in CDS coordinates. (needs to be divided by 3 to get AA positions)
     *
     * @param chromPos
     * @return
     */
    public static List<Range> getCDSExonRanges(GeneChromosomePosition chromPos){
        if ( chromPos.getOrientation() == '+')
            return getCDSExonRangesForward(chromPos);

        return getCDSExonRangesReverse(chromPos);
    }

    private static List<Range> getCDSExonRangesReverse(GeneChromosomePosition chromPos) {
        List<Integer> exonStarts = chromPos.getExonStarts();
        List<Integer> exonEnds   = chromPos.getExonEnds();

        List<Range> data = new ArrayList<Range>();
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
        if ( debug)
            s = new StringBuffer();

        int lengthExons = 0;

        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {

            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            lengthExons += end - start;
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
                    s.append("     UTR         :" + format(cdsEnd + 1) + " | " + format(end));
                    s.append(newline);
                    if (tmpstart == start)
                        s.append(" ->  ");
                    else
                        s.append(" <-> ");
                    s.append("Exon        :" + format(tmpstart + 1) + " - " + format(cdsEnd) + " | " + (cdsEnd - tmpstart) + " | " + codingLength + " | " + (codingLength % 3));
                    s.append(newline);
                    // single exon with UTR on both ends
                    if (tmpstart != start)
                        s.append("     UTR         :" + format(cdsStart) + " - " + format(start + 1));
                    s.append(newline);
                }
                Range r = Range.closed(0,codingLength);
                data.add(r);

            } else if (start <= cdsStart && end >= cdsStart) {
                inCoding = false;

                Range r = Range.closed(codingLength,codingLength+(end-cdsStart));
                data.add(r);

                codingLength += (end - cdsStart);
                if  (debug) {
                    s.append(" <-  Exon        : " + format(cdsStart + 1) + " - " + format(end) + " | " + (end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));
                    s.append(newline);
                    s.append("     UTR         : " + format(start + 1) + " - " + format(cdsStart));
                    s.append(newline);
                }


            } else if (inCoding) {
                // full exon is coding

                Range r = Range.closed(codingLength,codingLength+(end-start));
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
            s.append("CDS length: " + (codingLength - 3));
            s.append(newline);
            System.out.println(s.toString());
        }

        return data;
    }


    private static List<Range> getCDSExonRangesForward(GeneChromosomePosition chromPos) {

        List<Range> data = new ArrayList<Range>();
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
//
                Range r = Range.closed(0,codingLength);
                data.add(r);

            } else if (start <= cdsEnd && end >= cdsEnd) {
                //System.out.println(" <-- CDS end at: " + cdsEnd );
                inCoding = false;

                Range r = Range.closed(codingLength,codingLength+(cdsEnd-start));
                data.add(r);
                codingLength += (cdsEnd - start);

            } else if (inCoding) {
                // full exon is coding

                Range r = Range.closed(codingLength,codingLength+(end-start));
                data.add(r);
                codingLength += (end - start);

            }

        }

        return data;
    }
//

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


    public static int getCDSPosReverse(int chromPos, List<Integer> exonStarts, List<Integer> exonEnds,
                                       int cdsStart, int cdsEnd) {
        boolean inCoding = false;
        int codingLength = 0;

        if (cdsEnd < cdsStart) {
            int tmp = cdsEnd;
            cdsEnd = cdsStart;
            cdsStart = tmp;
        }

        if ( debug )
            System.out.println("looking for CDS position for " +format(chromPos));


        if ( chromPos <  cdsStart+1 ) {
            // this is not in a coding region!
            if (debug)
                System.out.println(chromPos + " < " + cdsStart+1 );
            return -1;
        }

        if ( chromPos >  cdsEnd+1 ) {
            // this is not in a coding region!
            if (debug)
                System.out.println(chromPos + " > " + cdsEnd+1 );
            return -1;
        }

        int lengthExons = 0;

        // map reverse
        for (int i = exonStarts.size() - 1; i >= 0; i--) {
            if (debug)
                System.out.println("Reverse Exon #" + (i+1) + "/" + exonStarts.size());
            int end = exonStarts.get(i);
            int start = exonEnds.get(i);

            if (end < start) {
                int tmp = end;
                end = start;
                start = tmp;
            }
            lengthExons += end - start;

            if (debug) {
                System.out.println("     is " + format(chromPos) + " part of Reverse exon? s:" + format(start+1) + " - e:" + format(end) + " | " + (end - start+1));
                System.out.println("     CDS start: " + format(cdsStart+1) + "-" + format(cdsEnd) + " coding length counter:" + codingLength);

            }

            if (start+1 <= cdsEnd && end >= cdsEnd ) {

                // first exon with UTR

                inCoding = true;

                int tmpstart = start;
                if (start < cdsStart) {
                    tmpstart = cdsStart;
                }

                if ( debug) {
                    System.out.println(" --- codingLength " + codingLength +
                            " s:" +
                            format(tmpstart+1) +
                            " e:" +
                            format(cdsEnd) +
                            " p:" +
                            format(chromPos) + " tmp: " + (chromPos - cdsStart));

                    System.out.println("check: " + (codingLength + cdsEnd - tmpstart+1) + " ==?? " + format(chromPos));
                }
                int tmp = cdsEnd - chromPos ;
                // if (codingLength + cdsEnd - tmpstart >= chromPos) {
                //if (end >= chromPos && start + (end-start) >= chromPos) {
                // if (codingLength + cdsEnd - tmpstart >= chromPos) {
                if ( chromPos >= start +1 && chromPos <= end){

                    if ( debug)
                        System.out.println(" -> found position in UTR exon:  P: " + format(chromPos) + " s:" + format(tmpstart+1) + " l:" + format(tmp) + " cdsS:" + format(cdsStart+1) + " cdsE:" + format(cdsEnd)  + " codingL:" + codingLength);
                    return codingLength + tmp;
                }

                if ( debug)
                    System.out.println("     codinglength " + codingLength + " + " + (cdsEnd - tmpstart ) );

                // do not add 1 here
                codingLength += (cdsEnd - tmpstart );
                if (debug) {
                    System.out.println("     UTR         :" + format(cdsEnd + 1) + " - " + format(end));
                    if (tmpstart == start)
                        System.out.print(" ->  ");
                    else
                        System.out.print(" <-> ");
                    System.out.println("Reverse Exon        :" + format(tmpstart+1) + " - " + (cdsEnd) + " | " + format(cdsEnd - tmpstart) + " - " + codingLength + " | " + (codingLength % 3));

                    // single exon with UTR on both ends
                    if (tmpstart != start)
                        System.out.println("     UTR         :" + format(cdsStart - 1) + " - " + format(start));
                }
            } else if (start <= cdsStart && end >= cdsStart) {

                // terminal exon with UTR
                inCoding = false;

                if ( debug)
                    System.out.println(format(start  + codingLength + end - cdsStart) + " ?? " + format(chromPos));
                // (start  + codingLength + end - cdsStart >= chromPos &&
                if (( start+1 <= chromPos) && ( end >= chromPos)) {

                    //int tmp =  end - cdsStart ;
//                    int tmp =  chromPos - cdsStart ;
//                    int l = end - cdsStart;
                    int tmp = end-chromPos ;
                    if  ( tmp > end -cdsStart) {
                        tmp = end-cdsStart ;
                        if ( debug)
                            System.out.println("Adjust tmp to " + tmp);
                    }

                    if ( debug) {

                        System.out.println(  codingLength + " | " + (end -chromPos) + " | " + (end - cdsStart) );
                        System.out.println(" <-  Exon        : " + format(cdsStart) + " - " + format(end) + " | " + format(end - cdsStart +1) + " | ");
                        System.out.println("     UTR         : " + format(start+1) + " - " + format(cdsStart ));
                        System.out.println(" <- YYY found position noncoding exon:  #" + (i+1) + " "  + format(chromPos) + " s:" + format(start) + " tmp: " + format(tmp) + " cdsStart" + format(cdsStart) + " cdl:" + codingLength + " " + format(cdsEnd));
                    }
                    return codingLength + tmp;
                }

                if ( debug )
                    System.out.println("     codinglength " + codingLength + " + " + (end - cdsStart) );
                codingLength += (end - cdsStart+1);
                if (debug) {
                    System.out.println(" <-  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + format(end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));
                    System.out.println("     UTR         : " + format(start+1) + " - " + format(cdsStart ));
                }
            } else if (inCoding) {
                // standard coding exon
                // if (codingLength + end - start >= chromPos) {
                if ( chromPos >= start+1 && chromPos <= end) {

                    int tmp = end -chromPos  ;
                    if ( tmp > (end-start+1)) {

                        tmp = (end - start+1);
                        if ( debug)
                            System.out.println("Adjusting tmp to " + tmp);
                    }
                    if ( debug) {

                        System.out.println(" -> found position in reverse coding exon:  #" + (i+1) + " chromPos:"  + format(chromPos) + " start:" + format(start+1) + " end:" + format(end) + " tmp:" + tmp + " cdsStart:" + cdsStart + " codingLength:" + codingLength);
                    }
                    return codingLength+tmp;
                }

                // full exon is coding
                if ( debug )
                    System.out.println("     codinglength " + codingLength + " + " + (end - start) );
                // don't add 1
                codingLength += (end - start);
                if (debug)
                    System.out.println("     Exon        : " + format(start+1) + " - " + format(end) + " | " + format(end - start) + " | " + codingLength + " | " + (codingLength % 3));
            } else {
                // e.g. see UBQLN3
                if (debug) {
                    System.out.println(" no translation! cdl:" + codingLength);
                }
            }

            //if ( inCoding )
            //	System.out.println("     exon phase at end:" + ((codingLength) % 3));
            if (debug)
                System.out.println("     coding length: " + codingLength + "(phase:" + (codingLength % 3) + ") CDS POS trying to map:" + chromPos);


        }
        if (debug)
            System.out.println("length exons: " + lengthExons);
        // could not map, or map over the full length??


        return -1;

    }

    /**
     * Get the chromosome position mapped onto the mRNA CDS transcript position (needs to be divided by 3 to get protein coordinate)
     *
     * @param exonStarts
     * @param exonEnds
     * @param cdsStart
     * @param cdsEnd
     * @return
     */
    public static int getCDSPosForward(int chromPos, List<Integer> exonStarts, List<Integer> exonEnds,
                                       int cdsStart, int cdsEnd) {
        boolean inCoding = false;
        int codingLength = 0;

        if ( debug)
            System.out.println("looking for CDS position for " +chromPos);

        int lengthExons = 0;
        // map forward
        for (int i = 0; i < exonStarts.size(); i++) {


            // start can include UTR
            int start = exonStarts.get(i);
            int end = exonEnds.get(i);

            lengthExons += end - start;

            if (debug)
                System.out.println("forward exon: " + (start+1) + " - " + end + " | " + (end - start) + " ? overlaps with " + format(chromPos));

            if (start +1 <= cdsStart +1 && end >= cdsStart+1) {

                if (end >= chromPos) {
                    // we are reaching our target position
                    // -1 is important  here...
                    int tmp = chromPos - cdsStart -1;

                    if ( debug) {
                        System.out.println("cdl:" + codingLength + " | " + tmp);
                        System.out.println(" -> found position in UTR exon:  " + chromPos + " " + format(start) + " " + format(tmp) + " " + cdsStart + " " + codingLength);
                    }
                    return codingLength + tmp;
                }
                inCoding = true;
                codingLength += (end - cdsStart);
                if (debug) {
                    System.out.println("     UTR         : " + format(start) + " - " + (cdsStart ));
                    System.out.println(" ->  Exon        : " + format(cdsStart+1) + " - " + format(end) + " | " + format(end - cdsStart) + " | " + codingLength + " | " + (codingLength % 3));
                }
            } else if (start <= cdsEnd && end >= cdsEnd) {
                //System.out.println(" <-- CDS end at: " + cdsEnd );
                inCoding = false;
                if (cdsEnd >= chromPos && (start +1 <= chromPos)) {
                    int tmp =  chromPos - start -1 ;

                    if ( debug)
                        System.out.println(" -> cdsForward found position in non coding exon#"+i+":  " + chromPos + " " + format(start+1) + " " + format(tmp) + " " + cdsStart   + " " + codingLength);
                    return codingLength + tmp ;
                }
                codingLength += (cdsEnd - start);
                if (debug) {
                    System.out.println(" <-  Exon        : " + format(start+1) + " - " + format(cdsEnd) + " | " + format(cdsEnd - start+1) + " | " + codingLength + " | " + (codingLength % 3));
                    System.out.println("     UTR         : " + format(cdsEnd + 1) + " - " + format(end));
                }

            } else if (inCoding) {

                if (end >= chromPos && (start +1 <=chromPos)) {

                    int tmp = chromPos-start-1 ;

                    if ( debug) {
                        System.out.println(codingLength + " | " + tmp);
                        System.out.println(" -> found position in coding exon #" + (i + 1) + ":  " + format(chromPos) + " " + format(start + 1) + " " + format(tmp) + " " + cdsStart + " " + codingLength);
                    }

                    return codingLength + tmp ;
                }
                // full exon is coding
                codingLength += (end - start);
                if (debug)
                    System.out.println("     Exon        :" + format(start) + " - " + format(end) + " | " + format(end - start) + " | " + codingLength + " | " + (codingLength % 3));
            }
            //			if ( inCoding )
            //				System.out.println("exon phase at end:" + (codingLength % 3));
            //
            //			System.out.println("   coding length: " + codingLength);


        }

        //System.out.println("length exons: " + lengthExons);
        //return codingLength - 3;

        // could not map!

        return -1;
    }


}
